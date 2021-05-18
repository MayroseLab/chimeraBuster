#!/usr/bin/env python

"""
chimeraBuster is a tool for the detection and correction
of chimeric gene annotations.
This script is the first step of the workflow, aimed
at detecting chimeric genes based on transcript mapping.
For further information please visit the tool Github:
https://github.com/MayroseLab/chimeraBuster
"""

# Built-in/Generic Imports
from argparse import ArgumentParser
import os
import sys
import logging
from copy import deepcopy
from itertools import combinations, chain
from multiprocessing import Pool

# Libs
import gffutils
from intervaltree import IntervalTree, Interval
from pafpy import PafFile
from Bio import SeqIO
from sklearn.cluster import DBSCAN
import pandas as pd

from util import *

__author__ = "Lior Glick"
__license__ = "MIT"
__version__ = "0.2.0"
__maintainer__ = "Lior Glick"
__status__ = "Dev"


# FUNCTIONS

def filter_paf(paf_path, min_qcov=0.95):
    """
    Filter a PAF file according to
    specified criteria. Return a list
    of PAF records
    """
    records = []
    with PafFile(paf_path) as paf:
        for record in paf:
            if record.mlen / record.qlen >= min_qcov:
                records.append(record)
    return records


def paf_records_to_intervals(paf_records):
    """
    Convert a list of PAF records
    into a dict where keys are
    chromosomes and values are
    IntervalTree objects. Coordinates
    are always of the target.
    """
    iv_dict = {}
    for rec in paf_records:
        if rec.tname not in iv_dict:
            iv_dict[rec.tname] = IntervalTree()
        iv_dict[rec.tname].add(Interval(rec.tstart, rec.tend, [rec.qname]))
    return iv_dict

def merge_overlaps_min_frac(self, data_reducer=None, data_initializer=None, min_frac=0.05):
    """
    Same as IntervalTree.merge_overlaps(), except intervals are only
    merged if the fraction of overlap in both intervals exceeds the
    cutoff specified by min_frac.
    """
    if not self:
        return

    sorted_intervals = sorted(self.all_intervals)  # get sorted intervals
    merged = []
    # use mutable object to allow new_series() to modify it
    current_reduced = [None]
    higher = None  # iterating variable, which new_series() needs access to

    def new_series():
        if data_initializer is None:
            current_reduced[0] = higher.data
            merged.append(higher)
            return
        else:  # data_initializer is not None
            current_reduced[0] = copy(data_initializer)
            current_reduced[0] = data_reducer(current_reduced[0], higher.data)
            merged.append(Interval(higher.begin, higher.end, current_reduced[0]))

    for higher in sorted_intervals:
        if merged:  # series already begun
            lower = merged[-1]
            if intervals_frac_overlap(lower,higher) >= min_frac and intervals_frac_overlap(higher,lower) >= min_frac:
                upper_bound = max(lower.end, higher.end)
                if data_reducer is not None:
                    current_reduced[0] = data_reducer(current_reduced[0], higher.data)
                else:  # annihilate the data, since we don't know how to merge it
                    current_reduced[0] = None
                merged[-1] = Interval(lower.begin, upper_bound, current_reduced[0])
            else:
                new_series()
        else:  # not merged; is first of Intervals to merge
            new_series()

    self.__init__(merged)

IntervalTree.merge_overlaps_min_frac = merge_overlaps_min_frac

def find_mapping_regions(iv_dict, overlaps_min_frac=0.05):
    """
    Detect rough regions in which
    transcripts were mapped. These
    regions will be refined later.
    """
    res = deepcopy(iv_dict)
    for chrom in res:
        res[chrom].merge_overlaps_min_frac(min_frac=overlaps_min_frac, data_reducer=lambda iv1, iv2: iv1 + iv2)
    n_mr = sum([len(res[chrom]) for chrom in res])
    logging.info("{} mapping regions detected.".format(n_mr))
    return res


def intervals_frac_overlap(iv1, iv2):
    """
    Given two intervals, calculates the
    fraction of iv1 covered by iv2
    """
    return iv1.overlap_size(iv2.begin, iv2.end) / iv1.length()



def remove_outliers(intervals, dbscan_eps=300, dbscan_minPts=4, max_outlier_frac=0.2):
    """
    Given a set of intervals, detect
    outliers based on start and end
    coordinates, using the DBSCAN algorithm.
    If the fraction of intervals detected as
    outliers exceeds max_outlier_frac, this is
    an unreliable case - return original intervals set.
    Otherwise - return set of non-outlier intervals
    *TODO*: automatically choose epsilon and minPts
    """
    # create DF of interval start-end positions
    df = pd.DataFrame([[iv.begin, iv.end, iv] for iv in intervals])
    df.columns = ['start', 'end', 'interval']
    # cluster intervals using DBSCAN
    clustering = DBSCAN(eps=dbscan_eps, min_samples=dbscan_minPts).fit(df[['start', 'end']])
    df['cluster'] = pd.Series(clustering.labels_).astype(str)
    # intervals with cluster '-1' are outliers - suspected bridges or other
    n_outliers = df.query('cluster == "-1"').shape[0]
    n_iv = len(intervals)
    outlier_frac = n_outliers / n_iv
    logging.debug("{} out of {} intervals (fraction = {}) detected as outliers".format(n_outliers, n_iv, outlier_frac))
    if outlier_frac > max_outlier_frac:
        return intervals
    return set(df.query('cluster != "-1"')['interval'])


def refine_mapping_regions(iv_tree, mapping_regions, min_transcripts=2, ncpu=1):
    """
    Works through rough mapping regions
    and refines each by removing misleading
    intervals
    """
    # check if enough transcripts in MR
    large = lambda mr: len(iv_tree[mr.begin:mr.end]) >= min_transcripts
    mapping_regions_filter = IntervalTree(filter(large, mapping_regions))
    # remove outliers
    transcripts_per_mr = [iv_tree[mr.begin:mr.end] for mr in mapping_regions_filter]
    with Pool(ncpu) as pool:
        transcripts_per_mr_rem_outliers = pool.map(remove_outliers, transcripts_per_mr)
        # transcripts_per_mr_rem_outliers is a list of sets of Intervals
        mapping_regions_rem_outliers = [IntervalTree(mr) for mr in transcripts_per_mr_rem_outliers]
        refined = IntervalTree()
        for mr in mapping_regions_rem_outliers:
            mr.merge_overlaps()
            refined = refined.union(mr)
        return refined

def mapping_regions_to_bed(mr_dict, out_file):
    """
    Write all mapping regions to BED file
    """
    with open(out_file, 'w') as fo:
        for chrom in mr_dict:
            for mr in mr_dict[chrom]:
                print("%s\t%s\t%s" %(chrom, mr.begin, mr.end), file=fo)

def ensure_id(gff_feature):
    """
    If a gff feature object has no ID
    attribute, add one. The ID is taken
    from the Name attribute, if it exists,
    otherwise - generated based on Parent.
    """
    if 'ID' in gff_feature.attributes:
        pass
    elif 'Name' in gff_feature.attributes:
        gff_feature['ID'] = gff_feature['Name']
    elif 'Parent' in gff_feature.attributes:
        gff_feature['ID'] = [gff_feature['Parent'][0] + '_%s' % gff_feature.featuretype]
    else:
        logging.error("Encountered a GFF feature with no ID, Name, or Parent attribute:\n{}".format(str(gff_feature)))
        sys.exit()
    return gff_feature


def detect_overlapping_genes(gff_db, min_overlap=0.1):
    """
    Scans the genome annotation for regions
    containing overlapping genes. If the
    fraction of overlap in both genes is
    lower than min_overlap, the genes are
    not considered overlapping.
    Returns a set of overlapping gene IDs.
    """
    ol_set = set()
    for gene in gff_db.features_of_type('gene', order_by='start'):
        region = (gene.seqid, gene.start, gene.end)
        overlapping_genes = list(gff_db.region(region, featuretype=['gene']))
        if len(overlapping_genes) > 1:
            for gene_pair in combinations(overlapping_genes, 2):
                gene1 = Interval(gene_pair[0].start, gene_pair[0].end)
                gene2 = Interval(gene_pair[1].start, gene_pair[1].end)
                if intervals_frac_overlap(gene1, gene2) >= min_overlap and intervals_frac_overlap(gene2,
                                                                                                  gene1) >= min_overlap:
                    ol_set = ol_set.union(set([ol_gene['ID'][0] for ol_gene in gene_pair]))
    return ol_set


def detect_chimeras(gff_db, mapping_regions, min_intersect_frac=0.5):
    """
    Go through all genes and detect
    putative chimeras.
    min_intersect_frac is the minimal
    fraction of the transcript length
    that overlaps with the gene to be
    considered.
    """
    overlapping_genes = detect_overlapping_genes(gff_db)
    chimeras = {}
    for gene in gff_db.features_of_type('gene'):
        gene_id = gene['ID'][0]
        if gene_id in overlapping_genes:
            continue
        chrom = gene.seqid
        if chrom not in mapping_regions:
            continue
        intersect_mappings = mapping_regions[chrom][gene.start:gene.end]
        gene_iv = Interval(gene.start, gene.end, gene_id)
        intersect_mappings = [im for im in intersect_mappings if
                              intervals_frac_overlap(im, gene_iv) >= min_intersect_frac]
        n_im = len(intersect_mappings)
        logging.debug("Gene {} : {} intersecting mapping regions".format(gene_id, n_im))
        if n_im > 1:
            chimeras[gene] = intersect_mappings
    return chimeras

def print_chimeric_genes(chimeras_dict, ofile):
    """
    Print a dictionary of chimeric genes as:
    gene_id    chrom    start1-end1,start2-end2,...
    """
    with open(ofile, 'w') as fo:
        for cg in chimeras_dict:
            gene_id = cg['ID'][0]
            chrom = cg.seqid
            map_regions = ['%s-%s' %(mr.begin,mr.end) for mr in chimeras_dict[cg]]
            print('\t'.join([gene_id,chrom,','.join(map_regions)]), file=fo)


def main():
    # command line arguments
    DESC = "detect_chimeric_genes.py: detect chimeric gene models based on transcript mapping"
    USAGE = "python detect_chimeric_genes.py -h (display full usage doc)"
    parser = ArgumentParser(description=DESC, usage=USAGE)
    parser.add_argument('-g', '--gff', help='Input GFF', type=lambda x: is_valid_file(parser, x), default=None)
    parser.add_argument('-f', '--genome_fasta', help='Input genome in FASTA format', required=True,
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument('-t', '--transcripts', default=None,
                        help='Transcripts to be mapped to the genome, in FASTA format',
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument('-m', '--mapping', default=None,
                        help='Transcripts mapping to the genome, in minimap2 PAF format',
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument('-d', '--gff_db', default=None, help='GFF DB from a previous run',
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument('-o', '--output', help='Output path', required=True)
    parser.add_argument('-r', '--do_not_refine',
                        help='Skip mapping region refining step for a quick-and-dirty analysis', action='store_true',
                        default=False)
    parser.add_argument('-n', '--min_transcripts', type=int, default=2,
                        help='Minimum number of transcripts in mapping region')
    parser.add_argument('-q', '--min_query_cov', type=float, default=0.95,
                        help='Minimum transcript query coverage (0-1)')
    parser.add_argument('-c', '--cpus', type=int, default=1, help='Number of CPUs to use')
    parser.add_argument('-v', '--verbose', action="store_true", default=False, help='Increase verbosity')
    args = parser.parse_args()

    # logging
    log_format = '(%(processName)-10s) %(asctime)s %(message)s'
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG, format=log_format)
    else:
        logging.basicConfig(level=logging.INFO, format=log_format)

    if (not (args.transcripts or args.mapping)) or (args.transcripts and args.mapping):
        logging.error("Must provide either transcripts FASTA or mapping PAF, but not both")
        sys.exit()

    if (not (args.gff or args.gff_db)) or (args.gff and args.gff_db):
        logging.error("Must provide either a GFF or GFF DB, but not both")
        sys.exit()

    if not os.path.exists(args.output):
        os.mkdir(args.output)

    # Read GFF
    if args.gff:
        logging.info("Loading GFF from file...")
        db_path = os.path.join(args.output, os.path.basename(args.gff) + '.db')
        gff_db = gffutils.create_db(args.gff, dbfn=db_path, keep_order=True, merge_strategy='create_unique',
                                    sort_attribute_values=True, verbose=args.verbose, force=True, transform=ensure_id)
    else:
        logging.info("Loading GFF from DB...")
        db_path = args.gff_db
    gff_db = gffutils.FeatureDB(db_path, keep_order=True)

    # Read genome fasta
    logging.info("Loading genome FASTA...")
    chrom_seq = SeqIO.to_dict(SeqIO.parse(args.genome_fasta, "fasta"))

    # Run minimap, if needed
    if args.transcripts:
        logging.info("Running Minimap2...")
        args.mapping = os.path.join(args.output, os.path.basename(args.transcripts) + '_map.paf')
        try:
            minimap_cmd = "minimap2 -x splice:hq -uf {} {} -t {}".format(args.genome_fasta, args.transcripts, args.cpus)
            run_external_cmd(minimap_cmd, args.mapping)
        except:
            logging.error("Minimap run failed. Terminating")
            sys.exit()
    logging.info("Using PAF file {}".format(args.mapping))

    # Filter PAF file for high quality mappings
    if not os.path.isfile(args.mapping):
        logging.error("PAF file {} does not exist. Terminating.".format(args.mapping))
    logging.info("Filtering PAF records...")
    paf_hq_records = filter_paf(args.mapping)
    logging.debug("{} HQ PAF records".format(len(paf_hq_records)))

    # Find and refine mapping regions
    logging.info("Computing transcript mapping regions...")
    mapping_intervals = paf_records_to_intervals(paf_hq_records)
    mapping_regions = find_mapping_regions(mapping_intervals)
    if not args.do_not_refine:
        logging.info("Refining mapping regions...")
        for chrom in mapping_regions:
            logging.debug("Refining mapping regions on chromosome {}...".format(chrom))
            mapping_regions[chrom] = refine_mapping_regions(mapping_intervals[chrom], mapping_regions[chrom],
                                                            min_transcripts=args.min_transcripts, ncpu=args.cpus)
    mapping_regions_bed = os.path.join(args.output, 'mapping_regions.bed')
    mapping_regions_to_bed(mapping_regions, mapping_regions_bed)

    # Detect chimeric genes
    logging.info("Detecting chimeric genes...")
    chimeric_genes = detect_chimeras(gff_db, mapping_regions)
    chimeras_list = os.path.join(args.output, 'chimeric_genes.list')
    print_chimeric_genes(chimeric_genes, chimeras_list)

    logging.info("Done!")


# MAIN
if __name__ == "__main__":
    main()
