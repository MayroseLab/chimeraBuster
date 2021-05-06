"""

"""

from argparse import ArgumentParser
import os
import sys
import logging
import gffutils
from intervaltree import IntervalTree, Interval
from pafpy import PafFile
from Bio import SeqIO
from copy import deepcopy
from itertools import combinations, chain
from sklearn.cluster import DBSCAN
import pandas as pd
from multiprocessing import Pool

def run_minimap(genome, transcripts):
  return

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
    iv_dict[rec.tname].add(Interval(rec.tstart,rec.tend,[rec.qname]))
  return iv_dict

def find_mapping_regions(iv_dict):
  """
  Detect rough regions in which
  transcripts were mapped. These
  regions will be refined later.
  """
  res = deepcopy(iv_dict)
  for chrom in res:
    res[chrom].merge_overlaps(data_reducer=lambda iv1,iv2: iv1+iv2)
  n_mr = sum([len(res[chrom]) for chrom in res])
  logging.info("{} mapping regions detected.".format(n_mr))
  return res

def intervals_frac_overlap(iv1, iv2):
  """
  Given two intervals, calculates the
  fraction of iv1 covered by iv2
  """
  return iv1.overlap_size(iv2.begin, iv2.end) / iv1.length()

def is_bridge(iv1, iv2, iv3):
  """
  Given 3 intervals, determines if
  iv1 is a bridge connecting iv2 and iv3.
  To be considered a bridge, iv2 and iv3
  must *not* overlap, and iv1 must
  connect them
  """
  # check if iv2 and iv3 overlap
  if intervals_frac_overlap(iv2, iv3) > 0:
    return False
  # check if iv1 is bridge
  if intervals_frac_overlap(iv1, iv2) > 0 and intervals_frac_overlap(iv1, iv3) > 0:
    return True
  return False

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
  clustering = DBSCAN(eps=dbscan_eps, min_samples=dbscan_minPts).fit(df[['start','end']])
  df['cluster'] = pd.Series(clustering.labels_).astype(str)
  # intervals with cluster '-1' are outliers - suspected bridges or other
  n_outliers = df.query('cluster == "-1"').shape[0]
  n_iv = len(intervals)
  outlier_frac = n_outliers / n_iv
  logging.debug("{} out of {} intervals (fraction = {}) detected as outliers".format(n_outliers, n_iv, outlier_frac))
  if outlier_frac > max_outlier_frac:
    return intervals
  return set(df.query('cluster != "-1"')['interval'])

#def detect_bridges(intervals, dbscan_eps=300, dbscan_minPts=4):
#  """
#  Given a set of intervals, detect
#  intervals which "bridge over"
#  other, non-connected intervals.
#  """
#  transcript_pairs = list(combinations(intervals, 2))
#  for pair in transcript_pairs:
#    if intervals_frac_overlap(pair[0], pair[1]):
#      continue
#    for interval in intervals:
#      if is_bridge(interval, pair[0], pair[1]):
#        pair_t = (pair[0], pair[1])
#        if pair_t not in bridges:
#          bridges[pair_t] = []
#        bridges[pair_t].append(interval)
#  return bridges

def refine_mapping_regions(iv_tree, mapping_regions, min_transcripts=2, ncpu=1):
  """
  Works through rough mapping regions
  and refines each by removing misleading
  intervals and 
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

def detect_chimeras(gff_db, mapping_regions, min_intersect_frac=0.5):
  """
  Go through all genes and detect
  putative chimeras.
  min_intersect_frac is the minimal
  fraction of the transcript length
  that overlaps with the gene to be
  considered.
  """
  chimeras = {}
  for gene in gff_db.features_of_type('gene'):
    gene_id = gene['ID'][0]
    chrom = gene.seqid
    if chrom not in mapping_regions:
      continue
    intersect_mappings = mapping_regions[chrom][gene.start:gene.end]
    gene_iv = Interval(gene.start, gene.end, gene_id)
    intersect_mappings = [im for im in intersect_mappings if intervals_frac_overlap(im,gene_iv) >= min_intersect_frac]
    n_im = len(intersect_mappings)
    logging.debug("Gene {} : {} intersecting mapping regions".format(gene_id, n_im))
    if n_im > 1:
      chimeras[gene_id] = intersect_mappings
  return chimeras
    
def is_valid_file(parser, arg):
  if not os.path.isfile(arg):
    parser.error("The file %s does not exist!" % arg)
  else:
    return arg

def closest(lst, K):
  """Find closest number to K in list"""
  return lst[min(range(len(lst)), key = lambda i: abs(lst[i]-K))]

def find_nearest_exon(gff_db, gene_feature, coordinate, search_direction=None):
  """
  Given a gene feature and a genomic
  coordinate, find the nearest exon
  of that gene. If search_direction
  is None, look in both directions.
  If 'upstrea' or 'downstream' - only
  search in this direction.
  """
  assert search_direction in {None, 'upstream', 'downstream'} , "search_direction must be None or 'upstream' or'downstream'"
  # find all gene exons - save as list of intervals
  gene_exons = [Interval(exon.start, exon.end, exon) for exon in gff_db.children(gene_feature, featuretype='exon', order_by='start')]
  # if search direction turned on
  if search_direction == 'downstream':
    gene_exons = [exon for exon in gene_exons if exon.end >= coordinate]
  elif search_direction == 'upstream':
    gene_exons = [exon for exon in gene_exons if exon.begin <= coordinate]
  # find nearest
  if (search_direction == 'downstream' and gene_exons[-1].end < coordinate) or (search_direction == 'upstream' and gene_exons[0].end > coordinate):
    return None
  exons_start_end = [[exon.begin, exon.end] for exon in gene_exons]
  exons_start_end = list(chain(*exons_start_end))	# list of all start and end positions
  nearest_point = closest(exons_start_end, coordinate)	# find closest start/end
  nearest_exon = IntervalTree(gene_exons)[nearest_point-1:nearest_point+1].pop().data	# +/- 1 to avoid issues with including interval end
  return nearest_exon

def correct_exon(exon_feature, chrom_seq, max_coord=None, min_coord=None):
  """
  Correct an exon by finding the earliest
  stop codon between the exon start and the
  gap end.
  exon_feature : gffutils exon
  max_coord : highest genomic position to extend
    the exon to when searching for stop codons
  min_coord : lowest genomic position to extend
    the exon to when searching for start codons
  chrom_seq : SeqRecord object of the chromosome
  Returns a modified exon feature or, if
  no stop codon found - the original exon.
  """
  assert (max_coord or min_coord) and not (max_coord and min_coord), "Must provide max_coord or min_coord, but not both"
  stop_codons = {'TAA','TAG','TGA'}
  start_codon = {'ATG'}
  if exon_feature.strand == '+':
    if max_coord:	# exon before gap
      search_region = (exon_feature.start, max_coord)
      search_for = stop_codons
    elif min_coord:	# exon after gap
      search_region = (min_coord, exon_feature.end)
      search_for = start_codon
  elif exon_feature.strand == '-':
    if min_coord:	# exon before gap
      search_region = (min_coord, exon_feature.end)
      search_for = stop_codons
    elif max_coord:	# exon after gap
      search_region = (exon_feature.start, max_coord)
      search_for = start_codon

  search_region_seq = chrom_seq[search_region[0]:search_region[1]+1]
  if exon_feature.strand == '-':
    search_region_seq = search_region_seq.reverse_complement()
  codons = [search_region_seq[i:i+3] for i in range(0,len(search_region_seq),3)]
  if min_coord:
    codons.reverse()

  extend_to_codon = None	# index of codon where the corrected exon should end
  i = 0
  for codon in codons:
    if str(codon) in search_for:
      extend_to_codon = i
      break
    i += 1
  if not extend_to_codon:
    return exon_feature
  exon_feature_corrected = deepcopy(exon_feature)
  if max_coord:
    exon_feature_corrected.end = exon_feature_corrected.start + i*3
  elif min_coord:
    exon_feature_corrected.start = exon_feature_corrected.end - 3*i
  return exon_feature_corrected

def split_gene(gff_db, gene_feature, split_coords=[]):
  """
  Splits a gene feature into multiple
  gene features by cutting at given
  genomic coordinates. split_coords
  is a list of (start,end) pairs.
  """
  new_features = []
  i = 1
  for pair in split_coords:
    start, end = pair
    # create the gene feature
    new_gene = deepcopy(gene_feature)
    new_gene_id = new_gene['ID'][0] + '_%s' % i
    new_gene['ID'][0] = new_gene_id
    new_gene.start = start
    new_gene.end = end
    new_gene.source = 'chimeraBuster'
    # create mRNA feature
    new_mrna_id = new_gene_id + '_mRNA_1'
    new_mrna = deepcopy(new_gene)
    new_mrna.featuretype = 'mRNA'
    new_mrna['ID'][0] = new_mrna_id
    new_mrna['Parent'] = [new_gene_id]
    # get exon and CDS features within the start-end coordinates
    all_features = [deepcopy(feat) for feat in gff_db.region(region=(gene_feature.seqid, start, end))]
    exons = list(filter(lambda feat: feat.featuretype == "exon", all_features))
    cds = list(filter(lambda feat: feat.featuretype == "CDS", all_features))
    exons.sort(key=lambda feat: feat.start)
    cds.sort(key=lambda feat: feat.start)
    # set start and end coordinates for first and last exon/CDS
    exons[0].start = start
    exons[-1].end = end
    cds[0].start = start
    cds[-1].end = end
    # modify exon and CDS features
    for j in range(len(exons)):
      exons[j]['ID'][0] = '%s_exon_%s' %(new_mrna_id, j+1)
      exons[j]['Parent'][0] = str(new_mrna_id)
    for j in range(len(cds)):
      cds[j]['ID'][0] = '%s_CDS_%s' %(new_mrna_id, j+1)
      cds[j]['Parent'][0] = str(new_mrna_id)

    # add new features
    new_features.extend([new_gene, new_mrna] + exons + cds)
    i += 1

  return new_features

def correct_chimeric_gene(gene_id, chimeric_genes, gff_db, mapping_regions, genome_seq_dict):
  """
  Break a chimeric gene according to
  mapping regions and fix the ends
  of resulting genes
  """
  mapping_regions = IntervalTree(chimeric_genes[gene_id])

  # find the gap region/s between MRs
  gene_feature = gff_db[gene_id]
  gene_iv = Interval(gene_feature.start, gene_feature.end, gene_id)
  gap_regions = IntervalTree([gene_iv])
  for mr in mapping_regions:
    gap_regions.chop(mr.begin, mr.end)
  # remove terminal gaps
  term = []
  for gap in gap_regions:
    if gap.begin == gene_feature.start or gap.end == gene_feature.end:
      term.append(gap)
  for term_gap in term:
    gap_regions.remove(term_gap)
  # sort gap regions by position (5' to 3')
  rev = (gene_feature.strand == '-')	# True or False
  gap_regions = sorted(list(gap_regions), key=lambda iv: iv.begin, reverse=rev)

  # split gene
  chrom_seq = genome_seq_dict[gene_feature.seqid]
  
  exons_before_gaps = []
  exons_after_gaps = []
  for gap in gap_regions:
    # find last exon before the gap and first exon after the gap
    # (these exons will be corrected)
    if gene_feature.strand == '+':
      exon_before_gap = find_nearest_exon(gff_db, gene_feature, gap.begin, search_direction='upstream')
      exon_after_gap = find_nearest_exon(gff_db, gene_feature, gap.end, search_direction='downstream')
    elif gene_feature.strand == '-':
      exon_before_gap = find_nearest_exon(gff_db, gene_feature, gap.end, search_direction='downstream')
      exon_after_gap = find_nearest_exon(gff_db, gene_feature, gap.begin, search_direction='upstream')     
    if not (exon_before_gap and exon_after_gap):
      # if for some reason no exon was found - skip and do not split on this gap
      continue
    # correct exons
    if gene_feature.strand == '+':
      exon_before_gap_corrected = correct_exon(exon_before_gap, chrom_seq, max_coord=gap.end)
      exon_after_gap_corrected = correct_exon(exon_after_gap, chrom_seq, min_coord=exon_before_gap_corrected.end)
    elif gene_feature.strand == '-':
      exon_before_gap_corrected = correct_exon(exon_before_gap, chrom_seq, min_coord=gap.begin)
      exon_after_gap_corrected = correct_exon(exon_after_gap, chrom_seq, max_coord=exon_before_gap_corrected.start)
    exons_before_gaps.append(exon_before_gap_corrected)
    exons_after_gaps.append(exon_after_gap_corrected) 
  # correct gene according to corrected exons
  # create start/end coordinates for gene splitting
  if gene_feature.strand == '+':
    splits_starts = [gene_feature.start] + [eag.start for eag in exons_after_gaps]
    splits_ends = [ebg.end for ebg in exons_before_gaps] + [gene_feature.end]
  elif gene_feature.strand == '-':
    splits_starts = [gene_feature.start] + [ebg.start for ebg in reversed(exons_before_gaps)]
    splits_ends = [eag.end for eag in reversed(exons_after_gaps)] + [gene_feature.end]
  split_coords = list(zip(splits_starts, splits_ends))
  # split gene (and adjust children features)
  new_features = split_gene(gff_db, gene_feature, split_coords=split_coords)
  # create list of features to remove (original gene + children features)
  remove_features = [gene_feature] + list(gff_db.children(gene_feature))

  return new_features, remove_features

def update_gff_db(gff_db, new_features, features_to_remove, new_db_path):
  """
  Create a new GFF DB. Add all features
  from gff_db, except those listed in
  features_to_remove, and add new features
  from new_features.
  This function is only needed because the
  gffutils.FeatureDB.update() method doesn't
  work (sqlite3 DB locked error)
  """
  features_from_orig = list(filter(lambda feat: feat not in features_to_remove, gff_db.all_features()))
  features_to_add = features_from_orig + new_features
  new_gff_db = gffutils.create_db(data=features_to_add, dbfn=new_db_path, merge_strategy='create_unique')
  new_gff_db = gffutils.FeatureDB(new_db_path)
  return new_gff_db


def main():

  # command line arguments
  DESC = "ChimeraBuster: break chimeric gene models based on transcript mapping"
  USAGE = "python chimeraBuster.py -h (display full usage doc)"
  parser = ArgumentParser(description=DESC, usage=USAGE)
  parser.add_argument('-g', '--gff', help='Input GFF', type=lambda x: is_valid_file(parser, x), default=None)
  parser.add_argument('-f', '--genome_fasta', help='Input genome in FASTA format', required=True, type=lambda x: is_valid_file(parser, x))
  parser.add_argument('-t', '--transcripts', default=None, help='Transcripts to be mapped to the genome, in FASTA format', type=lambda x: is_valid_file(parser, x))
  parser.add_argument('-m', '--mapping', default=None, help='Transcripts mapping to the genome, in minimap2 PAF format', type=lambda x: is_valid_file(parser, x))
  parser.add_argument('-d', '--gff_db', default=None, help='GFF DB from a previous run', type=lambda x: is_valid_file(parser, x))
  parser.add_argument('-o', '--output', help='Output path', required=True)
  parser.add_argument('-r', '--do_not_refine', help='Skip mapping region refining step for a quick-and-dirty analysis', action='store_true', default=False)
  parser.add_argument('-n', '--min_transcripts', type=int, default=2, help='Minimum number of transcripts in mapping region')
  parser.add_argument('-q', '--min_query_cov', type=float, default=0.95, help='Minimum transcript query coverage (0-1)')
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
  
  # Run minimap, if needed
  if args.transcripts:
    logging.info("Running Minimap2...")
    try:
      args.mapping = run_minimap(args.genome_fasta, args.transcripts)
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
      mapping_regions[chrom] = refine_mapping_regions(mapping_intervals[chrom], mapping_regions[chrom], min_transcripts=args.min_transcripts, ncpu=args.cpus)

  # Read GFF
  if args.gff:
    logging.info("Loading GFF from file...")
    db_path = os.path.join(args.output, os.path.basename(args.gff) + '.db')
    gff_db = gffutils.create_db(args.gff, dbfn=db_path, keep_order=True, merge_strategy='create_unique', sort_attribute_values=True, verbose=args.verbose)
  else:
    logging.info("Loading GFF from DB...")
    db_path = args.gff_db
  gff_db = gffutils.FeatureDB(db_path, keep_order=True)

  # Read genome fasta
  logging.info("Loading genome FASTA...")
  chrom_seq = SeqIO.to_dict(SeqIO.parse(args.genome_fasta, "fasta"))

  # Detect chimeric genes
  logging.info("Detecting chimeric genes...")
  chimeric_genes = detect_chimeras(gff_db, mapping_regions)
  chimeras_list = os.path.join(args.output, 'chimeric_genes.list')
  with open(chimeras_list, 'w') as cfo:
    print('\n'.join(chimeric_genes.keys()), file=cfo)

  # Correct chimeras
  logging.info("Correcting chimeric genes...")
  new_features = []
  remove_features = []
  for chimera in chimeric_genes.keys():
    chimera_split_new_features, chimera_split_remove_features = correct_chimeric_gene(chimera, chimeric_genes, gff_db, mapping_regions, chrom_seq)
    new_features.extend(chimera_split_new_features)
    remove_features.extend(chimera_split_remove_features)

  # Create GFF DB with corrected chimeras
  logging.info("Updating GFF DB...")
  corrected_db_path = os.path.join(args.output, os.path.basename(args.gff) + '.corrected.db')
  corrected_gff_db = update_gff_db(gff_db, new_features, set(remove_features), corrected_db_path)
  
  # Print output
  logging.info("Writing corrected GFF...")
  out_gff_path = os.path.join(args.output, os.path.basename(args.gff) + '.corrected.gff')
  with open(out_gff_path, 'w') as fo:
    for feat in corrected_gff_db.all_features():
      print(str(feat), file=fo)

  logging.info("chimeraBuster is done!")

if __name__ == "__main__":
  main()
