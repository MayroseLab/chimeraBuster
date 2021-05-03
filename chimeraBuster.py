"""

"""

from argparse import ArgumentParser
import os
import sys
import logging
import gffutils
from intervaltree import IntervalTree, Interval
from pafpy import PafFile
#from Bio import SeqIO
from copy import deepcopy
from itertools import combinations, chain
from sklearn.cluster import DBSCAN
import pandas as pd
from multiprocessing import Pool

def run_minimap(genome, transcripts):
  return

def filter_paf(paf_path):
  """
  *** TODO ***
  Filter a PAF file according to
  specified criteria. Return a list
  of PAF records
  """
  records = []
  with PafFile(paf_path) as paf:
    for record in paf:
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
  chimeras = []
  for gene in gff_db.features_of_type('gene'):
    gene_id = gene['ID'][0]
    chrom = gene.seqid
    intersect_mappings = mapping_regions[chrom][gene.start:gene.end]
    gene_iv = Interval(gene.start, gene.end, gene_id)
    intersect_mappings = [im for im in intersect_mappings if intervals_frac_overlap(im,gene_iv) >= min_intersect_frac]
    n_im = len(intersect_mappings)
    logging.debug("Gene {} : {} intersecting mapping regions".format(gene_id, n_im))
    if n_im > 1:
      chimeras.append(gene_id)
  return chimeras
    
def is_valid_file(parser, arg):
  if not os.path.isfile(arg):
    parser.error("The file %s does not exist!" % arg)
  else:
    return arg

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
  #chrom_seq = SeqIO.to_dict(SeqIO.parse(args.genome_fasta, "fasta"))

  # Detect chimeric genes
  logging.info("Detecting chimeric genes...")
  chimeric_genes = detect_chimeras(gff_db, mapping_regions)
  chimeras_list = os.path.join(args.output, 'chimeric_genes.list')
  with open(chimeras_list, 'w') as cfo:
    print('\n'.join(chimeric_genes), file=cfo)

  # Correct chimeras
  #logging.info("Correcting chimeric genes...")
  
  # Print output

if __name__ == "__main__":
  main()
