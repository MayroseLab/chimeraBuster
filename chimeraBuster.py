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

def detect_bridges(iv_tree):
  """
  *** TODO ***
  Given an interval tree, detect
  intervals which "bridge over"
  other, non-connected intervals.
  """
  return []

def refine_mapping_regions(iv_tree, mapping_regions, max_bridges_frac=0.1, min_transcripts=2):
  """
  *** TODO ***
  Works through rough mapping regions
  and refines each by removing misleading
  intervals and 
  """
  #for mr in mapping_regions:
  #  mr_intervals = IntervalTree(iv_tree[mr.begin:mr.end])
  #  mr_bridges = detect_bridges(mr_intervals)
  refined = IntervalTree()
  for mr in mapping_regions:
    mr_transcripts = iv_tree[mr.begin:mr.end]
    if len(mr_transcripts) < min_transcripts:
      continue
    refined.add(mr)
  return refined

def intervals_frac_overlap(iv1, iv2):
  """
  Given two intervals, calculates the
  fraction of iv1 covered by iv2
  """
  return iv1.overlap_size(iv2.begin, iv2.end) / iv1.length()

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
  parser.add_argument('-b', '--max_bridges_frac', type=float, default=0.1, help='')
  parser.add_argument('-n', '--min_transcripts', type=int, default=2, help='Minimum number of transcripts in mapping region')
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
  logging.info("Analyzing transcript mapping...")
  mapping_intervals = paf_records_to_intervals(paf_hq_records)
  mapping_regions = find_mapping_regions(mapping_intervals)
  for chrom in mapping_regions:
    mapping_regions[chrom] = refine_mapping_regions(mapping_intervals[chrom], mapping_regions[chrom], min_transcripts=args.min_transcripts)

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
