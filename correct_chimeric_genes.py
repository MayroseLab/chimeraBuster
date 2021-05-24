#!/usr/bin/env python

"""
chimeraBuster is a tool for the detection and correction
of chimeric gene annotations.
This script is the second step of the workflow, aimed
at breaking and correcting chimeric genes by re-annotation
of relevant regions.
For further information please visit the tool Github:
https://github.com/MayroseLab/chimeraBuster
"""

# Built-in/Generic Imports
from argparse import ArgumentParser
import os
import sys
import logging
import re
from collections import OrderedDict

# Libs
import gffutils

from util import *

__author__ = "Lior Glick"
__license__ = "MIT"
__version__ = "0.2.0"
__maintainer__ = "Lior Glick"
__status__ = "Dev"


# FUNCTIONS

def create_mapping_regions_bed(chimeric_list, bed, flank=50):
    """
    Convert chimeric genes list to bed format
    """
    with open(chimeric_list) as f, open(bed,'w') as fo:
        for line in f:
            gene_id, chrom, regions = line.strip().split('\t')
            i = 1
            for reg in regions.split(','):
                start, end = [int(x) for x in reg.split('-')]
                start = str(max(0, start-flank))
                end = str(end + flank)
                new_gene_id = '%s_breakChimera%s' %(gene_id, i)
                print('\t'.join([chrom, start, end, new_gene_id]), file=fo)
                i += 1

def create_PASA_conf(template, out_conf, replace_dict):
    """
    Replace keys with values of replace_dict
    in template, and write to out_conf
    """
    with open(template) as f, open(out_conf, 'w') as fo:
        text = f.read()
        for k,v in replace_dict.items():
            text = text.replace(k, v)
        print(text, file=fo)

def parse_seqid(feature):
    """
    Given a gff feature object with
    the expected seqid format: 
    <orig gene name>_breakChimera<mapping region index>::<chromosome>:<start>-<end>
    return gene_id, chrom, mapping_region_start
    """
    seqid_regex = re.compile(r'(.+_breakChimera\d+)::(\w+):(\d+)-(\d+)')
    res = seqid_regex.search(feature.seqid)
    if not res:
        logging.error("Can't convert GFF feature with seqid %s. Unexpected seqid format" % feature.seqid)
        sys.exit()
    gene_id, chrom, mr_start, mr_end = res.groups()
    mr_start = int(mr_start)
    return gene_id, chrom, mr_start

def convert_to_genome_coords(gff_db):
    """
    Convert from mapping region coordinates
    and PASA gene IDs to genome coordinates
    and chimeraBuster IDs.
    Returns a list of gff lines.
    """
    gff_lines = []
    gi = 1
    for gene in gff_db.features_of_type('gene', order_by='start'):
        # convert gene feature
        gene_id, chrom, mr_start = parse_seqid(gene)
        gene_id += '_%s' % gi
        gi += 1
        gene.seqid = chrom
        gene.start += mr_start
        gene.end += mr_start
        gene.attributes = {'ID': [gene_id]}
        gene.source = 'chimeraBuster'
        gff_lines.append(str(gene))
        # convert mRNA - only use longest mRNA
        mrna_features = list(gff_db.children(gene, featuretype='mRNA', order_by='start'))
        mrna = sorted(mrna_features, key=lambda feat: feat.end - feat.start)[-1]
        mrna.seqid = chrom
        mrna.source = 'chimeraBuster'
        mrna.start += mr_start
        mrna.end += mr_start
        mRNA_id = gene_id + '_mRNA'
        mrna.attributes = {'ID': [mRNA_id], 'Parent': [gene_id]}
        gff_lines.append(str(mrna))
        # convert other features
        feture_index = {'exon': 1, 'CDS': 1}
        for feature in gff_db.children(mrna, featuretype=['exon','CDS','five_prime_UTR','three_prime_UTR'], order_by='start'):
            feature.seqid = chrom
            feature.source = 'chimeraBuster'
            feature.start += mr_start
            feature.end += mr_start
            if feature.featuretype in feture_index:
                feature_id = gene_id + '_%s_%s' %(feature.featuretype, feture_index[feature.featuretype])
                feture_index[feature.featuretype] += 1
            else:
                feature_id = gene_id + '_%s' % feature.featuretype
            feature.attributes = {'ID': [feature_id], 'Parent': [mRNA_id]}
            gff_lines.append(str(feature))
    return gff_lines

def write_final_gff(orig_gff_db, chimeric_genes_list, new_lines, final_gff):
    """
    Remove chimeric genes and add
    new lines to create final gff
    """
    chimeric_gene_ids = set([line.strip().split('\t')[0] for line in open(chimeric_genes_list)])
    with open(final_gff, 'w') as fo:
        print('##gff-version 3', file=fo)
        for gene in orig_gff_db.features_of_type('gene'):
            if gene['ID'][0] not in chimeric_gene_ids:
                print(str(gene), file=fo)
                for feat in orig_gff_db.children(gene):
                    print(str(feat), file=fo)
        print('\n'.join(new_lines), file=fo)
    

def main():
    # command line arguments
    DESC = "correct_chimeric_genes.py: break chimeric gene models by re-annotation."
    USAGE = "python correct_chimeric_genes.py -h (display full usage doc)"
    parser = ArgumentParser(description=DESC, usage=USAGE)
    parser.add_argument('-g', '--gff', help='Input GFF', type=lambda x: is_valid_file(parser, x), default=None)
    parser.add_argument('-f', '--genome_fasta', help='Input genome in FASTA format', required=True,
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument('-t', '--transcripts', required=True,
                        help='Transcripts to be mapped to the genome, in FASTA format',
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument('-l', '--chimeric_list', required=True,
                        help='List of chimeric genes, created by detect_chimeric_genes.py',
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument('-b', '--flank_size', help='Flank size (bp) around mapping regions', type=int, default=50)
    parser.add_argument('-o', '--output', help='Output path', required=True)
    parser.add_argument('-c', '--cpus', type=int, default=1, help='Number of CPUs to use')
    parser.add_argument('-v', '--verbose', action="store_true", default=False, help='Increase verbosity')
    args = parser.parse_args()

    # logging
    log_format = '(%(processName)-10s) %(asctime)s %(message)s'
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG, format=log_format)
    else:
        logging.basicConfig(level=logging.INFO, format=log_format)

    if not os.path.exists(args.output):
        os.mkdir(args.output)
    
    # Create mapping regions bed
    logging.info("Creating mapping regions bed...")
    mapping_regions_bed = os.path.join(args.output, 'chimeric_mapping_regions.bed')
    create_mapping_regions_bed(args.chimeric_list, mapping_regions_bed, args.flank_size)

    # Extract mapping regions sequences
    logging.info("Extracting mapping regions sequences...")
    mapping_regions_fasta = os.path.join(args.output, 'chimeric_mapping_regions.fasta')
    run_external_cmd("bedtools getfasta -fi {} -bed {} -name".format(args.genome_fasta, mapping_regions_bed), mapping_regions_fasta)

    # Predict genes in mapping regions (using PASA)
    logging.info("Predicting gene structures...")
    # file paths required by PASA
    pasa_db = os.path.realpath(os.path.join(args.output, 'pasa_db.sqlite'))
    pasa_align_conf = os.path.join(args.output, 'pasa_align.conf')
    pasa_annotate_conf = os.path.join(args.output, 'pasa_annotate.conf')
    pasa_fake_gff = os.path.join(args.output, 'fake.gff3')
    # PASA scripts and configs
    conda_env_path = os.environ['CONDA_PREFIX']
    pasa_script = os.path.join(conda_env_path, 'opt/pasa-2.4.1/Launch_PASA_pipeline.pl')
    pasa_conf_templates = os.path.join(conda_env_path, 'opt/pasa-2.4.1/pasa_conf')
    pasa_load_script = os.path.join(conda_env_path, 'opt/pasa-2.4.1/scripts/Load_Current_Gene_Annotations.dbi')
    # create pasa alignment conf
    pasa_align_conf_template = os.path.join(pasa_conf_templates, 'pasa.alignAssembly.Template.txt')
    pasa_conf_replace = {'<__DATABASE__>':pasa_db, '<__MIN_PERCENT_ALIGNED__>':'75', '<__MIN_AVG_PER_ID__>': '90'}
    create_PASA_conf(pasa_align_conf_template, pasa_align_conf, pasa_conf_replace)
    # run PASA transcript align and assembly
    pasa_align_out = os.path.join(args.output, 'pasa_align.out')
    pasa_align_cmd = "{} -c {} -C -R -g {} -t {} --ALIGNERS blat --TRANSDECODER --CPU {}".format(pasa_script, os.path.realpath(pasa_align_conf), os.path.realpath(mapping_regions_fasta), os.path.realpath(args.transcripts), args.cpus)
    err = run_external_cmd(pasa_align_cmd, pasa_align_out, cwd=args.output)
    if err:
        logging.error("An error occurred while running PASA. Error message: {}".format(err))
        sys.exit()
    # create fake gff towards PASA annotation
    with open(pasa_fake_gff, 'w') as fo:
        print('# fake annotation', file=fo)
    # load fake gff
    load_cmd = "{} -c {} -g {} -P {}".format(pasa_load_script, pasa_align_conf, mapping_regions_fasta, pasa_fake_gff)
    load_out = os.path.join(args.output, 'pasa_load.out')
    err = run_external_cmd(load_cmd, load_out)
    if err:
        logging.error("An error occurred while running PASA. Error message: {}".format(err))
        sys.exit()
    # create PASA annotation conf
    pasa_annotate_conf_template = os.path.join(pasa_conf_templates, 'pasa.annotationCompare.Template.txt')
    pasa_conf_replace = {'<__DATABASE__>':pasa_db}
    create_PASA_conf(pasa_annotate_conf_template, pasa_annotate_conf, pasa_conf_replace)
    # run PASA annotation
    pasa_annotate_out = os.path.join(args.output, 'pasa_annotate.out')
    pasa_annotate_cmd = "{} -c {} -A -g {} -t {} --CPU {}".format(pasa_script, os.path.realpath(pasa_annotate_conf), os.path.realpath(mapping_regions_fasta), os.path.realpath(args.transcripts), args.cpus)
    err = run_external_cmd(pasa_annotate_cmd, pasa_annotate_out, cwd=args.output)
    if err:
        logging.error("An error occurred while running PASA. Error message: {}".format(err))
        sys.exit()
    
    # Create GFF of new gene models
    logging.info("Writing corrected gene models to GFF...")
    # locate and load PASA output gff
    pasa_mr_gff = [f for f in os.listdir(args.output) if f.startswith('pasa_db.sqlite.gene_structures_post_PASA_updates') and f.endswith('.gff3')][0]
    pasa_mr_gff = os.path.join(args.output, pasa_mr_gff)
    # PASA sometimes includes duplicate records, to avoid this, pass the gff lines through an OrderedDict
    pasa_mr_gff_lines = OrderedDict.fromkeys([line.strip() for line in open(pasa_mr_gff)])
    pasa_mr_gff_text = '\n'.join(pasa_mr_gff_lines.keys())
    db_path = pasa_mr_gff + '.db'
    gff_db = gffutils.create_db(pasa_mr_gff_text, dbfn=db_path, keep_order=True, merge_strategy='create_unique', from_string=True)
    gff_db = gffutils.FeatureDB(db_path, keep_order=True)
    # convert coordinates and IDs to create final GFF
    chimeric_corrected_gff = os.path.join(args.output, 'chimeric_genes_corrected.gff')
    chimeric_corrected_gff_lines = convert_to_genome_coords(gff_db)
    with open(chimeric_corrected_gff, 'w') as fo:
        print('##gff-version 3', file=fo)
        print('\n'.join(chimeric_corrected_gff_lines), file=fo)
    
    # Create final GFF (remove chimeric genes and insert the new ones)
    logging.info("Writing final annotation GFF...")
    # connect to original gff DB
    orig_gff_db_path = os.path.join(args.output, os.path.basename(args.gff) + '.db')
    orig_gff_db = gffutils.FeatureDB(orig_gff_db_path, keep_order=True)
    # write GFF
    out_gff = os.path.join(args.output, os.path.basename(args.gff) + '.chimeraBuster.corrected.gff')
    write_final_gff(orig_gff_db, args.chimeric_list, chimeric_corrected_gff_lines, out_gff)
    
    logging.info("Done!")

# MAIN
if __name__ == "__main__":
    main()
