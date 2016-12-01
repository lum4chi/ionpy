#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright (C) 2016 Francesco Lumachi <francesco.lumachi@gmail.com>
''' Goal of this script is to convert ThermoFisher Ion Torrent™ Variant Caller
    output (VCF file) in a TSV table, preserving much information as possible.
'''
#TODO: assuming only one sample x vcf. Hasn't tested on n_sample x vcf yet!!

import vcf as pyvcf
import pandas as pd
import numpy as np
import json, argparse, os, sys

def flattenRecord(r):
    ''' Extract from vcf interesting field: main goal is to exclude unwanted
        info but never include it explicity. This allow to NOT skip some
        unexpected info.
    '''
    # Import as possible
    variant = r.__dict__
    # Get sample info
    sample_as_val = {'GT':r.samples[0]['GT'], 'GQ':r.samples[0]['GQ']}
    # Remove unused information (TODO: there is something else to keep?)
    del(variant['FORMAT'])
    del(variant['samples'])
    del(variant['start'])
    del(variant['end'])
    del(variant['_sample_indexes'])
    del(variant['affected_end'])
    del(variant['affected_start'])
    del(variant['alleles'])
    # Pop info for flattening
    info = variant.pop('INFO', None)
    func_as_garbage = info.pop('FUNC', None)
    # Flatten FUNC
    j = json.JSONDecoder()
    func_as_json = ','.join(func_as_garbage) # pyvcf doesn't like json...
    func_as_jlist = j.decode(func_as_json.replace('\'','\"'))
    func_as_list = dict()
    for transcript in func_as_jlist:
        # Gather all the keys (need for preserving order)
        for k in transcript.keys():
            func_as_list[k] = []
    for transcript in func_as_jlist:
        # now append, if key missing on json, add ""
        for k in func_as_list.keys():
            try:
                func_as_list[k].append(transcript[k])
            except KeyError:
                func_as_list[k].append(".")
    func_as_val = {k:';'.join(v) for k,v in func_as_list.items()}
    # Flatten INFO
    info_as_list = {k:v for k,v in info.items() if type(v)==list}
    info_as_val = {k:v for k,v in info.items() if type(v)!=list}
    # Comprehension needed to stringify _Substitution type
    info_as_no_list = {k:';'.join([str(el) for el in v]) for k,v in \
                                                        info_as_list.items()}
    info_as_val.update(info_as_no_list)
    # Flatten variant info
    variant_as_list = {k:v for k,v in variant.items() if type(v)==list}
        # filter empty "FILTER"
    variant_as_list = {k:v for k,v in variant_as_list.items() if v!=[]}
    variant_as_val = {k:v for k,v in variant.items() if type(v)!=list}
    # Comprehension needed to stringify _Substitution type
    variant_as_no_list = {k:';'.join([str(el) for el in v]) for k,v in \
                                                        variant_as_list.items()}
    variant_as_val.update(variant_as_no_list)
    # Merge in a (almost) one-dim-dict
    variant_as_val.update(sample_as_val)
    variant_as_val.update(info_as_val)
    variant_as_val.update(func_as_val)
    return variant_as_val

COLUMNS_ORDER = [
'CHROM','POS','OPOS','REF','Ref','OREF','ALT','OALT','OMAPALT','TYPE','Type',
'FILTER','No Call Reason','gene','Genes','location','Location','LEN',
'Homopolymer Length','Length','HS','Info','HRUN','ID','OID','Variant ID',
'Variant Name','exon','Exon','transcript','Transcript','Strand','coding',
'Coding','codon','function','Variant Effect','Amino Acid Change','protein',
'QUAL','QD','p-value','Phred QUAL Score','DP','Coverage','FDP','AF','AO',
'% Frequency','FAO','Allele Coverage','Allele Ratio','Ref+/Ref-/Var+/Var-',
'FRO','FR','FSAF','FSAR','FSRF','FSRR','FWDB','FXX','MLLD','RBI','REFB','REVB',
'RO','SAF','SAR','SRF','SRR','SSEN','SSEP','SSSB','STB','STBP','VARB',
'grantham','Grantham','PhyloP','sift','SIFT','polyphen','PolyPhen','PFAM',
'dbSNP','DGV','MAF','EMAF','AMAF','GMAF','UCSC Common SNPs','COSMIC','OMIM',
'Gene Ontology','DrugBank','ClinVar','gt','GT','Genotype','GQ','normalizedAlt',
'normalizedPos','normalizedRef','origAlt','origPos','origRef'
]

def orderer(order, element):
    ''' Custom order for columns and push to the bottom unknown field'''
    try:
        return order.index(element)
    except:
        return np.inf

def vcf2df(a_vcf):
    ''' Take a vcf and parsing in a pandas DataFrame '''
    vcf = pyvcf.Reader(open(a_vcf))
    df = pd.DataFrame()
    # Loading
    for r in vcf:
        df = df.append(flattenRecord(r), ignore_index=True)
    # Columns ordering
    df=df[sorted(df.columns, key=lambda x:orderer(COLUMNS_ORDER, x))]
    # Type formatting
    df.POS = df.POS.astype(int)
    df.DP = df.DP.astype(int)
    df.FDP = df.FDP.astype(int)
    df.HS = df.HS.astype(bool)
    return df

# Merge information
#TODO Merge other tsv information by locus(CHROM:POS)
def mergeTable(df, table):
    tab = pd.read_table(table, comment='#')
    df['KEY'] = df.CHROM.str.cat(df.POS.map(str), sep=':')
    m = pd.merge(df, tab, how='left', left_on='KEY', right_on='Locus')
    m.drop(['KEY', 'Locus'], axis=1, inplace=True)
    return m[sorted(m.columns, key=lambda x:orderer(COLUMNS_ORDER, x))]

if __name__ == '__main__':
    # Read input
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs='+', help="List of vcf and a (optional) \
                        table to merge to it")
    parser.add_argument("-o", "--output_suffix", help="Specify a suffix to \
                        merged output files. If provided, files will be saved \
                        in same vcf input directory.")
    args = parser.parse_args()
    files = args.files
    out_suffix = args.output_suffix

    # Split files by extension and (hopefully) paired by lexicographic sorting
    vcfs = sorted([v for v in files if os.path.splitext(v)[1] == '.vcf'])
    tables = sorted([t for t in files if os.path.splitext(t)[1] != '.vcf'])

    # (soft) check on wrong pairing
    if len(tables)>0 and len(vcfs)!=len(tables):
        raise Error('Files mismatch: check correct pairing 1 vcf -> 1 table')

    # Process by pair
    for i, vcf in enumerate(vcfs):
        # Print on stdout (if multiple file, header will be repeated) or append
        # to vcf name the provided suffix and write to multiple file
        out = sys.stdout if out_suffix is None else \
                            os.path.splitext(vcf)[0] + out_suffix + '.tsv'
        print("Converting {}".format(vcf))
        vcf_df = vcf2df(vcf)                    # Convert vcf in a DataFrame
        if len(tables)>0:
            print(" + Merging with {}".format(tables[i]))
            vcf_df = mergeTable(vcf_df, tables[i])  # Enrich with paired table
        vcf_df.to_csv(out, sep='\t', index=False)
    print("Done!")
