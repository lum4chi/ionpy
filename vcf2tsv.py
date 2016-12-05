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

# Impose an order for known columns
KNOWN_COLUMNS = [
"Locus","CHROM","POS","OPOS","REF","Ref","OREF","ALT","OALT","OMAPALT","TYPE",
"Type","FILTER","No Call Reason","gene","Genes","location","Location","LEN",
"Homopolymer Length","Length","HS","Info","HRUN","ID","OID","Variant ID",
"Variant Name","exon","Exon","transcript","Transcript","Strand","coding",
"Coding","codon","function","Variant Effect","Amino Acid Change","protein",
"QUAL","QD","p-value","Phred QUAL Score","DP","Coverage","FDP","AF","AO",
"% Frequency","FAO","Allele Coverage","Allele Ratio","Ref+/Ref-/Var+/Var-",
"FRO","FR","FSAF","FSAR","FSRF","FSRR","FWDB","FXX","MLLD","RBI","REFB","REVB",
"RO","SAF","SAR","SRF","SRR","SSEN","SSEP","SSSB","STB","STBP","VARB",
"grantham","Grantham","PhyloP","sift","SIFT","polyphen","PolyPhen","PFAM",
"dbSNP","DGV","MAF","EMAF","AMAF","GMAF","UCSC Common SNPs","COSMIC","OMIM",
"Gene Ontology","DrugBank","ClinVar","gt","GT","Genotype","GQ","normalizedAlt",
"normalizedPos","normalizedRef","origAlt","origPos","origRef"
]

# Impose datatype when known
KNOWN_DATATYPE = {
"AF": float, "AO": int, "DP": int, "FAO": int, "FDP": int, "FR": str,
"FRO": int, "FSAF": int, "FSAR": int, "FSRF": int, "FSRR": int, "FWDB": float,
"FXX": float, "HRUN": int, "HS": bool, "LEN": int, "MLLD": float, "NS": int,
"OALT": str, "OID": str, "OMAPALT": str, "OPOS": int, "OREF": str, "PB": float,
"PBP": float, "QD": float, "RBI": float, "REFB": float, "REVB": float,
"RO": int, "SAF": int, "SAR": int, "SRF": int, "SRR": int, "SSEN": float,
"SSEP": float, "SSSB": float, "STB": float, "STBP": float, "TYPE": str,
"VARB": float, "FUNC": str, "LB": str, "exon": int, "Exon": int, "POS": int,
"Coverage": int, "GQ": int, "Grantham": int
}

def readabilify(df):
    # Support function to extract info by genotype (method change by cols)
    def chooser(position_list, genotype):
        genotype = genotype[2]
        extra_rules = {".","0"}
        if genotype in extra_rules: return position_list[0]
        return position_list[int(genotype)-1]

    def chooser3(position_list, genotype):
        genotype = genotype[2]
        extra_rules = {".","0"}
        if genotype in extra_rules: return position_list[0]
        return ', '.join([str(position_list[0]), \
                            str(position_list[int(genotype)])])
    # Start
    selection=[
    "CHROM","POS","Genotype","Ref","Type","Genes","Location","Length","Info",
    "Variant ID","Variant Name","AF","Strand","Exon","Transcript","Coding",
    "Amino Acid Change","Variant Effect","PhyloP","SIFT","Grantham","PolyPhen",
    "PFAM","dbSNP","DGV","MAF","EMAF","AMAF","GMAF","UCSC Common SNPs","COSMIC",
    "OMIM","Gene Ontology","DrugBank","ClinVar","Allele Coverage",
    "Allele Ratio","p-value","Phred QUAL Score","Coverage",
    "Ref+/Ref-/Var+/Var-","OPOS","OREF","OALT","GQ","GT"
    ]
    composite1=["OPOS","OREF","OALT","AF"]  # split by ";"
    composite3=["Allele Coverage","Allele Ratio", \
                "Ref+/Ref-/Var+/Var-"] # split by ", " and keep 1st
    # Columns subset
    df=df[selection]
    # Filtering REF&NOCALL
    df=df[(df.Type!="REF")&(df.Type!="NOCALL")]
    if len(df)==0: return df # skip if empty
    # Choise by GT
    df[composite1]=df[composite1].apply(lambda x: x.str.split(";"))
    df[composite1]=df[list(composite1)+['GT']].apply( \
        lambda x: x.iloc[:-1][~x.isnull()].apply( \
            lambda y: chooser(y,x.GT)), axis=1)
    df[composite3]=df[composite3].apply(lambda x: x.str.split(", "))
    df[composite3]=df[list(composite3)+['GT']].apply( \
        lambda x: x.iloc[:-1][~x.isnull()].apply( \
            lambda y: chooser3(y,x.GT)), axis=1)
    # Split Ref/Var
    df[composite3]=df[composite3].apply(lambda x: x.str.split(", "))
    # Separate "Ref+/Ref-/Var+/Var-"
    df[["Ref+/-","Var+/-"]]=df["Ref+/Ref-/Var+/Var-"].apply(pd.Series)
    df["Ref+/-"]=df["Ref+/-"].str.split("=").apply(lambda x: x[-1])
    df[["Ref+","Ref-"]]=df["Ref+/-"].str.split("/").apply(pd.Series)
    df["Var+/-"]=df["Var+/-"].str.split("=").apply(lambda x: x[-1])
    df[["Var+","Var-"]]=df["Var+/-"].str.split("/").apply(pd.Series)
    del df["Ref+/Ref-/Var+/Var-"]
    del df["Ref+/-"]
    del df["Var+/-"]
    # Separate Allele info
    df[["Allele Coverage Ref","Allele Coverage Alt"]]= \
        df["Allele Coverage"].apply(pd.Series)
    df[["Allele Ratio Ref","Allele Ratio Alt"]]= \
        df["Allele Ratio"].apply(pd.Series)
    df["Allele Coverage Ref"]= \
        df["Allele Coverage Ref"].str.split("=").apply(lambda x: x[-1])
    df["Allele Coverage Alt"]= \
        df["Allele Coverage Alt"].str.split("=").apply(lambda x: x[-1])
    df["Allele Ratio Ref"]= \
        df["Allele Ratio Ref"].str.split("=").apply(lambda x: x[-1])
    df["Allele Ratio Alt"]= \
        df["Allele Ratio Alt"].str.split("=").apply(lambda x: x[-1])
    del df["Allele Coverage"]
    del df["Allele Ratio"]
    # AF -> % Frequency
    df["% Frequency"]=df.AF
    del df["AF"]
    # Heterozigosity flagging
    df["Genotype Flag"]=df.GT.str.split("/|\|").apply( \
        lambda x: x[0]==x[1]).apply(lambda x: "Hom" if x else "Het")
    # End
    return df

def flattenRecord(r):
    ''' Extract from vcf interesting field: main goal is to exclude unwanted
        info but never include it explicity. This allow to NOT skip some
        unexpected info.
        Assuming only one sample, can be divide single array as a value and
        composite value as list of values separated by ";" (to increase
        readability and compatibility to work with other apps, i.e. R)
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

def orderer(order, element):
    ''' Custom order for columns, pushing to the bottom unknown field'''
    try:
        return order.index(element)
    except:
        return np.inf

# Try to impose a type where possible
def formattify(df):
    for col in df:
        try:
            df[col] = df[col].astype(KNOWN_DATATYPE[col])
        except: continue
    return df

def vcf2df(a_vcf):
    ''' Take a vcf and parsing in a pandas DataFrame '''
    vcf = pyvcf.Reader(open(a_vcf))
    df = pd.DataFrame()
    # Loading
    for r in vcf:
        df = df.append(flattenRecord(r), ignore_index=True)
    # Columns ordering
    df=df[sorted(df.columns, key=lambda x:orderer(KNOWN_COLUMNS, x))]
    return df

# Merge information
def mergeTable(df, table):
    t = pd.read_table(table, comment='#')
    df['KEY'] = df.CHROM.str.cat(df.POS.map(str), sep=':')
    m = pd.DataFrame(df.merge(t, how='left', left_on='KEY', right_on='Locus'))
    m.drop(['KEY'], axis=1, inplace=True)
    return m[sorted(m.columns, key=lambda x:orderer(KNOWN_COLUMNS, x))]

if __name__ == '__main__':
    # Read input
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs='+', help="List of vcf and a (optional) \
                        table to merge to it")
    parser.add_argument("-o", "--output_suffix", help="Specify a suffix to \
                        merged output files. If provided, files will be saved \
                        in same vcf input directory.")
    parser.add_argument("-f", "--filter", dest="toReadabilify",
                        action="store_true", help="Subset/Filter/Genotyping \
                        according to 'readabilify' function (see code)")
    parser.set_defaults(toReadabilify=False)
    args = parser.parse_args()
    files = args.files
    out_suffix = args.output_suffix
    toReadabilify = args.toReadabilify

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
        vcf_df = formattify(vcf_df)
        if len(tables)>0:                       # Enrich with paired table
            print(" + Merging with {}".format(tables[i]))
            vcf_df = mergeTable(vcf_df, tables[i])
        if toReadabilify:                       # "Clean" version
            print(" + Apply cleaning...")
            vcf_df = readabilify(vcf_df)
        vcf_df = formattify(vcf_df)
        vcf_df.to_csv(out, sep='\t', index=False)
    print("Done!")
