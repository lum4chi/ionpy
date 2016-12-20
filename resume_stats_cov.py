#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright (C) 2016 Francesco Lumachi <francesco.lumachi@gmail.com>
''' Goal of this script is to resume ThermoFisher Ion Torrent™ Coverage Analysis
    output (*.stats.cov.txt files) in a CSV table. Those files are usually
    stored in:
        /results/analysis/output/Home/<run>/plugin_out/<coverageAnalysis_out>/*
'''

import pandas as pd
import argparse, sys

COL2FLOAT = [
    "Amplicons reading end-to-end",
    "Amplicons with at least 1 read",
    "Amplicons with at least 100 reads",
    "Amplicons with at least 20 reads",
    "Amplicons with at least 500 reads",
    "Amplicons with no strand bias",
    "Average base coverage depth",
    "Average reads per amplicon",
    "Percent assigned amplicon reads",
    "Percent base reads on target",
    "Percent end-to-end reads",
    "Percent reads on target",
    "Target base coverage at 100x",
    "Target base coverage at 1x",
    "Target base coverage at 20x",
    "Target base coverage at 500x",
    "Target bases with no strand bias",
    "Uniformity of amplicon coverage",
    "Uniformity of base coverage"
    ]

COL2INT = [
    "Bases in target regions",
    "Number of amplicons",
    "Number of mapped reads",
    "Total aligned base reads",
    "Total assigned amplicon reads",
    "Total base reads on target"
    ]

def txt2series(f):
    f = open(f)
    rows = f.readlines()[2:]
    rows_as_kv = [r.strip().split(":") for r in rows]
    d = {k_v[0].strip():k_v[1].strip().strip("%") for k_v in rows_as_kv if len(k_v)>1}
    # Separate barcode_run
    d['Barcode'] = d['Alignments'][:13]
    d['Run'] = d['Alignments'][14:]
    del d['Alignments']
    s = pd.Series(d)
    return s

if __name__ == '__main__':
    # Read input
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs='+', help="List of stats.cov.txt")
    args = parser.parse_args()
    files = args.files
    out = sys.stdout

    # Assert input only files ending with ".stats.cov.txt"
    txts = sorted([t for t in files if t[-14:] == ".stats.cov.txt"])

    df = pd.DataFrame()
    for t in txts:
        df = df.append(txt2series(t), ignore_index=True)

    # Ordering
    ordered_index = list(df.columns.difference(COL2FLOAT+COL2INT)) \
                        + COL2INT \
                        + COL2FLOAT
    df = df[ordered_index]

    # Fix type
    df[COL2FLOAT] = df[COL2FLOAT].astype(float)
    df[COL2INT] = df[COL2INT].astype(int)
    
    # Write to output
    df.to_csv(out, index=False)
