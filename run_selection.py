#! /usr/bin/env python
from __future__ import print_function, division
from collections import defaultdict, OrderedDict
import os
import sys
import math
import concurrent.futures
import gzip
import pickle
import json
import time
#import numexpr
import array
from functools import partial
import re
import numpy as np
from IPython.display import display

if __name__ == "__main__":

    sys.path.append("/home/xyan13/Trijet/TrijetAna")

    from TrijetAna.filelists.filelists import filelist
    from TrijetAna.filelists.samplenames import samples
    import argparse
    parser = argparse.ArgumentParser(description="Make histograms for Trijet data")
    
    input_group = parser.add_mutually_exclusive_group()
    input_group.add_argument("--subsamples", "-d", type=str, help="List of subsamples to run (comma-separated)")
    input_group.add_argument("--allsamples", "-a", type=str, help="Run all subsamples of the given sample (comma-separated)")
    input_group.add_argument("--test", "-t", action="store_true", help="Run a small test job")
    parser.add_argument("--quicktest", "-q", action="store_true", help="Run a small test job on selected dataset")
    parser.add_argument("--year", "-y", type=str, help="Year: 2016, 2017, or 2018")
    parser.add_argument("--isMC", "-m", action="store_true", help="Set run over MC instead of collision data")
    args = parser.parse_args()
    
    samples2process = []

    if args.test:
        year = "2017"
        samples2process = ["Res1ToRes2GluTo3Glu_M1-3000_R-0p5"]
        isMC = True
    elif args.allsamples:
        year = args.year
        allsamples = args.allsamples.split(",")
        isMC = args.isMC
        for item in allsamples:
            if "QCD" in item:
                samples2process += samples[year]["QCD"]
            if "SingleMuon" in item:
                samples2process += samples[year]["SingleMuon"]
            if "Res1ToRes2GluTo3Glu" in item:               
                samples2process += samples[year]["Res1ToRes2GluTo3Glu"]
            if "Res1ToRes2QTo3Q" in item:
                samples2process += samples[year]["Res1ToRes2QTo3Q"]
            if "ZprimeTo3Gluon" in item:
                samples2process += samples[year]["ZprimeTo3Gluon"]
            break
    elif args.subsamples:
        year = args.year
        samples2process = args.subsamples.split(",")
        isMC = args.isMC
    
    print("Please check samples to process: ", samples2process)

    # Make dictionary of subsample : [files to run]
    subsample_files = {}
    for subsample_name in samples2process:
        if not subsample_name in filelist[year]:
            raise ValueError(f"Dataset {subsample_name} not in dictionary.")
        # Drop some abundant QCD MC files
        if "QCD_Pt_600to800" in subsample_name:
            subsample_files[subsample_name] = filelist[year][subsample_name][:13]
        elif "QCD_Pt_800to1000" in subsample_name:
            subsample_files[subsample_name] = filelist[year][subsample_name][:3]
        elif "QCD_Pt_300to470" in subsample_name:
            subsample_files[subsample_name] = filelist[year][subsample_name]
        elif "QCD_Pt_470to600" in subsample_name:
            subsample_files[subsample_name] = filelist[year][subsample_name]
        elif "QCD_Pt_" in subsample_name:
            subsample_files[subsample_name] = filelist[year][subsample_name][:1]
        else:
            subsample_files[subsample_name] = filelist[year][subsample_name]

        if args.quicktest or args.test:
            subsample_files[subsample_name] = subsample_files[subsample_name][:1]

    for key in subsample_files.keys():
        size = len(subsample_files[key])
        print(f"For {key}, {size} file(s) will be processed")
        
    for sample, files in subsample_files.items():
        nEvents = -1
        if "QCD_Pt_1800to2400" in sample:
            nEvents = 15000
        if "QCD_Pt_2400to3200" in sample:
            nEvents = 1500
        if "QCD_Pt_3200toInf" in sample:
            nEvents = 500
        print("Processing QCD samples from which we only need a small amount of events!")
        print("Events to be processed: ", nEvents)
        
        with open(f"{sample}.txt", 'w') as f:
            for item in files:
                f.write(f"{item}\n")
        f.close()
        
        R1M_str = sample.split("_")[1].split("-")[1]
        rho_str = sample.split("_")[2].split("p")[1]
        rho = float(rho_str) * 0.1
        os.system(r"root -l -q -b -x selection/select_ML.C+\(\"" + sample + r"\"," + R1M_str + "," + str(rho) + r",\"outputs\"," + str(nEvents) + "," + str(year) + "\)")
        os.system(f"rm {sample}.txt")