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
import ROOT

def readfiles(subsample, year):
    filelist = []
    inputfile = open(f'filelists/{year}_{subsample}.txt', 'r')
    for line in inputfile.readlines():
        filelist.append(line.rstrip('\n'))
    return filelist

if __name__ == "__main__":

    sys.path.append("/home/xyan13/Trijet/TrijetAna")

    from TrijetAna.filelists.listsamples import samples
    import argparse
    parser = argparse.ArgumentParser(description="Make histograms for Trijet data")
    
    input_group = parser.add_mutually_exclusive_group()
    input_group.add_argument("--subsamples", "-d", type=str, help="List of subsamples to run (comma-separated)")
    input_group.add_argument("--allsamples", "-a", type=str, help="Run all subsamples of the given sample (comma-separated)")
    input_group.add_argument("--test", "-t", action="store_true", help="Run a small test job")
    parser.add_argument("--year", "-y", type=str, help="Year: 2016, 2017, or 2018")
    parser.add_argument("--isMC", "-m", action="store_true", help="Set run over MC instead of collision data")
    args = parser.parse_args()
    
    samples2process = []

    if args.test:
        year = "2017"
        samples2process += samples[year]["ZprimeTo3Gluon"]
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
        # Drop some abundant QCD MC files
        if "QCD_Pt_600to800" in subsample_name:
            subsample_files[subsample_name] = readfiles(subsample_name, args.year)[:13]
        elif "QCD_Pt_800to1000" in subsample_name:
            subsample_files[subsample_name] = readfiles(subsample_name, args.year)[:3]
        else:
            subsample_files[subsample_name] = readfiles(subsample_name, args.year)

        if args.test:
            subsample_files[subsample_name] = readfiles(subsample_name, args.year)[:1]

    for key in subsample_files.keys():
        size = len(subsample_files[key])
        print(f"For {key}, {size} file(s) will be processed")
    
    print(subsample_files)
    
    for sample, files in subsample_files.items():
        nEvents = -1
        # if "QCD_Pt_1800to2400" in sample:
            # nEvents = 15000
            # print("Processing QCD samples from which we only need a small amount of events!")
        # if "QCD_Pt_2400to3200" in sample:
            # nEvents = 1500
            # print("Processing QCD samples from which we only need a small amount of events!")
        # if "QCD_Pt_3200toInf" in sample:
            # nEvents = 500
            # print("Processing QCD samples from which we only need a small amount of events!")
        print("Events to be processed: ", nEvents)
        total_events = 0
        trig_events = 0
        sel_events = 0
        with open(f"{sample}.txt", 'w') as f:
            for item in files:
                f.write(f"{item}\n")
                fdir = item.split("nanoskim")[0]
                filei = item.split("nanoskim")[1]
                tfile = ROOT.TFile.Open(fdir+"hists"+filei, "READ")
                ntott = tfile.Get("h_ProcessedEvents")
                total_events += ntott.GetEntries()
                ntrigt = tfile.Get("h_TriggeredEvents")
                trig_events += ntrigt.GetEntries()
                nselt = tfile.Get("h_SelectedEvents")
                sel_events += nselt.GetEntries()
        f.close()
        print(sample, total_events, trig_events, sel_events)
        
        # if "QCD" in sample:
            # os.system(r"root -l -q -b -x selection/3_jets/select_ML_QCD.C+\(\"" + sample + r"\"," + r"\"outputs_temp\"," + str(nEvents) + "," + str(year) + "\) 2>&1 | tee " + f"{sample}.log")
        # else:
            # os.system(r"root -l -q -b -x selection/3_jets/select_ML_Fullmatch.C+\(\"" + sample + r"\"," + r"\"outputs_temp\"," + str(nEvents) + "," + str(year) + "\) 2>&1 | tee " + f"{sample}.log")
        # os.system(f"rm {sample}.txt")
        
        os.system(r"root -l -q -b -x selection/EffHLT/HLTEff.C+\(\"" + sample + r"\"," + r"\"outputs_temp\"," + str(nEvents) + "," + str(year) + "\) 2>&1 | tee " + f"{sample}.log")