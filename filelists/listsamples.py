import os
from pprint import pprint

samples = {
	"2016": {
		"ZprimeTo3Gluon": ["ZprimeTo3Gluon_TuneCP5_13TeV_pythia8"]
	},
	"2017": {
		"JetHT"         : [f"JetHT_Run2017{x}" for x in ["B", "C", "D", "E", "F"]], 
		"SingleMuon"    : [f"SingleMuon_Run2017{x}" for x in ["B", "C", "D", "E", "F"]], 
		"QCD"           : ["QCD_Pt_300to470_TuneCP5_13TeV_pythia8", "QCD_Pt_470to600_TuneCP5_13TeV_pythia8", "QCD_Pt_600to800_TuneCP5_13TeV_pythia8", "QCD_Pt_800to1000_TuneCP5_13TeV_pythia8", "QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8", "QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8", "QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8", "QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8", "QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8"],
		"ZprimeTo3Gluon": ["ZprimeTo3Gluon_TuneCP5_13TeV_pythia8"]
	},
	"2018": {
		"ZprimeTo3Gluon": ["ZprimeTo3Gluon_TuneCP5_13TeV_pythia8"]
	}
}

prefixdir = {
	"2016": {
		"ZprimeTo3Gluon": "/isilon/hadoop/store/user/dryu/DijetSkim/v2_0_10/ZprimeTo3Gluon_2016"
	},
	"2017": {
		"JetHT"         : "/isilon/hadoop/store/user/dryu/DijetSkim/v2_0_7/JetHT2017/JetHT", 
		"SingleMuon"    : "/isilon/hadoop/store/user/dryu/DijetSkim/v2_0_7/SingleMuon2017/SingleMuon", 
		"QCD"           : "/isilon/hadoop/store/user/dryu/DijetSkim/v2_0_7/QCD_Pt_2017",
		"ZprimeTo3Gluon": "/isilon/hadoop/store/user/dryu/DijetSkim/v2_0_10/ZprimeTo3Gluon_2017"
	},
	"2018": {
		"ZprimeTo3Gluon": "/isilon/hadoop/store/user/dryu/DijetSkim/v2_0_10/ZprimeTo3Gluon_2018"
	}
}

if __name__ == "__main__":
	for year in ["2016","2017","2018"]:
		print(year)
		for sample, subsamples in samples[year].items():
			for subsample in subsamples:
				filelist = open(f"{year}_{subsample}.txt", "w")
				path = (prefixdir[year][sample])
				for path, subdirs, files in os.walk(f"{path}/{subsample}"):
					for filename in files:
						if ".root" not in filename or "hist" in filename:
							continue
						filelist.write(f"{os.path.join(path, filename)}\n")
				filelist.close()