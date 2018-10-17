#!/usr/bin/python

import os, sys, fnmatch
import math

samples = [
# '2tevsig','3tevsig','4tevsig','5tevsig','6tevsig',
# 'ttbar0','ttbar1','ttbar2','ttbar3',
'ttinc0','ttinc1','ttinc2','ttinc3','ttinc4','ttinc5','ttinc6',
# 'qcd0','qcd1','qcd2','qcd3','qcd4',
# 'qcdpt50t80','qcdpt80t120','qcdpt120t170','qcdpt170t300','qcdpt300t470',
# 'qcdpt470t600','qcdpt600t800','qcdpt800t1000','qcdpt1000tInf'
]

# samples = [
# 'qcdflat',
# # 'qcdpt300t470f',
# ]

for sample in samples:
	outFile = open('output/out.'+sample, 'rU')
	outLines = outFile.readlines()
	outFile.close()

	jobStatus = {}
	for line in outLines:
		if sample+'_' in line:
			fileIDind = 1
			if len(line.strip().split())==3: fileIDind = 0
			fileID = line.strip().split()[fileIDind]
			if sample not in fileID: print line, fileID
			jobStatus[fileID] = 'child_success' in line

	failedFiles = []
	succesFiles = []
	for job in sorted(jobStatus.keys()):
		if jobStatus[job]: succesFiles.append(int(job.replace(sample+'_','')))
		if not jobStatus[job]: 
			failedFiles.append(int(job.replace(sample+'_','')))
			if os.path.exists(os.getcwd()+'output/'+job.replace('_','')+'_ntuple.root'):
				os.system('rm output/'+job.replace('_','')+'_ntuple.root')
	#for fileId in failedFiles: print fileId+17

	outConfFile = open('config/Zprime14TeV_'+sample+'.txt', 'rU')
	outConfLines = outConfFile.readlines()
	outConfFile.close()

	outConfFileNew=open('config/Zprime14TeV_'+sample+'fc.txt', 'w')
	count = 0
	for line in outConfLines:
		count+=1
		if count<17 or count-17 in failedFiles or count-17 not in succesFiles: 
			outConfFileNew.write(line.replace(sample,sample+'fc'))
