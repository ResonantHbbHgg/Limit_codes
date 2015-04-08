#!/usr/bin/env python

#run with args <output from combine> <nonres senario> <output file for all limits>
#the goal is to extract the limits (6 numbers) from input #1 and the scenario (lambda, yt, c2) from input #2
#the output is appended to output #3

import sys
from string import atof

if len(sys.argv) < 4:
    print "Not enough arguments"
    exit(1)

    
combineFileName=sys.argv[1]
scenarioName=sys.argv[2]
outputFileName=sys.argv[3]

combineFile = open(combineFileName,'r')

obsLim = -1
expLim = [-1, -1, -1, -1, -1]

for line in combineFile:
    
    #Does the line contain Observerd?
    if line.find("Observed Limit:") >= 0 :
        obsLim=atof(line[line.find('r < ')+4:])
        continue

    #Does the line contain Expected?
    if line.find("Expected  2.5%:") >= 0 :
        expLim[0]=atof(line[line.find('r < ')+4:])
        continue
    if line.find("Expected 16.0%:") >= 0 :
        expLim[1]=atof(line[line.find('r < ')+4:])
        continue
    if line.find("Expected 50.0%:") >= 0 :
        expLim[2]=atof(line[line.find('r < ')+4:])
        continue
    if line.find("Expected 84.0%:") >= 0 :
        expLim[3]=atof(line[line.find('r < ')+4:])
        continue
    if line.find("Expected 97.5%:") >= 0 :
        expLim[4]=atof(line[line.find('r < ')+4:])
        continue

combineFile.close()

if obsLim < 0:
    obsLim = expLim[2]


#parse the scenario name
scenarioList=scenarioName.replace('d','.').replace('m','-').split('_')

if len(scenarioList) <6:
    print "Invalid scenario name"
    exit(1)

lambda_hhh=atof(scenarioList[1])
yt=atof(scenarioList[3])
c2=atof(scenarioList[5])


#append everything to outputFileName

outputfile = open(outputFileName,'a')
outputfile.write("%2d %4.2f %3d %.4f %.4f %.4f %.4f %.4f %.4f\n" % (c2,yt,lambda_hhh,expLim[0],expLim[1],expLim[2],expLim[3],expLim[4],obsLim) )
outputfile.close()

