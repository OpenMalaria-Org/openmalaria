#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of OpenMalaria.

Copyright (C) 2005-2014 Swiss Tropical Institute

OpenMalaria is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
"""

import sys
import math
from optparse import OptionParser
from approxEqual import ApproxSame
from readOutput import readEntries

REL_PRECISION=1e-6
ABS_PRECISION=1e-6

def charEqual (fn1,fn2):
    MAX=10*1024
    f1 = open(fn1,'r')
    f2 = open(fn2,'r')
    while True:
        s1 = f1.read(MAX)
        s2 = f2.read(MAX)
        if (len(s1)==0) or (len(s2)==0):
            # end of one or both files; equal if it's the end of both
            return len(s1) == len(s2)
        if s1 != s2:
            return False

def main(fn1,fn2,maxDiffsToPrint=6):
    """Takes names of the two files to compare and optionally an argument describing
the maximum number of differences to print directly (note: order is not intuitive).
Returns a tuple ret,ident; ret is 0 if test passes (output considered near-enough equal),
ident is 1 if files are binary-equal."""
    ret=0
    opt=""
    if REL_PRECISION!=1e-6:
        opt+=" --rel-prescision="+str(REL_PRECISION)
    if ABS_PRECISION!=1e-6:
        opt+=" --abs-prescision="+str(ABS_PRECISION)
    print(("\033[1;34m  compareOutput.py"+opt+" "+fn1+" "+fn2+" "+str(maxDiffsToPrint)+"\033[0;0m"))
    
    # Read both files and combine into a map of key to pairs (v1, v2)
    try:
        if charEqual (fn1,fn2):
            print("output.txt files are identical")
            return 0,True
        print("output.txt files aren't binary-equal")
        
        values1=readEntries(fn1)
        values2=readEntries(fn2)
    # python 3000 syntax is "except IOError as e", backported to 2.6 but not always supported. Old syntax:
    except IOError as e:
        print((str(e)))
        return 1,False
    values=dict()
    for (k,v1) in list(values1.items()):
        v2=None
        if (k in values2):
            v2=values2[k]
            del values2[k]
        values[k] = (v1,v2)
    for (k,v2) in list(values2.items()):
        values[k] = (None,v2)
    
    # Go through all values:
    numPrinted=0
    numDiffs=0
    numMissing1=0
    numMissing2=0
    perMeasureNum = dict()
    perMeasureTotal1 = dict()
    perMeasureTotal2 = dict()
    perMeasureNumDiff = dict()
    perMeasureDiffSum = dict()
    perMeasureDiffAbsSum = dict()
    approxSame = ApproxSame(REL_PRECISION, ABS_PRECISION)
    
    for (k,(v1,v2)) in list(values.items()):
        if v1==None:
            numMissing1 += 1
        elif v2==None:
            numMissing2 += 1
        else:
            perMeasureNum[k.a] = perMeasureNum.get(k.a, 0) + 1
            perMeasureTotal1[k.a] = perMeasureTotal1.get(k.a, 0.0) + v1
            perMeasureTotal2[k.a] = perMeasureTotal2.get(k.a, 0.0) + v2
            
            # Compare with relative precision
            if approxSame (v1, v2):
                continue
            
            numDiffs += 1
            # Sum up total difference per measure
            perMeasureDiffSum[k.a]    = perMeasureDiffSum.get(k.a,0.0)    + v2 - v1
            perMeasureDiffAbsSum[k.a] = perMeasureDiffAbsSum.get(k.a,0.0) + math.fabs(v2-v1)
        
        numPrinted += 1
        perMeasureNumDiff[k.a] = perMeasureNumDiff.get(k.a,0) + 1;
        if (numPrinted <= maxDiffsToPrint):
            print(("survey "+str(k.b)+", group "+str(k.c)+", measure "+str(k.a)+": "+str(v1)+" -> "+str(v2)))
            if (numPrinted == maxDiffsToPrint):
                print("[won't print any more line-by-line diffs]")
    
    if (numMissing1 > 0) or (numMissing2 > 0):
        print((str(numMissing1) + " entries missing from first file, " + str(numMissing2) +" from second"))
        ret = 3
    
    maxDiffSum=0.0
    maxAbsDiffSum=0.0
    for (k.a,absDiff) in list(perMeasureDiffAbsSum.items()):
        if not (absDiff <= 1e-6):   # handle NANs
            # standard division throws on divide-by-zero, which I don't want
            def div(x,y):
                try:
                    return x/y
                except ZeroDivisionError:
                    return 1e400 * 0 # nan
            
            diff=perMeasureDiffSum[k.a]
            sum1=perMeasureTotal1[k.a]
            sum2=perMeasureTotal2[k.a]
            diffSum=div(diff,sum1)
            maxDiffSum=max(maxDiffSum,math.fabs(diffSum))
            absDiffSum=div(absDiff,sum1)
            maxAbsDiffSum=max(maxAbsDiffSum,absDiffSum)
            print(("for measure "+str(k.a)+":\tsum(1st file):"+str(sum1)+"\tsum(2nd file):"+str(sum2)+"\tdiff/sum: "+str(diffSum)+"\t(abs diff)/sum: "+str(absDiffSum)))
    if maxDiffSum>0 or maxAbsDiffSum>0:
        print(("Max diff/sum:",maxDiffSum,"max (abs diff)/sum:",maxAbsDiffSum))
    
    if numDiffs == 0:
        print(("No significant differences (total relative diff: "+str(approxSame.getTotalRelDiff())+"), ok."))
        return ret,False
    else:
        print(("\033[1;31m"+str(numDiffs)+" significant differences (total relative diff: "+str(approxSame.getTotalRelDiff())+ ")!\033[0;0m"))
        return 1,False

# Test for options
def evalOptions (args):
    parser = OptionParser(usage="Usage: %prog [options] logfile1 logfile2 [max different lines to print]",
            # damn reformatting into a single paragraph: this doesn't get printed very nicely when --help is invoked
            description="""Compare logfile1 and logfile2 for differences, returning a measure of difference.
See https://github.com/OpenMalaria-Org/openmalaria/wiki/UtilsRunScripts#compareoutputpy for details on output.""")
    
    parser.add_option("-R","--rel-precision",
            action="store", dest="rel_precision", type="float",
            help="Set relative precision (default: 1.0e-6)")
    parser.add_option("-A","--abs-precision",
            action="store", dest="abs_precision", type="float",
            help="Set absolute precision (default: 1.0e-6)")
    (options, others) = parser.parse_args(args=args)
    
    return options,others

if __name__ == '__main__':
    (options,others) = evalOptions (sys.argv[1:])
    if options.rel_precision:
        REL_PRECISION=options.rel_precision
    if options.abs_precision:
        ABS_PRECISION=options.abs_precision
    
    if (len(others) == 3):
        ret,ident = main (others[0],others[1],int(others[2]))
    elif (len(others) == 2):
        ret,ident = main (others[0],others[1])
    else:
        print(("Usage: "+sys.argv[0]+" logfile1 logfile2 [max different lines to print]"))
        ret=-1
    sys.exit(ret)
