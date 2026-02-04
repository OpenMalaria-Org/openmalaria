#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This file is part of OpenMalaria.

Copyright (C) 2005-2026 Swiss Tropical and Public Health Institute
Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
Copyright (C) 2020-2026 University of Basel
Copyright (C) 2026 The Kids Research Institute Australia

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

import os

#excludedir = ["..\\Lib"]

oldheaders=set()
nonheaders=set()

def update_source(filename, copyright):
    global oldheaders
    global nonheaders
    
    utfstr = chr(0xef)+chr(0xbb)+chr(0xbf)
    fdata = open(filename,"r+").read()
    isUTF = False
    if (fdata.startswith(utfstr)):
        isUTF = True
        fdata = fdata[3:]
    
    i=0
    line_comment=False
    block_comment=False
    prev_slash=False # or prev_star
    want_block=True # false once we have a block comment or set of line comments
    for c in fdata:
        i+=1
        if block_comment:
            # use prev_slash as prev_star instead
            if prev_slash and c=='/':
                block_comment=False
                prev_slash=False
                want_block=False
                # or: break to stop now (don't eat any new-lines)
            elif c=='*':
                prev_slash=True
            else:
                prev_slash=False
        elif c in [' ','\t','\f','\v']:
            pass
        elif c in ['\r','\n']:
            line_comment=False # end of line
        elif want_block and c=='/':
            if prev_slash:
                line_comment=True
                prev_slash=False
            else:
                prev_slash=True
        elif want_block and c=='*' and prev_slash:
            prev_slash=False
            block_comment=True
        else:
            i -= 1  # go back (keep this character)
            break   # end of header
    header=fdata[0:i]
    
    if header in oldheaders:
        fdata = fdata[len(header):]
    elif header in nonheaders:
        pass
    else:
        print(header)
        r=input("Remove above header? (y/N): ")
        if r[0]=='y' or r[0]=='Y':
            oldheaders.add(header)
            fdata = fdata[len(header):]
        else:
            nonheaders.add(header)
    
    if not (fdata.startswith(copyright)):
        print(("updating "+filename))
        fdata = copyright + fdata
        if (isUTF):
            open(filename,"w").write(utfstr+fdata)
        else:
            open(filename,"w").write(fdata)

def recursive_traversal(dir, copyright):
    global excludedir
    fns = os.listdir(dir)
    #print "listing "+dir
    for fn in fns:
        fullfn = os.path.join(dir,fn)
        #if (fullfn in excludedir):
            #continue
        if (os.path.isdir(fullfn)):
            recursive_traversal(fullfn, copyright)
        else:
            if (fullfn.endswith(".h") or fullfn.endswith(".cpp")):
                update_source(fullfn, copyright)


cright = open("util/licence-template.txt","r+").read()
recursive_traversal("model", cright)
recursive_traversal("unittest", cright)

print('Remember to update the text for --version in model/CommandLine.cpp!')

exit()
