#!/usr/bin/env python3
#
# This tool evolves the cgyro restart file
# Given an original (input.gen,restart) pair
# and a new input.gen,it generates a new restart file
# that keeps the old info
#
# Currently limited to adding an additional species
#

import sys,os
import argparse
import libcgyrorestart

def get_arguments():

   parser=argparse.ArgumentParser(description="Utility get info of a CGYRO restart file")

   parser.add_argument('-d',
                       metavar='ORIG',
                       help="Directory containing bin.cgyro.restart",
                       type=str,
                       required=True)
    
   args=parser.parse_args()

   return args.d


d=get_arguments()
head=libcgyrorestart.CGyroRestartHeader()
head.load(d, allow_old=True)
print(head)
