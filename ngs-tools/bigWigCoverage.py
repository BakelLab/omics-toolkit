#! /usr/bin/env python

# IMPORT
import sys
import os
import pyBigWig as py
import csv
import pandas as pd
import numpy as np

#############
# FUNCTIONS #
#############

## Parse options
def parse_options(argv):
   from optparse import OptionParser, OptionGroup
   parser   = OptionParser()
   required = OptionGroup(parser, 'MANDATORY')
   required.add_option('-r', '--region_file', dest='bed_filename',    metavar='FILE.bed', help='Region list in bed6 format', default='-')
   required.add_option('-p', '--bw_pos_file', dest='bw_pos_filename', metavar='FILE.bw',  help='Bigwig file with postive strand coverage', default='-')
   required.add_option('-n', '--bw_neg_file', dest='bw_neg_filename', metavar='FILE.bw',  help='Bigwig file with negative strand coverage', default='-')
   
   parser.add_option_group(required)
   (options, args) = parser.parse_args()

   if len(argv) < 6:
      parser.print_help()
      sys.exit(2)

   options.parser = parser
   return options


########
# MAIN #
########
 
def main():
   # get command line options
   opt = parse_options(sys.argv)
   
   # Test if all files are accessible and readable
   d = vars(opt)
   for key in ["bed_filename", "bw_pos_filename", "bw_neg_filename"]:
      if not ( os.path.isfile(d[key]) and os.access(d[key], os.R_OK) ):
         print >> sys.stderr, 'ERROR: file "%s" does not exist or is not readable.' % d[key]
         sys.exit(1)

   # Open positive and negative strand bigwig files
   bw_pos = py.open( opt.bw_pos_filename )
   bw_neg = py.open( opt.bw_neg_filename )

   # Open bed file with regions and output count data per region position
   range_data = []
   df = pd.read_csv( opt.bed_filename, header=None, names=["chr","start","stop","name","score","strand"], sep="\t")
   for row in df.itertuples():
      if row[6] == "+":
         range_data = bw_pos.values(row[1], row[2], row[3])
      else:
         range_data = list(reversed( bw_neg.values(row[1], row[2], row[3]) ) )
      range_data = np.nan_to_num(range_data)
      for i in range(0, len(range_data) ):
         print "\t".join([str(row[4]), str(i), str(range_data[i])])
      
      
if __name__ == "__main__":
   main()
