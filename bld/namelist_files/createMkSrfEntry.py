#!/usr/bin/env python

import os, sys

class mksrfDataEntry_prog:

   # Class data
   year_start = 2016
   year_end   = 2100
   ssp_rcp    = "SSP5-8.5"
   subdir     = "pftcftdynharv.0.25x0.25.SSP5-8.5.simyr2016-2100.c171005"
   cdate      = 171005
   desc       = "SSP5RCP85_clm5"

   def parse_cmdline_args( self ):
      "Parse the command line arguments for create data entry list"
      from optparse import OptionParser, OptionGroup

      parser   = OptionParser( usage="%prog [options]" )
      options = OptionGroup( parser, "Options" )
      options.add_option( "-s", "--year_start", dest="year_start", default=self.year_start, \
                        help="Start year" )
      options.add_option( "-f", "--year_end", dest="year_end", default=self.year_end, \
                        help="End year" )
      options.add_option( "-d", "--subdir", dest="subdir", default=self.subdir, \
                        help="Subdirectory" )
      options.add_option( "--cdate", dest="cdate", default=self.cdate, \
                        help="Creation date" )
      options.add_option( "--desc", dest="desc", default=self.desc, \
                        help="Description string" )
      parser.add_option_group(options)
      (options, args) = parser.parse_args()
      if len(args) != 0:
          parser.error("incorrect number of arguments")

      self.year_start = options.year_start
      self.year_end   = options.year_end
      self.subdir     = options.subdir
      self.cdate      = options.cdate
      self.desc       = options.desc

   def printentry( self, year ):
      "Print a single entry"
      print '<mksrf_fvegtyp hgrid="0.25x0.25" ssp_rcp="%s" sim_year="%d" crop="on"' % (self.ssp_rcp, year)
      print '>lnd/clm2/rawdata/%s/mksrf_landuse_%s_%s.c%s.nc' % (self.subdir, self.desc, year, self.cdate)
      print '</mksrf_fvegtyp>\n'

entry = mksrfDataEntry_prog()
entry.parse_cmdline_args()

for year in range(entry.year_start, entry.year_end+1):
   entry.printentry( year )




