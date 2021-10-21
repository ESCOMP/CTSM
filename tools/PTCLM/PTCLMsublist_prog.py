#########################################################################################
#
# PTCLMsublist_prog
#
# Top level class to define the PTCLMsublist program. Parse's arguments and has Init,
# and run methods to submit a list of PTCLMmkdata sites to the batch queue.
#
#########################################################################################
import os, sys
from batchque import batchque

class PTCLMsublist_prog:
#----------------------------------------------------------------------------------------
# Class to handle command line input to the program
#----------------------------------------------------------------------------------------
   # Class data
   name         = "PTCLMsublist"
   cmdline      = ""
   account      = "P93300606"
   ctsmdir_def  = "../../"
   ctsmdir      = os.getenv("CTSM_ROOT", ctsmdir_def )
   inputdir_def = "/glade/p/cesmdata/cseg/inputdata"
   inputdir     = os.getenv("DIN_LOC_ROOT", inputdir_def )
   sitelistcsv  = "US-CHATS,US-FPe,CA-Let,US-NR1,CA-Man,BR-Sa1,BR-Sa3"
   mach         = "cheyenne"
   wall         = "02:00:00"
   parse_args   = False
   ptclm_opts   = ""
   que          = batchque()
   setup        = False

   # --  Error function ---------------------------------
   def error( self, desc ):
       "error function to abort with a message"
       print( "ERROR("+self.name+"):: "+desc )
       sys.exit(100)

   def parse_cmdline_args( self ):
      "Parse the command line arguments for the PTCLM batch submission script"
      from optparse import OptionParser, OptionGroup

      for arg in sys.argv:
          self.cmdline = self.cmdline+arg+" "
      parser   = OptionParser( usage="%prog [options]" )
      options = OptionGroup( parser, "Options" )
      options.add_option("-r", "--ctsm_root", dest="ctsm_root", default=self.ctsmdir, \
                        help="Location of CTSM root directory (also set with CTSM_ROOT env variable)")
      options.add_option("-d", "--inputdir", dest="inputdir", default=self.inputdir, \
                        help="Location of CESM inputdata directory (also set with CSMDATA env variable)")
      options.add_option("-o", "--PTCLM_options", dest="options", default=self.ptclm_opts, \
                        help="PTCLM options to run with")
      options.add_option("-l", "--list", dest="sitelist", default=self.sitelistcsv, \
                        help="Comma seperated list of PTCLM sites to submit to batch")
      options.add_option("--account", dest="account", default=self.account, \
                        help="Account number to use for batch queue")
      options.add_option("--wall", dest="wall", default=self.wall, \
                        help="Wall clock time to submit in queue for")
      options.add_option("--mach", dest="mach", default=self.mach, \
                        help="Machine name to use for batch submital")
      parser.add_option_group(options)
      stdout    = os.popen("pwd")
      cwd       = os.path.abspath( stdout.read().rstrip( ) )
      tagvers = ""
      clog = open( cwd+"/ChangeLog", "r" );
      for line in clog:
         if ( line.find("Tag: ",0) == 0 ): 
            n = line.count("")
            tagvers = line[5:n-2]
            break

      clog.close
      versiongroup    = OptionGroup( parser, tagvers )
      parser.add_option_group(versiongroup)
      (options, args) = parser.parse_args()
      if len(args) != 0:
          parser.error("incorrect number of arguments")
   
      self.mach         = options.mach
      self.sitelistcsv  = options.sitelist
      self.account      = options.account
      self.options      = options.options
      self.ctsmdir      = options.ctsm_root
      self.inputdir     = options.inputdir
      # Initialize batch que object, will abort if bad machine
      self.que.Initialize( self, mach=self.mach, account=self.account )
      # Error checking
      if ( not os.path.isdir(self.ctsmdir) ):
         self.error( "CTSM_root directory does NOT exist: "+self.ctsmdir )
      if ( not os.path.isdir(self.inputdir) ):
         self.error( "CESM inputdata directory does NOT exist: "+self.inputdir )

      # Get site list from csv formatted string list
      if ( self.sitelistcsv.find( " " ) != -1 ):
         self.error( "Site list has white space in it, just use comma's to seperate sites: "+self.sitelistcsv )
      if ( self.sitelistcsv.find( ",," ) != -1 or self.sitelistcsv.endswith( "," ) or self.sitelistcsv.startswith( "," ) ):
         self.error( "Site list has empty site names, make sure comma's do not go after each other: "+self.sitelistcsv )
      self.sitelist = self.sitelistcsv.split( "," )

      # Flag that parsing was accomplished
      self.parse_args   = True

   def ctsm_root( self ):
      "Return the CTSM_ROOT directory"
      if ( not self.parse_args ):
         self.error( "parse_cmdline_args was NOT run first" )
      return( self.ctsmdir )

   def get_SiteList( self ):
      "Return the Site list"
      if ( not self.parse_args ):
         self.error( "parse_cmdline_args was NOT run first" )
      return( self.sitelist )

   def Initialize( self ):
      "Initialize the PTCLM batch submission"
      if ( not self.parse_args ):
         self.error( "parse_cmdline_args was NOT run first" )

      self.que.Initialize( self, self.mach, self.account )
      self.setup = True


   def Submit( self, site, submit=True ):
      "Submit the PTCLMmkdata job to the batch queue"
      if ( not self.setup ):
         self.error( "Initialize was NOT run first" )

      jobcommand = "./PTCLMmkdata --ctsm_root "+self.ctsmdir+" -s "+site+" -d "+self.inputdir+" "+self.options
      print( jobcommand );
      bsub = self.que.Submit( self, jobcommand, jobname="PTCLM_"+site, submit=submit, wall=self.wall )
      return( bsub )

#
# Unit testing for above classes
#
import unittest

class test_PTCLMsublist_prog(unittest.TestCase):

   def setUp( self ):
     self.prog = PTCLMsublist_prog()

   def test_badinit( self ):
     # Bad option will fail
     self.prog = PTCLMsublist_prog()
     sys.argv[1:] = [ "--zztop" ]
     self.assertRaises(SystemExit, self.prog.parse_cmdline_args )
     # Test that doing stuff before parse_args fails
     self.prog = PTCLMsublist_prog()
     self.assertRaises(SystemExit, self.prog.ctsm_root )
     self.assertRaises(SystemExit, self.prog.Initialize )
     self.assertRaises(SystemExit, self.prog.Submit, "US-UMB" )
     # Test that doing stuff after parse_args before Initialize fails
     self.prog = PTCLMsublist_prog()
     sys.argv[1:] = [ ]
     self.prog.parse_cmdline_args( )
     self.assertRaises(SystemExit, self.prog.Submit, "US-UMB" )
     # Test that a non existant directory for ctsm_root fails
     self.prog = PTCLMsublist_prog()
     sys.argv[1:] = [ "--ctsm_root", "zztop" ]
     self.assertRaises(SystemExit, self.prog.parse_cmdline_args )
     # Test that a non existant directory for inputdata fails
     self.prog = PTCLMsublist_prog()
     sys.argv[1:] = [ "-d", "inpzztop" ]
     self.assertRaises(SystemExit, self.prog.parse_cmdline_args )
     # Test that a bad site list fails
     self.prog = PTCLMsublist_prog()
     sys.argv[1:] = [ "-l", "thing thing2 thing3" ]
     self.assertRaises(SystemExit, self.prog.parse_cmdline_args )
     self.prog = PTCLMsublist_prog()
     sys.argv[1:] = [ "-l", "thing,thing2,,thing3" ]
     self.assertRaises(SystemExit, self.prog.parse_cmdline_args )
     self.prog = PTCLMsublist_prog()
     sys.argv[1:] = [ "-l", "thing,thing2,thing3," ]
     self.assertRaises(SystemExit, self.prog.parse_cmdline_args )
     self.prog = PTCLMsublist_prog()
     sys.argv[1:] = [ "-l", ",thing,thing2,thing3" ]
     self.assertRaises(SystemExit, self.prog.parse_cmdline_args )

   def test_init( self ):
     # check that setting ctsm_root works
     sys.argv[1:] = [ ]
     self.prog.parse_cmdline_args( )
     ctsmdir_def = os.getenv("CTSM_ROOT", self.prog.ctsmdir_def )
     self.assertTrue( self.prog.ctsm_root( ) == ctsmdir_def )
     cwd = os.getcwd()
     sys.argv[1:] = [ "--ctsm_root", cwd ]
     self.prog.parse_cmdline_args( )
     self.assertTrue( self.prog.ctsm_root( ) == cwd )
     # Initialize and submit
     self.prog.Initialize( )
     site = "US-UMB"
     bsub = self.prog.Submit( site, submit=False )
     jobid = str(os.getpid())
     checkstring = "qsub  -o PTCLM_US-UMB."+jobid+".stdout.out  -N PTCLM_US-UMB  -l walltime=02:00:00  " + \
                   "-A P93300606 -l select=3:ncpus=1:mpiprocs=1:mem=109GB -q regular " + \
                   "-V -m ae -j oe  PTCLM_"+site+"."+jobid+".job"
     print( "\n" )
     print( "bsubcm:"+bsub+":end" )
     print( "expect:"+checkstring+":end" )
     self.assertTrue( bsub == checkstring )

   def test_sitelist( self ):
     sitelistcsv = "US-UMB,US-Ha1"
     sitelist    = [ "US-UMB", "US-Ha1" ]
     sys.argv[1:] = [ "-l", sitelistcsv ]
     self.prog.parse_cmdline_args( )
     self.prog.Initialize( )
     slist = self.prog.get_SiteList( )
     self.assertTrue( slist == sitelist )
       

if __name__ == '__main__':
     unittest.main()
