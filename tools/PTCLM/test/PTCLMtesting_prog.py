#########################################################################################
#
# PTCLMtesting_prog
#
# Top level class to define the PTCLMtesting program. Parse's arguments and has Init,
# run, and final methods to setup, execute, and then display test results.
#
#########################################################################################
from PTCLMtestlist import PTCLMtestlist
import os, sys

class PTCLMtesting_prog:
#----------------------------------------------------------------------------------------
# Class to handle command line input to the program
#----------------------------------------------------------------------------------------
   # Class data
   name         = "PTCLMtesting"
   cmdline      = ""
   redo_compare = False
   ctsmdir_def  = "../../../"
   ctsmdir      = os.getenv("CTSM_ROOT", ctsmdir_def )
   cimedir_def  = ctsmdir_def + "/cime"
   cimedir      = os.getenv("CIME_ROOT", cimedir_def )
   parse_args   = False
   setup        = False
   n_tests      = {}
   testlist     = PTCLMtestlist()

   # --  Error function ---------------------------------
   def error( self, desc ):
       "error function to abort with a message"
       print( "ERROR("+self.name+"):: "+desc )
       sys.exit(100)

   def parse_cmdline_args( self ):
      "Parse the command line arguments for the PTCLM testing script"
      from optparse import OptionParser, OptionGroup

      for arg in sys.argv:
          self.cmdline = self.cmdline+arg+" "
      parser   = OptionParser( usage="%prog [options]" )
      options = OptionGroup( parser, "Options" )
      options.add_option("-r", "--ctsm_root", dest="ctsm_root", default=self.ctsmdir, \
                        help="Location of CTSM root directory (also set with CTSM_ROOT env variable)")
      options.add_option("--cime_root", dest="cime_root", default=self.cimedir, \
                        help="Location of CIME root directory (also set with CIME_ROOT env variable)")
      options.add_option("--redo_compare_files", dest="redo_compare", action="store_true", default=self.redo_compare, \
                        help="Redo the compare files")
      parser.add_option_group(options)
      tagvers = ""
      clog = open( "../ChangeLog", "r" );
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
   
      self.redo_compare = options.redo_compare
      self.ctsmdir      = options.ctsm_root
      self.cimedir      = options.cime_root
      self.parse_args   = True

   def redo_compare_files( self ):
      "Check if should redo the comparison files"

      if ( not self.parse_args ):
         self.error( "parse_cmdline_args was NOT run first" )
      return( self.redo_compare )

   def ctsm_root( self ):
      "Return the CTSM_ROOT directory"
      if ( not self.parse_args ):
         self.error( "parse_cmdline_args was NOT run first" )
      return( self.ctsmdir )

   def cime_root( self ):
      "Return the CIME_ROOT directory"
      if ( not self.parse_args ):
         self.error( "parse_cmdline_args was NOT run first" )
      return( self.cimedir )

   def Initialize( self ):
      "Initialize the PTCLM tests"
      if ( not self.parse_args ):
         self.error( "parse_cmdline_args was NOT run first" )
      self.testlist.Setup(self.ctsm_root())
      self.testlist.Read()


      self.n_tests['PASS']          = 0;
      self.n_tests['FAIL']          = 0;

      self.n_tests['PASS-COMP']     = 0;
      self.n_tests['FAIL-COMP']     = 0;
      self.n_tests['NO-COMPS-DONE'] = 0;
      self.itest = 0
      self.setup = True

   def Run( self ):
      "Run the testing"
      if ( not self.setup ):
         self.error( "Initialize was NOT run first" )
      for test in self.testlist.get_testlist():
         self.itest += 1
         (stat,comp) = self.testlist.run_PTCLMtest( test, self.redo_compare_files() )
         self.n_tests[stat] += 1
         if   ( comp == "PASS" ):          self.n_tests['PASS-COMP']  += 1
         elif ( comp == "NO-COMPS-DONE" ): self.n_tests[comp]         += 1
         else:                             self.n_tests['FAIL-COMP']  += 1

      for test in self.testlist.get_failtestlist():
         self.itest += 1
         (stat,comp) = self.testlist.run_PTCLMtest( test, self.redo_compare_files() )
         self.n_tests[stat] += 1
         if   ( comp == "PASS" ):          self.n_tests['PASS-COMP']  += 1
         elif ( comp == "NO-COMPS-DONE" ): self.n_tests[comp]         += 1
         else:                             self.n_tests['FAIL-COMP']  += 1


   def Finalize( self ):
      "Finalize and print out results"
      if ( not self.setup ):
         self.error( "Initialize was NOT run first" )
      print( "Total number of tests             = "+str(self.itest) )
      print( "Number of tests that PASS         = "+str(self.n_tests['PASS']) )
      print( "Number of tests that Fail         = "+str(self.n_tests['FAIL']) )
      print( "Number of compare tests that PASS = "+str(self.n_tests['PASS-COMP']) )
      print( "Number of compare tests that Fail = "+str(self.n_tests['FAIL-COMP']) )
      print( "Number of tests without compare   = "+str(self.n_tests['NO-COMPS-DONE']) )


#
# Unit testing for above classes
#
import unittest

class test_PTCLMtesting_prog(unittest.TestCase):

   def setUp( self ):
     self.prog = PTCLMtesting_prog()

   def test_badinit( self ):
     # Bad option will fail
     sys.argv[1:] = [ "--zztop" ]
     self.assertRaises(SystemExit, self.prog.parse_cmdline_args )
     # Test that doing stuff before parse_args fails
     self.prog = PTCLMtesting_prog()
     self.assertRaises(SystemExit, self.prog.redo_compare_files )
     self.assertRaises(SystemExit, self.prog.ctsm_root )
     self.assertRaises(SystemExit, self.prog.cime_root )
     self.assertRaises(SystemExit, self.prog.Initialize )
     self.assertRaises(SystemExit, self.prog.Run )
     self.assertRaises(SystemExit, self.prog.Finalize )
     # Test that doing stuff after parse_args before Initialize fails
     sys.argv[1:] = [ ]
     self.prog.parse_cmdline_args( )
     self.assertRaises(SystemExit, self.prog.Run )
     self.assertRaises(SystemExit, self.prog.Finalize )

   def test_init( self ):
     # check that setting redo_compare_files works
     sys.argv[1:] = [ ]
     self.prog.parse_cmdline_args( )
     self.assertFalse( self.prog.redo_compare_files( ) )
     sys.argv[1:] = [ "--redo_compare_files" ]
     self.prog.parse_cmdline_args( )
     self.assertTrue( self.prog.redo_compare_files( ) )
     # check that setting ctsm_root works
     sys.argv[1:] = [ ]
     self.prog.parse_cmdline_args( )
     ctsmdir_def = os.getenv("CTSM_ROOT", self.prog.ctsmdir_def )
     self.assertTrue( self.prog.ctsm_root( ) == ctsmdir_def )
     checkstring  = "A_string_to_check_to_make_sure_this_works"
     sys.argv[1:] = [ "--ctsm_root", checkstring ]
     self.prog.parse_cmdline_args( )
     self.assertTrue( self.prog.ctsm_root( ) == checkstring )
     # check that setting cime_root works
     sys.argv[1:] = [ ]
     self.prog.parse_cmdline_args( )
     cimedir_def = os.getenv("CIME_ROOT", self.prog.cimedir_def )
     self.assertTrue( self.prog.cime_root( ) == cimedir_def )
     checkstring  = "A_string_to_check_to_make_sure_setting_cimeroot_works"
     sys.argv[1:] = [ "--cime_root", checkstring ]
     self.prog.parse_cmdline_args( )
     self.assertTrue( self.prog.cime_root( ) == checkstring )

if __name__ == '__main__':
     unittest.main()
