#----------------------------------------------------------------------------------
# Class to read in a test list XML file for PTCLM
#----------------------------------------------------------------------------------
from   xml.sax.handler import ContentHandler
class PTCLMtestlistXML( ContentHandler ):

   def startDocument(self):
     self.testlist     = [];
     self.failtestlist = [];

   def startElement(self, name, attrs):

     attributes = [ 'id', 'type', 'site', 'opts', 'resultfile', 'compdir' ]
     testmap    = {}
     for key in attributes:
       testmap[key] = str( attrs.get(key,"") )

     if name == 'test':
       for test in self.testlist:
          if ( testmap['compdir'] != "" and testmap['compdir'] == test['compdir'] ):
             print( "compdir is repeated: "+test['compdir'] )
             sys.exit(100)
          if ( testmap['id'] == test['id'] and testmap['site'] == test['site']):
             print( "id and site is duplicated: "+test['id']+" "+test['site'] )
             sys.exit(100)

       self.testlist.append( testmap )

     if name == 'failtest':
       self.failtestlist.append( testmap )

#----------------------------------------------------------------------------------
# Class to read in a test list for PTCLM and do operations on it
#----------------------------------------------------------------------------------
from   xml.sax         import make_parser
import os
import sys
class PTCLMtestlist:
   # Class data
   testlist     = "PTCLMtestlist.xml"
   testing_dir  = "testing_dir"
   inputdatadir = "/glade/p/cesmdata/cseg/inputdata"
    
   # Construct the class
   def Setup( self, ctsmdir ):
     self.testXML = make_parser()
     self.xml     = PTCLMtestlistXML()
     self.testXML.setContentHandler(self.xml)
     # Get the CLM root directory
     self.ctsmdir = ctsmdir
     if ( not os.path.exists(self.ctsmdir) or not os.path.isdir(self.ctsmdir) or not os.path.exists(self.ctsmdir+"/doc/ChangeLog") ):
         self.error("CTSM_ROOT does NOT exist or NOT a directory (set with CTSM_ROOT env variable):"+self.ctsmdir)

   # Set the input-dir
   def SetDataDir( self, inputdatadir ):
     self.inputdatadir = inputdatadir
     if ( not os.path.exists(self.inputdatadir) or \
          not os.path.isdir(self.inputdatadir)    ):
         self.error("Inputdatadir does NOT exist or NOT a directory:"+self.inputdatadir)

   # Read in the testlist file
   def Read( self ):
     if ( not os.path.exists( self.testlist) ): self.error("File does NOT exist:"+self.testlist)
     print( "Open file: "+self.testlist )
     self.testXML.parse( self.testlist )

   # Replace specific identifyable information with generic
   def ReplaceIDInfoInFile( self, filename, curdates=False ):
     if ( not os.path.exists( filename ) ): self.error("File does NOT exist to do replace operation in:"+filename)
     outfilename = filename+".tmp"
     infile  = open( filename,    "r" );
     outfile = open( outfilename, "w" );
     stdout  = os.popen( "cd "+self.ctsmdir+"; pwd" )
     absctsmdir = stdout.read().rstrip( );
     stdout  = os.popen( "date +%y%m%d" );
     sdate   = stdout.read().rstrip( );
     stdout  = os.popen( "../PTCLMmkdata --version" )
     version = stdout.read().rstrip()
     stdout  = os.popen( "which ncl" )
     nclpath = stdout.read().rstrip( );
     for line in infile:
         line = line.replace( self.ctsmdir,      "$CTSMROOTDIR" );
         line = line.replace( absctsmdir,        "$CTSMROOTDIR" );
         line = line.replace( nclpath,           "$NCLPATH"     );
         line = line.replace( self.inputdatadir, "$DIN_LOC_ROOT" );
         line = line.replace( version,           "$PTCLMVERSION" );
         if ( curdates ):
            line = line.replace( sdate,     "$YYMMDD" );
         outfile.write( line )

     infile.close()
     outfile.close()
     os.system( "mv "+outfilename+" "+filename )
     
      

   # --  Error function ---------------------------------
   def error( self, desc ):
       "error function"
       print( "ERROR(PTCLMtestlist):: "+desc )
       sys.exit(100)

   # Get the list of tests to do
   def get_testlist( self ):
     return( self.xml.testlist )

   # Get the list of fail tests to do
   def get_failtestlist( self ):
     return( self.xml.failtestlist )

   # Get an Identifier name for this test
   def get_testID( self, test ):
     tid = ""
     for att in ["id","opts","site"]:
        val  = test[att].replace(" ","+")
        if ( val != "" ):
           tid  = tid + val + "."

     return(  tid )

   # Fail test or not...
   def IsNOTFailTest( self, test ):
     if ( test['type'] != "Fail" ):
        return( True  )
     else:
        return( False )
     
   # Get the PTCLM command line options to use
   def get_PTCLMoptions( self, test ):
     opts = test['opts']
     opts = opts + " --ctsm_root "+self.ctsmdir
     if ( self.IsNOTFailTest( test ) ):
        opts = opts + " -s " + test['site']
        opts = opts + " -d "+self.inputdatadir
        opts = opts + " --debug"
        opts = opts + " --sdate 140204"
        opts = opts + " --map_gdate 140204"
        if ( test['type'] == "RUN" ):
           opts = opts + " --mydatadir "+self.testing_dir+"/"+test['id']


     return( opts )

   # run the test
   def run_PTCLMtest( self, test, redo_comp_files ):
     opts = self.get_PTCLMoptions( test )
     tid  = self.get_testID( test )
     tlog = "run.log"
     if ( self.IsNOTFailTest( test ) ):
        errcode =  os.system( "../PTCLMmkdata "+opts+" > "+tlog  )
     else:
        errcode =  os.system( "../PTCLMmkdata "+opts )
     if ( errcode == 0 ):
        stat = True
     else:
        stat = False

     teststatus = ""
     if ( stat == self.IsNOTFailTest( test ) ):
        teststatus = "PASS"
     else:
        teststatus = "FAIL"

     print( teststatus+" "+tid )

     overallcompstatus = "NO-COMPS-DONE"
     if ( os.path.exists(tlog) ): self.ReplaceIDInfoInFile( tlog, curdates=True )   # Replace identifing info in the log file
     if ( test['resultfile'] != "" ):
        dstat = os.system( "diff "+tlog+" "+test['resultfile']  )
        if ( dstat == 0 ):
           compstatus = "PASS"
           desc       = " compare to result file"
        else:
           compstatus = "FAIL-COMP"
           desc       = " different from result file: "+test['resultfile']

        if ( redo_comp_files ):
           if ( not os.path.exists( test['resultfile'] ) ):
              cmd = "mkdir "+os.path.dirname(test['resultfile'])
              print( "Create new file directory that does NOT exist\n" )
              print( cmd+"\n" )
              os.system( cmd )
           os.system( "cp "+tlog+" "+test['resultfile']  )
        overallcompstatus = compstatus

        print( compstatus+" "+tid+" "+desc )

     if ( test['compdir'] != "" ):
        testdir  = os.path.abspath(self.testing_dir)+"/"+test['id']
        stdout   = os.popen("cd "+testdir+"/*"+test['site']+"; pwd")
        testdir  = os.path.abspath( stdout.read().rstrip( ) )
        os.system( "mv "+tlog+" "+testdir )
        filelist = ["README.PTCLM", tlog, "user_nl_clm", "shell_commands" ]
        for file in filelist:
           srcfile = testdir+"/"+file
           cmpfile = "compdirs/"+test['compdir']+"/"+file
           if ( not os.path.exists( srcfile ) ):
                 compstatus = "FAIL-DNE"
                 desc       = "source compare file does NOT exist: "+srcfile
           elif ( not os.path.exists( cmpfile ) ):
                 compstatus = "BFAIL"
                 desc       = "compare file does NOT exist: "+cmpfile
           else:
              self.ReplaceIDInfoInFile( srcfile )   # Replace identifing info in the log file
              dstat = os.system( "diff "+srcfile+" "+cmpfile )

              if ( dstat == 0 ):
                 compstatus = "PASS"
                 desc       = "same as comp directory file: "+cmpfile
              else:
                 compstatus = "FAIL-COMP"
                 desc       = "different from comp directory file: "+cmpfile

           if ( redo_comp_files ):
              if ( not os.path.exists( cmpfile ) ):
                 cmd = "mkdir "+os.path.dirname(cmpfile)
                 print( "Create new file directory that does NOT exist\n" )
                 print( cmd+"\n" )
                 os.system( cmd )
              os.system( "cp "+srcfile+" "+cmpfile )
                 
           if ( overallcompstatus == "PASS" or overallcompstatus == "NO-COMPS-DONE" ): overallcompstatus = compstatus

           print( compstatus+" "+tid+" "+desc+" "+file )

     if ( teststatus == "PASS" and overallcompstatus == "PASS" and os.path.exists( tlog ) ): 
        os.system( "/bin/rm "+tlog )

     return( [teststatus, overallcompstatus] )
#
# Unit testing for above classes
#
import unittest

class test_PTCLMtestlist(unittest.TestCase):

   def setUp( self ):
     self.test = PTCLMtestlist()
     ctsmdir = os.getenv("CTSM_ROOT", "../../../")
     if ( not os.path.exists( ctsmdir+"/doc/ChangeLog" ) ):
         print( "CTSM_ROOT NOT input\n" )
         sys.exit( 200 )

     self.test.Setup(ctsmdir)

   def test_read( self ):
     inpd_orig = self.test.inputdatadir
     self.assertRaises(SystemExit, self.test.SetDataDir, "nodir" )
     self.assertRaises(SystemExit, self.test.SetDataDir, "PTCLMtestlist.py" )
     self.test.Read()
     self.test.SetDataDir( inpd_orig )
     self.test.Read()
     print( "\nValid Tests: " )
     for test in self.test.get_testlist():
        print( str(test) )
        self.assertTrue( self.test.IsNOTFailTest( test ) )
     print( "\n\n" )
     print( "\nValid Fail Tests: " )
     for test in self.test.get_failtestlist():
        print( str(test) )
        self.assertTrue( not self.test.IsNOTFailTest( test ) )
     print( "\n\n" )

   def test_getopts( self ):
     self.test.Read()
     for test in self.test.get_testlist():
       print( self.test.get_testID( test )+" = "+ self.test.get_PTCLMoptions( test ) )

     for test in self.test.get_failtestlist():
       print( self.test.get_testID( test )+" = "+ self.test.get_PTCLMoptions( test ) )

   def test_run( self ):
     self.test.Read()
     testlist     = self.test.get_testlist()
     failtestlist = self.test.get_failtestlist()

     self.test.run_PTCLMtest( testlist[0],     False )
     self.test.run_PTCLMtest( testlist[5],     False )
     self.test.run_PTCLMtest( failtestlist[0], False )

   def test_ReplaceIDInfoInFile( self ):
     self.assertRaises(SystemExit, self.test.ReplaceIDInfoInFile, "non-existant-file" )
     filename = "testreplace.tmp"
     tfile = open( filename, "w" )
     tfile.write( "file = '"+self.test.ctsmdir+"/myfile'\n" )
     stdout     = os.popen( "cd "+self.test.ctsmdir+" ; pwd" )
     absctsmdir = stdout.read().rstrip()
     tfile.write( "file2 = '"+absctsmdir+"/myfile'\n" )
     tfile.write( "inpfile = '"+self.test.inputdatadir+"/myfile'\n" )
     tfile.write( "inpfile2 = '$DIN_LOC_ROOT/myfile'\n" )
     nclpath    = "/glade/u/apps/ch/opt/ncl/6.6.2/intel/18.0.5/bin/ncl"
     tfile.write( nclpath+"\n" )
     stdout        = os.popen( "date +%y%m%d" );
     sdate         = stdout.read().rstrip( );
     tfile.write( self.test.ctsmdir+"/surfdata_"+sdate+".namelist\n" )
     stdout  = os.popen( "../PTCLMmkdata --version" )
     version = stdout.read().rstrip()
     tfile.write( "++++ "+version+" +++++\n" )
     tfile.close()
     self.test.ReplaceIDInfoInFile( filename )
     tfile = open( filename, "r" )
     line = tfile.readline()
     self.assertEqual( line, "file = '$CTSMROOTDIR/myfile'\n" )
     line = tfile.readline()
     self.assertEqual( line, "file2 = '$CTSMROOTDIR/myfile'\n" )
     line = tfile.readline()
     self.assertEqual( line, "inpfile = '$DIN_LOC_ROOT/myfile'\n" )
     line = tfile.readline()
     self.assertEqual( line, "inpfile2 = '$DIN_LOC_ROOT/myfile'\n" )
     line = tfile.readline()
     self.assertEqual( line, "$NCLPATH\n" )
     line = tfile.readline()
     self.assertEqual( line, "$CTSMROOTDIR/surfdata_"+sdate+".namelist\n" )
     line = tfile.readline()
     self.assertEqual( line, "++++ $PTCLMVERSION +++++\n" )
     tfile.close()
     self.test.ReplaceIDInfoInFile( filename, curdates=True )
     tfile = open( filename, "r" )
     line = tfile.readline()
     self.assertEqual( line, "file = '$CTSMROOTDIR/myfile'\n" )
     line = tfile.readline()
     self.assertEqual( line, "file2 = '$CTSMROOTDIR/myfile'\n" )
     line = tfile.readline()
     self.assertEqual( line, "inpfile = '$DIN_LOC_ROOT/myfile'\n" )
     line = tfile.readline()
     self.assertEqual( line, "inpfile2 = '$DIN_LOC_ROOT/myfile'\n" )
     line = tfile.readline()
     self.assertEqual( line, "$NCLPATH\n" )
     line = tfile.readline()
     self.assertEqual( line, "$CTSMROOTDIR/surfdata_$YYMMDD.namelist\n" )
     line = tfile.readline()
     self.assertEqual( line, "++++ $PTCLMVERSION +++++\n" )
     tfile.close()
     os.system( "/bin/rm "+filename )


if __name__ == '__main__':
     unittest.main()
