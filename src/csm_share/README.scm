Filename:	 csm_share/README
Original Author: Erik Kluzek
Date:            Dec/12/1997
Description:     Description of the csm_share  module Configuration Management.
		 (code shared between the CSM models)
		 Description of the Software Configuration Management (SCM)
		 issue's unique to this module, and the procedures to deal
		 with them.
Version-Control: 

CVS: $Id$
CVS: $Source$
CVS: $Name$

############################################################################
			Contents:

I. 	SCM (Software Configuration Management) issue's:

	1.) csm_share codes MUST be identical between models: 

	2.) Different models MUST access the same version of "csm_share".

	3.) Versions of each model MUST document exactly what version of 
	"csm_share" was used.

	4.) The message passing version number file MUST be backwards compatable!

	5.) Changes to public "csm_share" should be infrequent and restricted.

II. 	SCM Procedures:

1.) Revision repository under offical CSM CVS repository..

2.) Checked out versions, under official CSM frozen code repository..

3.) Version naming conventions. (ie. share1_0, and model tag names)

4.) Restriced set of people allowed to make changes to csm_share

5.) Procedures for a change_request to public version of the csm_share codes.

6.) Private changes, and branch development changes.

III. 	Software Development coding conventions:

1.) File and subroutine naming conventions:

2.) Required documentation in the code:



############################################################################

I. 	SCM (Software Configuration Management issue's):

	There are several issue's of relevance here.

	1.) csm_share codes MUST be identical between models: 

	The reason the codes were put here in the first place is to ensure
	that all models are consistent and reproduce EXACTLY the same
	results between models.  Having the same equations in different
	codes is likely to create slightly different answers that sometimes
	produce problems, it also makes management and changes to the codes
	more difficult.

	2.) Different models MUST access the same version of "csm_share".

	Basically the same reason as above.  Using different versions of
	"csm_share" is the same as using different codes in the first place.

	3.) Versions of each model MUST document exactly what version of 
	"csm_share" was used.

	To reproduce old CSM runs with specific versions of models, the
	version of "csm_share" used must be documented as well as the 
	model version.

	4.) The message passing version number file (msg/cpl_compat.h) MUST be 
	backwards compatable!

	Details of this are described in the comments of the file itself.
	In short models reference the old names contained therein and 
	as such the existing names CAN NOT be removed or changed.

	5.) Changes to "csm_share" should be infrequent and restricted.

	Since changing the codes herein has impact on all the CSM models
	there can not be frequent changes to these codes.  Changes that
	do occur must be coordinated between all the models.  And in many
	cases changes in the codes here will require changes in all model
	components (dummy and active).

II. 	SCM Procedures:

1.) Revision repository under offical CSM CVS repository..

A.) Repository: /fs/cgd/csm/models/CVS.REPOS under CVS.

The csm_share codes will be archived in CVS under the official CSM CVS 
repository at /fs/cgd/csm/models/CVS.REPOS.  The code will reside under
shared/csm_share.  Codes managed under the CSM CVS repository can also
access the csm_share module as a submodule.  This way codes checked out
from the CVS repository will include the associated csm_share directly.
Without having to reference outside libraries or outside source codes.

B.) ChangeLog and ChangeSum files:

The files ChangeLog and ChangeSum will be kept in the module as a method of
keeping track of changes to the module.  ChangeLog is a long description of 
changes while ChangeSum is a one-line description.

When csm_share is tagged these files will be automatically updated.  The user
will be prompted for certain questions and the changes to the files will occur
automatically.  Also e-mail about the tag will automatically be sent to csm-progs.

2.) Checked out versions, under official CSM frozen code repository..

As well as accessing the code via. CVS developers can use frozen code
from /fs/cgd/csm/models/csm_share.

3.) Version naming conventions. (ie. share1_0, and model tag names)

Valid tags for the csm_share module are:

share#_#   		Where #_# refers to a major/minor revision number 
			(ie. share1_0).
share#_#_brnch_*  	For development branches.
share#_#_brnchT_*  	For frozen releases along a given development branch.

and global CSM tags can be placed on csm_share to specify releases of CSM.
Global CSM tags are applied to all CSM model components for a given CSM 
release.

csm#_#		#_# refers to the major/minor revision number of CSM.
		(ie. csm2_0).

Also, each model SHOULD place it's public release tags on the csm_share
module as well.  All tagnames should start with the component model and
the end with the version number.

For example,

ccm3_3_3
cpl4_alpha1
csim2_2_1
ncom1_3

Are all valid tags to place on the csm_share module.

By placing model tagnames as well as csm and share tag names, the codes that
make up each component model are exactly specified.  Without this there is
some ambiguity as to which version of csm_share goes with a given version of
a CSM model component.

When a version is tagged, the ChangeLog/ChangeSum files will automatically be
updated, and e-mail will be automatically sent out to csm-progs to keep people
informed of these changes.

4.) Restriced set of people allowed to make changes to csm_share

Only one software developer per model component will be added to the access 
list to change csm_share.  The people in this list should be those heavy
into CSM development and aware of global issues across CSM codes.

The set of people that are allowed to make changes are listed in the file:

/fs/cgd/csm/models/CVS.SCRIPTS/shared/csm_share/lists/checkin.list

5.) Procedures for a change_request to public version of the csm_share codes.

	A.) Agreement to the change must come from those responsible for 
	each component model as well as the CSM project leaders.

	B.) Changes should be immediately tagged with the next csh_share 
	tag name.  This will automatically update the ChangeLog/ChangeSum
	files and send mail to csm-progs about the change.  This keeps the
	CSM community aware of the changes.

	C.) A new version should immediately be checked out to 

	/fs/cgd/csm/models/csm_share

6.) Private changes, and branch development changes.

	One individual can of course make whatever changes they want to the code.
It's recommended that the CVS: tags be left in place so that it's documented
where the changes came from.  A note in the same location may be helpful to
note that files were changed.  Also CVS can be used to display the differences.
If CVS is used, CVS can also be used to apply public changes to the user modified
version.

	If several people need to apply the same sort of changes, a development
branch in CVS should be done.  This allows a method to archive the branch 
development history, as well as keep track of the differences, and also easily
merge in development changes along the public version.


II. 	Software Development coding conventions:

For the code that is placed in the csm_share module, a well defined set of

coding conventions (such as CCM3, NCOM, etc.) does not apply.  Also in most
cases the coding conventions for different CSM components is quite a bit
different.  Therefore, there are only a few restrictions we wish to make on 
coding.  The goal is to make codes readable, and well documented, and use
conventions that apply to the bulk of the CSM models.  And furthermore, as
these are shared codes we wish to simple naming conventions to make it
obvious what codes are, and what they belong with.  So the only hard
restrictions we want to place are documentation requirements, and naming
requirements.  We want the shared codes to be clean, well-written, organized,
modular, and hopefully optimized, but how those goals are achieved we leave
up to the individual developer.

1.) File and subroutine naming conventions:

* Source and associated header files should match.  

Ex. cpl_compat.F and cpl_compat.h go together.

* Subroutines/function names and filenames go together.

The names for the file should go with the names of the subroutines it contains.
But there maybe more than one subroutine per file.

Ex. cpl_compat.F contains subroutines: cpl_compat and cpl_compat_check_spval.
	
2.) Required documentation in the code:

	A.) Input and output arguments to subroutines/functions should
	be exactly and explicitly documented!  Including units!  And
	whether it's an input or output argument or input/output.

	B.) Files should list the subroutines they contained in a header for
	the file.  

	And each file should have a general header for the file as
	well as headers for each subroutine contained therein.

c Subroutines included:
c
c     cpl_compat -------------- Check that the version numbers of the messages
c                                  sent agree with what's expected.
c     cpl_compat_check_spval -- Check that the real message sent is not
c                                  set to the special value.


	C.) Name and date for the original author.

c Original Author: Erik Kluzek
c Date:            Dec/97

	D.) Version information (via RCS keywords)

	All files should have the following comment fields

c File Version information:
c
c CVS: $Id$
c CVS: $Source$
c CVS: $Name$
c

	which details the RCS revision number for the file, author, and date 
	of last change, the full name and path of the RCS file, and the 
	symbolic tag that the code was checked out with (if any).  It's 
	important to have these in here so it's exactly specified what 
	version of the codes are being used.  This is also important if 
	someone makes changes as then the differences to the standard 
	version can be exactly specified.

	Optionally, subroutines that are called on initialization may contain a 
	print statment that specifies the RCS-symbolic name, or the revision 
	name of the file.  For example...

      print *, 'This is revision: $Revision$ Tag: $Name$ of the message'
     $, ' compatability interface:"

	The above exactly detail the versions of the codes actually being used
	and regardless of the mechanism used to extract them.  And additionally
	they make it easy to check that different model component use the same
	codes.
