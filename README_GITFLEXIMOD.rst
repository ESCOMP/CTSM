Obtaining the full model code and associated scripting infrastructure
=====================================================================

CTSM is released via GitHub. You will need some familiarity with git in order
to modify the code and commit these changes. However, to simply checkout and run the
code, no git knowledge is required other than what is documented in the following steps.

To obtain the CTSM code you need to do the following:

#. Clone the repository. ::

      git clone https://github.com/escomp/ctsm.git my_ctsm_sandbox

   This will create a directory ``my_ctsm_sandbox/`` in your current working directory.

#. Run **./bin/git-fleximod update**. ::

      cd my_ctsm_sandbox
      ./bin/git-fleximod update
      ./bin/git-fleximod --help  # for a user's guide

   **git-fleximod** is a package manager that will
   populate the ctsm directory with the relevant versions of each of the
   components along with the CIME infrastructure code.
   Additional documentation for git-fleximod appears here:
   https://github.com/ESMCI/git-fleximod?tab=readme-ov-file#git-fleximod

At this point you have a working version of CTSM.

To see full details of how to set up a case, compile and run, see the CIME documentation at http://esmci.github.io/cime/ .

More details on git-fleximod
----------------------------

The file **.gitmodules** in your top-level CTSM directory tells
**git-fleximod** which tag/branch of each component should be
brought in to generate your sandbox.

NOTE: If you manually modify an external without updating .gitmodules,
e.g. switch to a different tag, then rerunning git-fleximod will warn you of
local changes you need to resolve.
git-fleximod will not change a modified external back to what is specified in
.gitmodules without the --force option.
See below documentation `Customizing your CTSM sandbox`_ for more details.

**You need to rerun git-fleximod whenever .gitmodules has
changed** (unless you have already manually updated the relevant
external(s) to have the correct branch/tag checked out). Common times
when this is needed are:

* After checking out a new CTSM branch/tag

* After merging some other CTSM branch/tag into your currently
  checked-out branch

Customizing your CTSM sandbox
=============================

There are several use cases to consider when you want to customize or modify your CTSM sandbox.

Switching to a different CTSM branch or tag
-------------------------------------------

If you have already checked out a branch or tag and **HAVE NOT MADE ANY
MODIFICATIONS** it is simple to change your sandbox. Say that you
checked out ctsm1.0.0 but really wanted to have ctsm1.1.0;
you would simply do the following::

  git checkout ctsm1.1.0
  ./bin/git-fleximod update

You should **not** use this method if you have made any source code
changes, or if you have any ongoing CTSM cases that were created from
this sandbox. In these cases, it is often easiest to do a second **git
clone**.

Pointing to a different version of a component
----------------------------------------------

Each entry in **.gitmodules** has the following form (we use CIME as an
example below)::

  [submodule "cime"]
  path = cime
  url = https://github.com/ESMCI/cime
  fxtag = cime6.0.246
  fxrequired = ToplevelRequired
  fxDONOTUSEurl = https://github.com/ESMCI/cime

Each entry specifies either a tag or a hash. To point to a new tag or hash:

#. Modify the relevant entry/entries in **.gitmodules** (e.g., changing
   ``cime6.0.246`` to ``cime6.0.247`` above)

#. Checkout the new component(s)::

     ./bin/git-fleximod update <component>

Keep in mind that changing individual components from a tag may result
in an invalid model (won't compile, won't run, not scientifically
meaningful) and is unsupported.

Committing your change to .gitmodules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After making this change, it's a good idea to commit the change in your
local CTSM git repository. First create a branch in your local
repository, then commit it. (Unlike with subversion, branches are stored
locally unless you explicitly push them up to GitHub. Feel free to
create whatever local branches you'd like.) For example::

  git checkout -b my_ctsm_branch
  git add .gitmodules
  git commit -m "Update CIME to cime5.4.0-alpha.20"

