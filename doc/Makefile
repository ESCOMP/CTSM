# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = clmdoc
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# 'make fetch-images' should be run before building the documentation. (If building via
# the build_docs command, this is run automatically for you.) This is needed because we
# have configured this repository (via an .lfsconfig file at the top level) to NOT
# automatically fetch any of the large files when cloning / fetching.
fetch-images:
	git lfs install
	git lfs pull --exclude="" --include=""

.PHONY: help fetch-images Makefile

# For the pdf, we only build the tech note, but use the conf.py file in
# the top-level source directory (the -c option says where to find
# conf.py). Note that we also override the setting of
# numfig_secnum_depth in order to have figure numbering as desired in
# the pdf, given that the pdf just contains the tech note, so doesn't
# have the top-level numbering present in the web documentation (where
# top-level section 1 is the User's Guide and section 2 is the Tech
# Note).
#
# The use of $(0) is as in the catch-all target.
latexpdf:
	$(SPHINXBUILD) -M $@ "$(SOURCEDIR)/tech_note" "$(BUILDDIR)" -c "$(SOURCEDIR)" -D numfig_secnum_depth=1 $(SPHINXOPTS) $(O)

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
