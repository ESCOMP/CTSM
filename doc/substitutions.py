"""
Substitutions for Sphinx
"""

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.

# pylint: disable=invalid-name

#################################
### Standard Sphinx variables ###
#################################

# General information about the project.
project = "ctsm"
copyright = "2020, UCAR"  # pylint: disable=redefined-builtin
author = ""

# The short X.Y version.
version = "CLM5.0"

# The full version, including alpha/beta/rc tags.
release = "release-clm5.0"

#####################################################
### Custom variables needed for doc-builder setup ###
#####################################################

# Version label used at the top of some pages.
version_label = "CLM5.0 (CESM2.1)"

#######################################################
### Custom variables optional for doc-builder setup ###
#######################################################

tex_category = "Miscellaneous"

# Used by HTML help builder
htmlhelp = {
    "basename": "clmdocdoc", # Output file base name
}

# Used for LaTeX output
latex = {
    "target_name": "clmdoc.tex",
    "title": "CLM Documentation",
    "documentclass": "manual", # howto, manual, or own class
    "category": tex_category,
}

# Used for man_pages and texinfo_documents
mantex = {
    "name": "clmdoc",
    "title": "clmdoc Documentation",
}

# Used for texinfo_documents
tex = {
    "dirmenu_entry": "clmdoc",
    "description": "One line description of project.",
    "category": tex_category,
}