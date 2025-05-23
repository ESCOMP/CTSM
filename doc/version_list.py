"""
Define the versions we want to build
"""
import sys
import os
dir2add = os.path.join(os.path.dirname(__file__), "doc-builder")
if not os.path.exists(dir2add):
    raise FileNotFoundError(dir2add)
sys.path.insert(0, dir2add)
# pylint: disable=wrong-import-position
from doc_builder.docs_version import DocsVersion  # pylint: disable=import-error,no-name-in-module
from doc_builder.sys_utils import get_git_head_or_branch  # pylint: disable=import-error,no-name-in-module

# Branch name, tag, or commit SHA whose version of certain files we want to preserve
LATEST_REF = get_git_head_or_branch()

# List of version definitions
VERSION_LIST = [
    DocsVersion(
        short_name="latest",
        display_name="Latest development code",
        landing_version=True,
        ref=LATEST_REF,
    ),
    DocsVersion(
        short_name="release-clm5.0",
        display_name="CLM5.0",
        ref="release-clm5.0",
    ),
]
# End version definitions (keep this comment; Sphinx is looking for it)
