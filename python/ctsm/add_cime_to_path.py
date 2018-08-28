"""Adds cime lib to path

Any file that can potentially be run as a top-level script (with an if __name__ ==
'__main__' block) that needs cime should import this. This includes unit test modules.

This should be the very first ctsm module imported. That way, cime will be added to your
path before other imports. That is, your script should have:

# Standard library imports
# Then something like this:
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'python')
sys.path.insert(1, _CTSM_PYTHON)
# Then:
from ctsm import add_cime_to_path

NOTE: For top-level scripts that want to pass some argument(s) to add_cime_lib_to_path,
they should include copies of what's in this module (adding the appropriate argument(s)),
rather than actually importing this module.
"""

from ctsm.path_utils import add_cime_lib_to_path
add_cime_lib_to_path()
