#!/usr/bin/env python3
"""add_super_testlists.py

Utilities to ensure tests in `testlist_clm.xml` are included in
super-testlists (for example, `ctsm_release`).

Reads `testlist_clm.xml`, inserts missing `<machine category="..."/>`
entries into `<machines>` blocks, and writes a modified file
`testlist_clm.xml.modified`.
"""

import re
import sys
from pathlib import Path

#
# This script reads in the testlist and determines if new lines for testlists that are
# supersets of other testlists need to be added.
#
# For example all tests need to be part of the ctsm_release test list.
#
# NOTE: This script was created by Co-Pilot (GPT-5 mini) in VS-Code with Erik Kluzek
#

inp = Path("testlist_clm.xml")
out = Path("testlist_clm.xml.modified")
text = inp.read_text()

machines_pattern = re.compile(r"<machines>(.*?)</machines>", re.S)


def make_new_machine_line_for_supertestlist(attrs, indent, testlist=None, super_testlist="ctsm_release"):
    """Create a `<machine/>` XML line using attributes for super_testlist from an existing machine.

    Args:
        attrs (str): Attribute string extracted from an existing `<machine/>` tag
            (e.g. ' name="derecho" compiler="intel" category="aux_clm"').
        indent (str): The whitespace indentation to prefix the new line with.
        super_testlist (str): The testlist (i.e. category) to use for the new machine line.

    Returns:
        str or None: The formatted `<machine .../>` line (without a trailing
        newline) or `None` if required attributes are missing.
    """
    # attrs is string like ' name="derecho" compiler="intel" category="aux_clm"'
    name = re.search(r'name\s*=\s*"([^\"]+)"', attrs)
    comp = re.search(r'compiler\s*=\s*"([^"]+)"', attrs)
    if not name or not comp:
        return None
    attrs_str = f'name="{name.group(1)}" compiler="{comp.group(1)}" category="{super_testlist}"'
    return f"{indent}<machine {attrs_str}/>"


def process_machine_block(block, testlist=None, not_in_testlist=None, super_testlist="ctsm_release"):
    """Ensure a machine block contains an entry for the given super testlist.

    The function inspects the inner contents of a `<machines>...</machines>`
    block and, if an entry for `super_testlist` is missing, creates and
    inserts a new `<machine/>` line using the first machine entry's attributes
    and indentation. The closing indentation before `</machines>` is preserved.

    When testlist is entered, check that the block has the testlist, before adding
       the super_testlist machine for it to the block
    When testlist is None, ignore it and add the super_testlist machine line for it
       to every block.

    Args:
        block (str): The inner content between `<machines>` and `</machines>`.
        testlist: The testlist to check for existance of, before adding super_testlist to the block
        not_in_testlist: Another testlist to verify it's not in before adding super_testlist to the block
        super_testlist (str): The testlist (i.e. category) to ensure is present.

    Returns:
        str: The modified block with the added `<machine/>` entry when needed.
    """
    # block is the inner content between <machines> and </machines>

    # If the4 super_testlist is already in the block return the original block
    if super_testlist in block:
        return block

    # Check if the input testlist is in the block and if not return
    if testlist is None:
        pass
    elif testlist not in block:
        return block

    if not_in_testlist is None:
        pass
    elif not_in_testlist in block:
        return block

    # find first machine line
    m = re.search(r"(\n\s*)(<machine\s+([^/>]+)/>)", block)
    if not m:
        # cannot find machine, skip
        return block
    indent = m.group(1)
    # strip the leading newline so indent contains only spaces
    if indent.startswith("\n"):
        indent = indent[1:]
    attrs = m.group(3)
    new_line = make_new_machine_line_for_supertestlist(attrs, indent, super_testlist=super_testlist)
    if not new_line:
        return block

    # Preserve the existing XML formatting: add the new machine entry as a
    # separate line with the same indentation as the surrounding entries,
    # and preserve the original indentation that preceded the closing
    # </machines> tag so the closing tag lines up as before.
    closing_indent_m = re.search(r"(\n[ \t]*)\Z", block)
    closing_indent = closing_indent_m.group(1) if closing_indent_m else "\n"
    block = block.rstrip()
    return f"{block}\n{new_line}{closing_indent}"


# Substitute the addition of the new machine lines for the super_testlist in each machines block
# machines_pattern catches the machines blocks (between <machines> and </machines>)
# The replacement is the argument sent to .sub before the text argument
# lambda runs process_machine_block on each machines_block labeled mo, to create the replacement text

orig_text = text
# First run for the ctsm_release super testlist which should be in every block
new_text = machines_pattern.sub(
    lambda mo: "<machines>" + process_machine_block(mo.group(1)) + "</machines>", text
)

# Now run for the ctsm_sci super testlist for each of the testlists that should be subsets of it
# Don't add to ctsm_sci if already in the aux_clm testlist either
ctsm_sci_testlists = ["hillslope", "fire", "ssp", "crop_calendars", "interim_restart", "subset_data"] 
for testlist in ctsm_sci_testlists:
   new_text = machines_pattern.sub(
       lambda mo: "<machines>" + process_machine_block(mo.group(1), testlist=testlist, not_in_testlist="aux_clm", super_testlist="ctsm_sci") + "</machines>", text
   )
   text = new_text

# Now run for the aux_clm super testlist for each of the testlists that should be subsets of it
aux_clm_testlists = ["clm_pymods", "prealpha", "prebeta", "aux_cime_baselines", "clm_short", "matrixcn", "aux_clm_mpi_serial"]
for testlist in aux_clm_testlists:
   new_text = machines_pattern.sub(
       lambda mo: "<machines>" + process_machine_block(mo.group(1), testlist=testlist, super_testlist="aux_clm") + "</machines>", text
   )
   text = new_text
# write only if changed
if new_text == orig_text:
    print("No changes needed")
    out.write_text(text)
    sys.exit(0)

out.write_text(new_text)
print("WROTE", out)
