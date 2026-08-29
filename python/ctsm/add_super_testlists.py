"""add_super_testlists.py

Utilities to ensure tests in a testlist XML file (e.g. `testlist_clm.xml`) are
included in super-testlists (for example, `ctsm_release`).

Reads the testlist file, inserts missing `<machine category="..."/>` entries
into `<machines>` blocks, and writes a modified file (by default
`<testlist-file>.modified`).

This script reads in the testlist and determines if new lines for testlists that are
supersets of other testlists need to be added.

For example all tests need to be part of the ctsm_release test list.

NOTE: This script was created by Co-Pilot (GPT-5 mini) in VS-Code with Erik Kluzek

Usage:
    add_super_testlists.py [--testlist-file testlist_clm.xml] [--output FILE] [--verbose]
"""

import argparse
import logging
import re
from pathlib import Path

from ctsm.ctsm_logging import (
    setup_logging_pre_config,
    add_logging_args,
    process_logging_args,
)
from ctsm.utils import abort

logger = logging.getLogger(__name__)

MACHINES_PATTERN = re.compile(r"<machines>(.*?)</machines>", re.S)

# testlists that should be treated as a subset of the ctsm_sci super-testlist,
# unless they're already part of aux_clm
CTSM_SCI_TESTLISTS = ["hillslope", "fire", "ssp", "crop_calendars", "interim_restart"]

# testlists that should be treated as a subset of the aux_clm super-testlist
AUX_CLM_TESTLISTS = [
    "clm_pymods",
    "prealpha",
    "prebeta",
    "aux_cime_baselines",
    "clm_short",
    "matrixcn",
    "aux_clm_mpi_serial",
    "subset_data",
]

# testlists that should be treated as a subset of the crop_calendars super-testlist
CROP_CAL_TESTLISTS = ["rxcropmaturity"]

# testlists that should be treated as a subset of the fates super-testlist
FATES_TESTLISTS = ["fates-landuse"]


def get_parser():
    """
    Get the parser object for add_super_testlists.py.

    Returns:
        parser (ArgumentParser):
            ArgumentParser which includes all the parser information.
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "--testlist-file",
        help="Testlist XML file to read in and check for missing super-testlist entries.",
        action="store",
        dest="testlist_file",
        type=str,
        default="testlist_clm.xml",
    )
    parser.add_argument(
        "--output",
        help="File to write the (possibly modified) testlist to."
        " Defaults to '<testlist-file>.modified'.",
        action="store",
        dest="output_file",
        type=str,
        default=None,
    )

    add_logging_args(parser)
    return parser


def process_and_check_args(args):
    """Process and check the arguments"""
    args.testlist_file = Path(args.testlist_file)
    if not args.testlist_file.is_file():
        abort(f"testlist file not found: {args.testlist_file}")

    if args.output_file is None:
        args.output_file = Path(str(args.testlist_file) + ".modified")
    else:
        args.output_file = Path(args.output_file)

    return args


def make_new_machine_line_for_supertestlist(attrs, indent, super_testlist="ctsm_release"):
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
    name = re.search(r'name\s*=\s*"([^"]+)"', attrs)
    comp = re.search(r'compiler\s*=\s*"([^"]+)"', attrs)
    if not name or not comp:
        logger.warning(
            "Could not find name and/or compiler attributes in machine line, "
            "so cannot add a %s entry for it: %s",
            super_testlist,
            attrs,
        )
        return None
    attrs_str = f'name="{name.group(1)}" compiler="{comp.group(1)}" category="{super_testlist}"'
    return f"{indent}<machine {attrs_str}/>"


def process_machine_block(
    block, testlist=None, not_in_testlist=None, super_testlist="ctsm_release"
):
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
        not_in_testlist: testlist to verify it's NOT in before adding super_testlist to the block
        super_testlist (str): The testlist (i.e. category) to ensure is present.

    Returns:
        str: The modified block with the added `<machine/>` entry when needed.
    """
    # block is the inner content between <machines> and </machines>

    # If the super_testlist is already in the block return the original block
    if super_testlist in block:
        return block

    # Check if the input testlist is in the block and if not return
    if testlist is not None and testlist not in block:
        return block

    if not_in_testlist is not None and not_in_testlist in block:
        return block

    # find first machine line
    match = re.search(r"(\n\s*)(<machine\s+([^/>]+)/>)", block)
    if not match:
        # cannot find machine, skip
        logger.debug(
            "No <machine/> line found in block, so cannot add a %s entry to it", super_testlist
        )
        return block
    indent = match.group(1)
    # strip the leading newline so indent contains only spaces
    if indent.startswith("\n"):
        indent = indent[1:]
    attrs = match.group(3)
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


def add_super_testlist(text, super_testlist, testlists=(None,), not_in_testlist=None):
    """Add missing super_testlist machine entries to every `<machines>` block in text.

    Runs process_machine_block over every `<machines>` block in text once per entry in
    testlists, accumulating the changes from each pass.

    Args:
        text (str): Full contents of the testlist XML file.
        super_testlist (str): The testlist (i.e. category) to ensure is present.
        testlists: Iterable of testlists to check for existence of (one at a time) before
            adding super_testlist to a given block. A single `(None,)` (the default) adds
            super_testlist to every block.
        not_in_testlist: Another testlist to verify a block is not in before adding
            super_testlist to it.

    Returns:
        str: The full, updated contents of the testlist XML file.
    """
    for testlist in testlists:
        text = MACHINES_PATTERN.sub(
            lambda mo, testlist=testlist: "<machines>"
            + process_machine_block(
                mo.group(1),
                testlist=testlist,
                not_in_testlist=not_in_testlist,
                super_testlist=super_testlist,
            )
            + "</machines>",
            text,
        )
    return text


def main():
    """Main function: read a testlist file and add missing super-testlist machine entries"""
    setup_logging_pre_config()
    parser = get_parser()
    args = parser.parse_args()

    process_logging_args(args)

    args = process_and_check_args(args)

    orig_text = args.testlist_file.read_text()
    text = orig_text

    # ctsm_release should be in every block
    text = add_super_testlist(text, "ctsm_release")

    # ctsm_sci: don't add it to a block already in the aux_clm testlist
    text = add_super_testlist(
        text, "ctsm_sci", testlists=CTSM_SCI_TESTLISTS, not_in_testlist="aux_clm"
    )

    text = add_super_testlist(text, "aux_clm", testlists=AUX_CLM_TESTLISTS)
    text = add_super_testlist(text, "crop_calendars", testlists=CROP_CAL_TESTLISTS)
    text = add_super_testlist(text, "fates", testlists=FATES_TESTLISTS)

    if text == orig_text:
        logger.info("No changes needed")
        args.output_file.write_text(orig_text)
        return

    args.output_file.write_text(text)
    logger.info("WROTE %s", args.output_file)


if __name__ == "__main__":
    main()
