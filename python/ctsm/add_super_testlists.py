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
        "--check-only",
        help="Only check if the testlist is correct, don't write a new file "
        + "even if it needs modifications",
        action="store_true",
        dest="check_only",
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
        help="Filename to write the modified testlist to"
        + " (only writes it if there are changes required)."
        " Defaults to '<testlist-file>.modified'."
        " Set to None if the --check-only option is used",
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


def get_machine_compiler_category_from_machine_line(attrs):
    """Get the machine, compiler and category from a machine line
    Args:
        attrs (str): Attribute string extracted from an existing `<machine/>` tag
            (e.g. ' name="derecho" compiler="intel" category="aux_clm"').

    Returns:
        a tuple of name, compiler, category as strings
    """
    # attrs is string like ' name="derecho" compiler="intel" category="aux_clm"'
    name = re.search(r'name\s*=\s*"([^"]+)"', attrs)
    comp = re.search(r'compiler\s*=\s*"([^"]+)"', attrs)
    cat = re.search(r'category\s*=\s*"([^"]+)"', attrs)
    if not name or not comp or not cat:
        logger.error(
            "Could not find name, compiler and/or category attributes in machine line: %s",
            attrs,
        )
        abort("This machine block doesn't have a name, compiler or category attributes for it")
    return (name.group(1), comp.group(1), cat.group(1))


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
    (name, comp, cat) = get_machine_compiler_category_from_machine_line(attrs)
    attrs_str = f'name="{name}" compiler="{comp}" category="{super_testlist}"'
    return f"{indent}<machine {attrs_str}/>"


def process_machines_block(
    block, testlist=None, not_in_testlist=None, super_testlist="ctsm_release"
):
    """Ensure a machines block contains an entry for the given super testlist.

    The function inspects the inner contents of a `<machines>...</machines>`
    block and, if an entry for `super_testlist` is missing, creates and
    inserts a new `<machine/>` line using the machine entry's attributes
    and indentation. The closing indentation before `</machines>` is preserved.

    When testlist is entered, check that the block has the testlist, before adding
       the super_testlist machine for it to the block
    When testlist is None, ignore it and add the super_testlist machine line for it
       to every block and machine/compiler combination.

    Args:
        block (str): The inner content between `<machines>` and `</machines>`.
        testlist: The testlist to check for existance of, before adding super_testlist to the block
        not_in_testlist: testlist to verify it's NOT in before adding super_testlist to the block
        super_testlist (str): The testlist (i.e. category) to ensure is present.

    Returns:
        str: The modified block with the added `<machine/>` entries when needed.
    """
    # block is the inner content between <machines> and </machines>
    # find the matches for each machine line in the block
    MACHINE_LINE_MATCH = r"(\n\s*)(<machine\s+([^/>]+)/>)"
    machine_lines = re.findall(MACHINE_LINE_MATCH, block)
    if not machine_lines:
        # cannot find machine, die with an error as this is a problem
        logger.error(
            "No <machine/> line found in block, so cannot add a %s entry to it", super_testlist
        )
        logger.error(block)
        abort("block doesn't have a machine line for it, so aborting")

    # Ensure no duplicates
    if len(machine_lines) != len(set(machine_lines)):
        logger.error("There are duplicated tests in a machines block")
        logger.error(block)
        abort("block has duplicated tests, so aborting")
    #
    # Iterate over the machine lines to figure out, for each machine/compiler
    # combination, the full set of categories (testlists) it's already in.
    # Also remember one representative (indent, attrs) line per machine/compiler
    # to use as a template if we need to add a new entry for it.
    #
    mach_comp_cat = {}
    representative_line = {}
    for this_machine_line in machine_lines:
        indent = str(this_machine_line[0])
        attrs = str(this_machine_line[1])
        (name, comp, cat) = get_machine_compiler_category_from_machine_line(attrs)
        mach_comp_cat.setdefault(name, {}).setdefault(comp, set()).add(cat)
        representative_line.setdefault((name, comp), (indent, attrs))
    #
    # Iterate over each unique machine/compiler combination to process the block
    #
    for (name, comp), (indent, attrs) in representative_line.items():
        categories = mach_comp_cat[name][comp]

        # If the super_testlist is already in the block for this machine/compiler, skip it
        if super_testlist in categories:
            continue

        # Check if the input testlist is in the block for this machine/compiler and if not skip it
        if testlist is not None and testlist not in categories:
            continue

        if not_in_testlist is not None and not_in_testlist in categories:
            continue

        # strip the leading newline so indent contains only spaces
        if indent.startswith("\n"):
            indent = indent[1:]
        new_line = make_new_machine_line_for_supertestlist(
            attrs, indent, super_testlist=super_testlist
        )
        if not new_line:
            return block

        # Log about this block needing to be modified
        logger.warning(
            "Found a block that needs to be modified and %s added to it for %s_%s",
            super_testlist,
            name,
            comp,
        )
        logger.warning(block)

        # Preserve the existing XML formatting: add the new machine entry as a
        # separate line with the same indentation as the surrounding entries,
        # and preserve the original indentation that preceded the closing
        # </machines> tag so the closing tag lines up as before.
        closing_indent_m = re.search(r"(\n[ \t]*)\Z", block)
        closing_indent = closing_indent_m.group(1) if closing_indent_m else "\n"
        block = block.rstrip()
        block = f"{block}\n{new_line}{closing_indent}"

    return block


def add_super_testlist(text, super_testlist, testlists=(None,), not_in_testlist=None):
    """Add missing super_testlist machine entries to every `<machines>` block in text.

    Runs process_machines_block over every `<machines>` block in text once per entry in
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
            + process_machines_block(
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
        logger.info("No changes needed, so output file NOT written to")
        logger.info("Successfully validated that the testlist.xml file is correct")
        return
    #
    # If changes are needed report on how it failed
    #

    # Output modified file if the --check-only option wasn't used
    if not args.check_only:
        args.output_file.write_text(text)
        logger.warning("WROTE %s", args.output_file)
        logger.warning("Use %s to correct the testlist_clm.xml file", args.output_file)
    # Otherwise: Document the steps to do to fix the problem
    else:
        logger.warning("Modified file was NOT written, because the --check-only option was used")
        logger.warning("If this failed from a github workflow run -- do the following steps")
        logger.warning(
            "1.) Run %s To both see the fails and also create a modified file with the fixes",
            __file__,
        )
        logger.warning(
            "2.) Copy the modified file to testlist_clm.xml and verify the changes are correct"
        )
        logger.warning("3.) git commit and git push the changes")
        logger.warning("4.) Verify the workflow runs correctly now -- or repeat the process")

    # Exit with an error
    logger.warning("The testlist file had problems and needs some updates")
    abort("The testlist didn't validate")


if __name__ == "__main__":
    main()
