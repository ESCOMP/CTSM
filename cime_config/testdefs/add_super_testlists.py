#!/usr/bin/env python3
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

inp = Path('testlist_clm.xml')
out = Path('testlist_clm.xml.modified')
text = inp.read_text()

pattern = re.compile(r'<machines>(.*?)</machines>', re.S)

def make_new_machine(attrs, indent, super_testlist="ctsm_release"):
    # attrs is string like ' name="derecho" compiler="intel" category="aux_clm"'
    name = re.search(r'name\s*=\s*"([^"]+)"', attrs)
    comp = re.search(r'compiler\s*=\s*"([^"]+)"', attrs)
    if not name or not comp:
        return None
    return f"{indent}<machine name=\"{name.group(1)}\" compiler=\"{comp.group(1)}\" category=\"{super_testlist}\"/>"

def process_block(block, super_testlist="ctsm_release"):
    # block is the inner content between <machines> and </machines>
    if super_testlist in block:
        return block
    # find first machine line
    m = re.search(r'(\n\s*)(<machine\s+([^/>]+)/>)', block)
    if not m:
        # cannot find machine, skip
        return block
    indent = m.group(1)
    # strip the leading newline so indent contains only spaces
    if indent.startswith('\n'):
        indent = indent[1:]
    attrs = m.group(3)
    new_line = make_new_machine(attrs, indent)
    if not new_line:
        return block

    # Preserve the existing XML formatting: add the new machine entry as a
    # separate line with the same indentation as the surrounding entries,
    # and preserve the original indentation that preceded the closing
    # </machines> tag so the closing tag lines up as before.
    closing_indent_m = re.search(r'(\n[ \t]*)\Z', block)
    closing_indent = closing_indent_m.group(1) if closing_indent_m else '\n'
    block = block.rstrip()
    return f"{block}\n{new_line}{closing_indent}"

new_text = pattern.sub(lambda mo: '<machines>' + process_block(mo.group(1)) + '</machines>', text)

# write only if changed
if new_text == text:
    print('No changes needed')
    out.write_text(text)
    sys.exit(0)

out.write_text(new_text)
print('WROTE', out)
