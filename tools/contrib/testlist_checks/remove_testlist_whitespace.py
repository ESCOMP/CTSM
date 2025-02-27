# %%
import xml.etree.ElementTree as ET
from xml.dom import minidom

# Read the input file and preserve the stuff at the top
def get_front_matter(input_file):
    with open(input_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    front_matter = []
    for line in lines:
        front_matter.append(line)
        if line.strip().startswith("-->"):
            break
    return front_matter

# Parse the input XML file and preserve comments
def parse_xml(input_file):
    parser = ET.XMLParser(target=ET.TreeBuilder(insert_comments=True))
    tree = ET.parse(input_file, parser=parser)
    return tree

# Write the preserved comments and parsed XML tree to the output file
def write_xml(tree, f):
    rough_string = ET.tostring(tree.getroot(), 'utf-8')
    reparsed = minidom.parseString(rough_string)
    pretty_xml = reparsed.toprettyxml(indent="", newl="")
    # Replace &quot; with "
    pretty_xml = pretty_xml.replace("&quot;", '"')
    # Skip version, which we did separately above the comment
    pretty_xml = pretty_xml.replace('<?xml version="1.0" ?>', '')
    f.write(pretty_xml)

# Define the input and output file paths
input_file = "cime_config/testdefs/testlist_clm.xml"
output_file = input_file

# Read the input file and preserve the stuff at the top
front_matter = get_front_matter(input_file)

# Parse the input XML file
tree = parse_xml(input_file)

# Write
with open(output_file, 'w', encoding='utf-8') as f:
    f.writelines(front_matter)
with open(output_file, 'a', encoding='utf-8') as f:
    write_xml(tree, f)