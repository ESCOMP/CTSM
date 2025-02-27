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
    # Replace this weirdness
    pretty_xml = pretty_xml.replace('<?xml version="1.0" ?>', '')
    # Fix indentation difference introduced by removal of duplicates
    pretty_xml = pretty_xml.replace("      </machines>", "    </machines>")
    f.write(pretty_xml)

# Print all the information associated with one test request
def print_test_info(test):
    print("Tests:")
    for machine in test.find('machines').findall('machine'):
        this_str = "  " + " ".join([machine.get('name'), machine.get('compiler'), machine.get('category')])
        print(this_str)

# Print all the information associated with each test
def print_tests_info(tree):
    root = tree.getroot()
    for test in root.findall('test'):
        print(f"Test name: {test.get('name')}")
        print(f"Grid: {test.get('grid')}")
        print(f"Compset: {test.get('compset')}")
        print(f"Testmods: {test.get('testmods')}")
        print_test_info(test)
        print("Options:")
        if test.find('options') is None:
            continue
        for option in test.find('options').findall('option'):
            print(f"  {option.get('name')}: {option.text}")
        print()
    return

# Identify and remove duplicate machines
def remove_duplicate_machines(tree):
    root = tree.getroot()
    seen = set()
    for test in root.findall('test'):
        test_tuple = (test.get('name'), test.get('grid'), test.get('compset'), test.get('testmods'))
        machines = test.find('machines')
        duplicates = []
        for machine in machines.findall('machine'):
            machine_tuple = test_tuple + (machine.get('name'), machine.get('compiler'), machine.get('category'))
            if machine_tuple in seen:
                duplicates.append(machine)
            else:
                seen.add(machine_tuple)
        for duplicate in duplicates:
            machines.remove(duplicate)

    return tree

# Define the input and output file paths
input_file = "cime_config/testdefs/testlist_clm.xml"
# output_file = "test.xml"
output_file = input_file

# Read the input file and preserve the stuff at the top
front_matter = get_front_matter(input_file)

# Parse the input XML file
tree = parse_xml(input_file)

# Remove duplicate machines
tree = remove_duplicate_machines(tree)

# # Print all the information associated with each test
# print_tests_info(tree)

# Write
with open(output_file, 'w', encoding='utf-8') as f:
    f.writelines(front_matter)
with open(output_file, 'a', encoding='utf-8') as f:
    write_xml(tree, f)