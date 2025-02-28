"""
Remove duplicate tests from testlist_clm.xml
"""

import os
import xml.etree.ElementTree as ET
from xml.dom import minidom
import argparse


def get_front_matter(file_in):
    """
    Read the input file and preserve the stuff at the top.

    Args:
        input_file (str): Path to the input XML file.

    Returns:
        list: List of lines representing the front matter.
    """
    with open(file_in, "r", encoding="utf-8") as file:
        lines = file.readlines()
    front_matter = []
    for line in lines:
        front_matter.append(line)
        if line.strip().startswith("-->"):
            break
    return front_matter


def parse_xml(file_in):
    """
    Parse the input XML file and preserve comments.

    Args:
        input_file (str): Path to the input XML file.

    Returns:
        ElementTree: Parsed XML tree.
    """
    xml_parser = ET.XMLParser(target=ET.TreeBuilder(insert_comments=True))
    tree = ET.parse(file_in, parser=xml_parser)
    return tree


def reparse(tree):
    """
    Convert the XML tree to a minidom Document for pretty printing.

    Args:
        tree (ElementTree): Parsed XML tree.

    Returns:
        Document: Minidom Document object.
    """
    rough_string = ET.tostring(tree.getroot(), "utf-8")
    reparsed = minidom.parseString(rough_string)
    return reparsed


def xml_to_text(tree):
    """
    Convert the XML tree to a formatted string with preserved comments and fixed indentation.

    Args:
        tree (ElementTree): Parsed XML tree.

    Returns:
        str: Formatted XML string.
    """
    reparsed = reparse(tree)
    pretty_xml = reparsed.toprettyxml(indent="", newl="")
    # Replace &quot; with "
    pretty_xml = pretty_xml.replace("&quot;", '"')
    # Skip version, which we did separately above the comment
    pretty_xml = pretty_xml.replace('<?xml version="1.0" ?>', "")
    # Fix indentation difference introduced by removal of duplicates
    pretty_xml = pretty_xml.replace("      </machines>", "    </machines>")
    return pretty_xml


def print_test_info(test):
    """
    Print all the information associated with one test request.

    Args:
        test (Element): XML element representing a test.
    """
    print("Tests:")
    for machine in test.find("machines").findall("machine"):
        this_str = "  " + " ".join(
            [machine.get("name"), machine.get("compiler"), machine.get("category")]
        )
        print(this_str)


def print_tests_info(tree):
    """
    Print all the information associated with each test.

    Args:
        tree (ElementTree): Parsed XML tree.
    """
    root = tree.getroot()
    for test in root.findall("test"):
        print(f"Test name: {test.get('name')}")
        print(f"Grid: {test.get('grid')}")
        print(f"Compset: {test.get('compset')}")
        print(f"Testmods: {test.get('testmods')}")
        print_test_info(test)
        print("Options:")
        if test.find("options") is None:
            continue
        for option in test.find("options").findall("option"):
            print(f"  {option.get('name')}: {option.text}")
        print()


def remove_duplicate_machines(tree):
    """
    Identify and remove duplicate machines from the XML tree.

    Args:
        tree (ElementTree): Parsed XML tree.

    Returns:
        ElementTree: XML tree with duplicate machines removed.
    """
    root = tree.getroot()
    seen = set()
    empty_tests = []
    for test in root.findall("test"):
        test_tuple = (test.get("name"), test.get("grid"), test.get("compset"), test.get("testmods"))
        machines = test.find("machines")
        duplicates = []
        for machine in machines.findall("machine"):
            machine_tuple = test_tuple + (
                machine.get("name"),
                machine.get("compiler"),
                machine.get("category"),
            )
            if machine_tuple in seen:
                duplicates.append(machine)
            else:
                seen.add(machine_tuple)
        for duplicate in duplicates:
            machines.remove(duplicate)
        if not machines.findall("machine"):
            empty_tests.append(test)
    for test in empty_tests:
        root.remove(test)

    return tree


def main():
    """
    Main function to remove duplicate tests from the input XML file and write the result to the
    output file.

    Args:
        input_file (str): Path to the input XML file.
        output_file (str): Path to the output XML file.
    """
    args = get_args()

    # Read the input file and preserve the stuff at the top
    xml_front_matter = get_front_matter(args.input_file)

    # Parse the input XML file
    xml_tree = parse_xml(args.input_file)

    # Remove duplicate machines
    xml_tree = remove_duplicate_machines(xml_tree)

    # # Print all the information associated with each test
    # print_tests_info(tree)

    # Write
    with open(args.output_file, "w", encoding="utf-8") as f_out:
        f_out.writelines(xml_front_matter)
    with open(args.output_file, "a", encoding="utf-8") as f_out:
        pretty_xml = xml_to_text(xml_tree)
        f_out.write(pretty_xml)


def get_args():
    parser = argparse.ArgumentParser(description="Remove duplicate tests from testlist_clm.xml")
    parser.add_argument(
        "-i", "--input-file", type=str, default=None, help="Path to the input XML file"
    )
    parser.add_argument(
        "-o", "--output-file", type=str, default=None, help="Path to the output XML file"
    )
    args = parser.parse_args()

    if not args.input_file:
        args.input_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            os.pardir,
            os.pardir,
            "cime_config",
            "testdefs",
            "testlist_clm.xml",
        )
    if not os.path.exists(args.input_file):
        raise FileNotFoundError(args.input_file)
    if not args.output_file:
        args.output_file = args.input_file

    if not os.path.exists(args.input_file):
        raise FileNotFoundError(args.input_file)
    return args

if __name__ == "__main__":
    main()
