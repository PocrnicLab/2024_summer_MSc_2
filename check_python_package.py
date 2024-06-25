#!/usr/bin/env python3

import os
import re

def find_python_files(directory):
    """Recursively find all Python files in the directory."""
    python_files = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.py'):
                python_files.append(os.path.join(root, file))
    return python_files

def extract_packages(file_path):
    """Extract imported packages from a Python file."""
    packages = set()
    import_pattern = re.compile(r'^\s*(import|from)\s+([a-zA-Z0-9_\.]+)')
    with open(file_path, 'r', encoding='utf-8') as file:
        for line in file:
            match = import_pattern.match(line)
            if match:
                package = match.group(2).split('.')[0]
                packages.add(package)
    return packages

def save_packages(packages, output_file):
    """Save the list of packages to a file."""
    with open(output_file, 'w', encoding='utf-8') as file:
        for package in sorted(packages):
            file.write(package + '\n')

def main(directory, output_file):
    """Main function to find and save packages."""
    python_files = find_python_files(directory)
    all_packages = set()
    for py_file in python_files:
        packages = extract_packages(py_file)
        all_packages.update(packages)
    save_packages(all_packages, output_file)

if __name__ == '__main__':
    current_directory = os.getcwd()
    output_filename = 'python_packages_list.txt'
    main(current_directory, output_filename)
