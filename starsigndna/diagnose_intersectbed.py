#!/usr/bin/env python3
"""
Diagnostic script to check intersectBed output format and identify potential issues.
Run this on the computer where the script produces empty files.
"""

import sys
import os

def diagnose_file(filename):
    """Diagnose the format and content of an intersectBed output file."""
    print(f"=== Diagnosing file: {filename} ===")
    
    if not os.path.exists(filename):
        print(f"ERROR: File {filename} does not exist!")
        return False
    
    file_size = os.path.getsize(filename)
    print(f"File size: {file_size} bytes")
    
    if file_size == 0:
        print("ERROR: File is empty!")
        return False
    
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        print(f"Total lines: {len(lines)}")
        
        if len(lines) == 0:
            print("ERROR: No lines found in file!")
            return False
        
        # Check first few lines
        print("\n=== First 5 lines ===")
        for i, line in enumerate(lines[:5]):
            print(f"Line {i+1}: {repr(line.strip())}")
        
        # Check line structure
        print("\n=== Line structure analysis ===")
        for i, line in enumerate(lines[:10]):  # Check first 10 lines
            if line.strip() == "":
                print(f"Line {i+1}: Empty line")
                continue
            
            fields = line.strip().split('\t')
            print(f"Line {i+1}: {len(fields)} fields")
            
            if len(fields) < 10:
                print(f"  WARNING: Line {i+1} has fewer than 10 fields!")
                print(f"  Fields: {fields}")
            else:
                # Check if it looks like intersectBed output
                chr_col = fields[0]
                start_col = fields[1]
                end_col = fields[2]
                name_col = fields[3]
                attributes = fields[-1] if len(fields) > 0 else ""
                
                print(f"  Chr: {chr_col}, Start: {start_col}, End: {end_col}")
                print(f"  Name: {name_col}")
                print(f"  Attributes (last field): {attributes[:100]}...")
                
                # Check for gene_name in attributes
                if 'gene_name=' in attributes:
                    print(f"  ✓ Contains gene_name attribute")
                else:
                    print(f"  ✗ No gene_name attribute found")
                
                # Check name format
                if '|' in name_col:
                    parts = name_col.split('|')
                    print(f"  ✓ Name contains '|' separators ({len(parts)} parts)")
                else:
                    print(f"  ✗ Name does not contain '|' separators")
        
        # Check for common issues
        print("\n=== Common issues check ===")
        
        # Check for different line endings
        cr_count = sum(1 for line in lines if '\r' in line)
        lf_count = sum(1 for line in lines if '\n' in line)
        crlf_count = sum(1 for line in lines if '\r\n' in line)
        
        print(f"Lines with CR: {cr_count}")
        print(f"Lines with LF: {lf_count}")
        print(f"Lines with CRLF: {crlf_count}")
        
        # Check encoding issues
        try:
            with open(filename, 'r', encoding='utf-8') as f:
                f.read()
            print("✓ File can be read with UTF-8 encoding")
        except UnicodeDecodeError as e:
            print(f"✗ UTF-8 encoding error: {e}")
            try:
                with open(filename, 'r', encoding='latin-1') as f:
                    f.read()
                print("✓ File can be read with latin-1 encoding")
            except Exception as e2:
                print(f"✗ latin-1 encoding also failed: {e2}")
        
        return True
        
    except Exception as e:
        print(f"ERROR reading file: {e}")
        return False

def check_python_environment():
    """Check Python environment and dependencies."""
    print("=== Python Environment ===")
    print(f"Python version: {sys.version}")
    print(f"Python executable: {sys.executable}")
    
    # Check required modules
    required_modules = ['csv', 'collections']
    for module in required_modules:
        try:
            __import__(module)
            print(f"✓ {module} module available")
        except ImportError as e:
            print(f"✗ {module} module not available: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python diagnose_intersectbed.py <intersectbed_output_file>")
        sys.exit(1)
    
    filename = sys.argv[1]
    
    check_python_environment()
    print()
    
    success = diagnose_file(filename)
    
    if success:
        print("\n=== Recommendations ===")
        print("1. If the file looks correct, try running the main script with verbose output")
        print("2. Check if the file has different line endings (Windows vs Unix)")
        print("3. Verify that the intersectBed output format matches expectations")
        print("4. Ensure the file has at least 10 tab-separated fields per line")
        print("5. Check that the last field contains gene_name attributes")
    else:
        print("\n=== File has issues ===")
        print("Please fix the identified issues before running the main script.") 