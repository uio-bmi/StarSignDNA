import csv
from collections import defaultdict

def filter_breakpoints(input_file, output_file, paired_output_file):
    with open(input_file, newline='') as fin:
        reader = csv.DictReader(fin, delimiter='\t')
        rows = list(reader)
        header = reader.fieldnames

    # Group by ID
    groups = defaultdict(list)
    for row in rows:
        key = row['ID']
        groups[key].append(row)

    # Find IDs with both start and end
    paired_ids = set()
    for key, group in groups.items():
        breakings = set(row['breaking'] for row in group)
        if 'start' in breakings and 'end' in breakings:
            paired_ids.add(key)

    # Write output: keep only IDs with both start and end
    with open(output_file, 'w', newline='') as fout:
        writer = csv.DictWriter(fout, fieldnames=header, delimiter='\t')
        writer.writeheader()
        for row in rows:
            if row['ID'] in paired_ids:
                writer.writerow(row)

    # Write paired output: keep only IDs that do NOT have both start and end
    with open(paired_output_file, 'w', newline='') as fout2:
        writer2 = csv.DictWriter(fout2, fieldnames=header, delimiter='\t')
        writer2.writeheader()
        for row in rows:
            if row['ID'] not in paired_ids:
                writer2.writerow(row)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Filter breakpoints: output only IDs with both start and end, or only those without.")
    parser.add_argument('--in', dest='input_file', required=True, help='Input file')
    parser.add_argument('--out', dest='output_file', required=True, help='Output file (IDs with both start and end)')
    parser.add_argument('--paired', dest='paired_output_file', required=True, help='Output file (IDs without both start and end)')
    args = parser.parse_args()
    filter_breakpoints(args.input_file, args.output_file, args.paired_output_file) 