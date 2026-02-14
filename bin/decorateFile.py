import sys
import csv

def read_pairs(filename):
    # Reads two-column file into a dict
    with open(filename, newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        return {row[0]: row[1] for row in reader if len(row) >= 2}

def process_file(source, gene_id_map=None, gene_name_map=None):
    with open(source, newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        rows = list(reader)
    
    header = rows[0]
    data_rows = rows[1:]

    # Transform header
    if gene_id_map:
        header = ['Gene ID', 'Gene Name'] + header
    else:
        header = [header[0], 'Gene Name'] + header[1:]
    
    print('\t'.join(header))

    for row in data_rows:
        gene_id = row[0]
        gene_name = gene_name_map.get(gene_id, '') if gene_name_map else ''
        if gene_id_map:
            mapped_id = gene_id_map.get(gene_id, '')
            new_row = [mapped_id, gene_name] + row
        else:
            new_row = [gene_id, gene_name] + row[1:]
        print('\t'.join(new_row))

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Decorate file with gene IDs and names.')
    parser.add_argument('--source', required=True, help='Source file to decorate')
    parser.add_argument('--geneIdFile', help='Gene ID to array design file')
    parser.add_argument('--geneNameFile', required=True, help='Gene ID to gene name file')
    args = parser.parse_args()

    gene_id_map = read_pairs(args.geneIdFile) if args.geneIdFile else None
    gene_name_map = read_pairs(args.geneNameFile)

    process_file(args.source, gene_id_map, gene_name_map)
