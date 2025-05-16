import json
import logging
import os
import sys
from pprint import pprint

project_dir = os.path.dirname(os.getcwd())
print('project direcotry is ', project_dir)
if project_dir not in sys.path:
    sys.path.append(project_dir)
print(sys.path)

from src.parse_abag import ParseAbAg


def main(pdb_file, outputs_dir):
    data = {}

    p = ParseAbAg(pdb_file)
    # get chains
    p.get_chains()

    try:
        # filter: is antibody-antigen complex
        info = p.filter_antibody_antigen(outputs_dir)
    except Exception as e:
        logging.error(f"{pdb_file} | filter abag complex | error={e}")

    try:
        if info:
            # update self.chain from self.info
            p.update_chains('chain_type')
    except Exception as e:
        logging.error(f"{pdb_file} | update chains | error={e}")

    try:
        if info:
            # retrieve uniprot data
            p.save_annot()
    except Exception as e:
        logging.error(f"{pdb_file} | retrieve annotation | error={e}")


if __name__ == '__main__':
    outdir = sys.argv[1]
    pdb_list_file = sys.argv[2]
    
    logging.basicConfig(
        filename=pdb_list_file + '.err',
        filemode='a', # Append mode 'w' or 'A'
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    errors, n = [], 0
    with open(pdb_list_file, 'r') as f:
        for line in f:
            pdb_file = line.rstrip()
            main(pdb_file, outdir)
            n += 1