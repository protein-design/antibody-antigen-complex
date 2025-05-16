import logging
import os
import sys

project_dir = os.path.dirname(os.getcwd())
print('project direcotry is ', project_dir)
if project_dir not in sys.path:
    sys.path.append(project_dir)

from pprint import pprint

from src.parse_abag import ParseAbAg
from src.analyze_sasa import AnalyzeSasa
from src.analyze_dist import AnalyzeDist

def main(pdb_file, outputs_dir):
    # get chains
    p = ParseAbAg(pdb_file)
    # get chains
    p.get_chains()

    # filter: is antibody-antigen complex
    info = p.filter_antibody_antigen(outputs_dir)
    if not info:
        return None

    pprint(f'Try to detect complex from {pdb_file}: {info}')
    p.update_chains('chain_type')
    # split chains
    p.chains_to_pdb()
    # retrieve uniprot data
    p.save_annot()

    # calculate SASA
    freesasa_file = None
    try:
        freesasa_file = p.build_freesasa()
    except Exception as e:
        logging.error(f"{pdb_file} | {p.__str__} | build_freesasa | error={e}")
    
    if freesasa_file:
        # filter outputs with various superparameters
        for max_span in (10, 20, 30, 40, 50):
            for delta_threshold in (0, 5, 10):
                try:
                    p2= AnalyzeSasa(
                        freesasa_file,
                        max_span=max_span,
                        delta_threshold = delta_threshold
                    )
                    p2.retrieve_seq()
                except Exception as e:
                    logging.error(f"{pdb_file} | {p2.__str__} | error={e}")

    # calculate distance
    dist_file = None
    try:
        dist_file = p.build_dist()
    except Exception as e:
        logging.error(f"{pdb_file} | {p.__str__} | build_dist | error={e}")
    
    if dist_file:
        # filter outputs with various superparameters
        for max_span in (10, 20, 30, 40, 50):
            for dist_cutoff in (1,2,3,4,5):
                try:
                    p3 = AnalyzeDist(
                        dist_file,
                        max_span=max_span,
                        dist_cutoff=dist_cutoff
                    )
                    p3.retrieve_seq()
                except Exception as e:
                    logging.error(f"{pdb_file} | {p3.__str__} | error={e}")


##############################################
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