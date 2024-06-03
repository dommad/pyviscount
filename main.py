import argparse
import toml
from src.partipy import postsearch, presearch

def run(args):

    config = _read_config(args.configuration_file)

    if config.get('mode') == 'postsearch':
        postsearch.PostSearchPartition(config, target_file = args.target_file, td_file = args.target_decoy_file, decoy_file = args.decoy_file).run()
    else:
        raise ValueError('the mode is not supported!')

    
    
def _read_config(file_path):
    """Read parameters from the TOML configuration file"""
    with open(file_path, 'r') as file:
        config = toml.load(file)
    return config


if __name__ == "__main__":

    arg_parser = argparse.ArgumentParser(description='Validation by partition')
    arg_parser.add_argument('-conf', '--configuration_file', required=True, type=str,
                            help="Configuration file in TOML format")
    arg_parser.add_argument('-t',  '--target_file', required=True,
                        type=str, help='file(s) with results of target-only search (accepted format: pep.xml, tsv, txt, mzid)')
    arg_parser.add_argument('-td', '--target_decoy_file', required=True,
                        type=str, help="file(s) with results of target-decoy search (accepted format: pep.xml, tsv, txt, mzid")
    arg_parser.add_argument('-d',  '--decoy_file', required=True,
                        type=str, help='file(s) with results of decoy-only search (accepted format: pep.xml, tsv, txt, mzid)')
    arg_parser.add_argument('-o', '--output', type=str,
                            help="core name of all output files")
    args = arg_parser.parse_args()

    run(args)