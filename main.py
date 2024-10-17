"""Main script running PyViscount"""

import argparse
from src.pyviscount import run
from src.pyviscount.utils import open_config


def execute(arguments):

    config = open_config(arguments.configuration_file)

    mode = config.get('validation.general', 'mode', fallback='postsearch').strip()

    if mode == 'postsearch':
        run.run_postsearch_validation(arguments.configuration_file,
                                          arguments.target_file,
                                          arguments.target_decoy_file,
                                          arguments.decoy_file)
    else:
        raise ValueError('the mode is not supported!')


if __name__ == "__main__":

    arg_parser = argparse.ArgumentParser(description='Validation by partition')
    arg_parser.add_argument('-conf', '--configuration_file',
                            required=True, type=str,
                            help="Configuration file in TOML format")
    arg_parser.add_argument('-t',  '--target_file',
                            required=True, type=str,
                            help='file(s) with results of target-only search')
    arg_parser.add_argument('-td', '--target_decoy_file',
                            required=True, type=str,
                            help="file(s) with results of target-decoy search")
    arg_parser.add_argument('-d',  '--decoy_file',
                            required=True, type=str,
                            help="file(s) with results of decoy-only search")

    args = arg_parser.parse_args()

    execute(args)
