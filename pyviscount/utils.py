"""Utility functions"""
import time
import logging
import configparser
from tqdm import tqdm
import numpy as np
from .constants import *


# Configure the logging settings
logging.basicConfig(
    level=logging.INFO,  # Set the logging level to INFO or desired level
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)


def get_thresholds_dictionary(num_thresholds):


    threshold_dict = {
            'tev': np.linspace(TEV_LOW, TEV_HIGH, num_thresholds),
            'refactored_xcorr': np.linspace(REF_XCORR_LOW, REF_XCORR_HIGH, num_thresholds),
            'xcorr': np.linspace(XCORR_LOW, XCORR_HIGH, num_thresholds),
            'hyperscore': np.linspace(HYPER_LOW, HYPER_HIGH, num_thresholds),
            'fval': np.linspace(SPECTRAST_FVAL_LOW, SPECTRAST_FVAL_HIGH, num_thresholds),
            'KS_score': np.linspace(SPECTRAST_KS_LOW, SPECTRAST_KS_HIGH, num_thresholds),
            'dot': np.linspace(SPECTRAST_DOT_LOW, SPECTRAST_DOT_HIGH, num_thresholds),
            'percolator_score': np.linspace(PERCOLATOR_LOW, PERCOLATOR_HIGH, num_thresholds),
            'msgf_raw_score': np.linspace(MSGF_RAW_LOW, MSGF_RAW_HIGH, num_thresholds),
            'msgf_denovo_score': np.linspace(MSGF_DENOVO_LOW, MSGF_DENOVO_HIGH, num_thresholds),
            'msgf_spec_e_value': np.linspace(MSGF_SPEC_EV_LOW, MSGF_SPEC_EV_HIGH, num_thresholds),
            }
    
    return threshold_dict



def open_config(config_file_path: str):
    
    with open(config_file_path, 'r', encoding='utf-8') as config_file:
        config = configparser.ConfigParser()
        config.read_file(config_file)
    
    return config


def fetch_instance(class_name, attribute_name):
    """general fetches for class attributes by name and possibly initializing them"""
    try:
        return getattr(class_name, attribute_name)

    except AttributeError as exc:
        raise ValueError(f"Unsupported or invalid instance type: {class_name}") from exc


def timeit(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"{func.__name__} took {end_time - start_time} seconds.")
        return result
    return wrapper


def log_function_call(func):

    def wrapper(*args, **kwargs):
        #print(f"Calling {func.__name__} with args {args} and kwargs {kwargs}")
        result = func(*args, **kwargs)
        print(f"{func.__name__} returned {result}")
        return result

    return wrapper



def _is_numeric(value):
    if not isinstance(value, str):
        return False
    try:
        float(value)
        return True

    except ValueError:
        return False


def largest_factors(n):
    for i in range(n // 2, 0, -1):
        if n % i == 0:
            return n // i, i


class StrClassNameMeta(type):

    def __str__(cls):
        return cls.__name__


class ParserError(Exception):
    """A custom exception class."""
    def __init__(self, message="An error occurred."):
        self.message = message
        super().__init__(self.message)



def _is_numeric(value):
        if not isinstance(value, str):
            return False
        try:
            float(value)
            return True

        except ValueError:
            return False
        

def count_pattern_occurrences(sequence, patterns):
    counts = dict([(pattern, sequence.count(pattern)) for pattern in patterns])
    return counts


def tqdm_decorator(total_steps=None, desc="Progress"):

    def decorator(func):

        def wrapper(*args, **kwargs):
            # Log the start of the function
            logging.info("Starting the function...")

            if total_steps is not None:
                progress_bar = tqdm(total=total_steps, desc=desc)
            else:
                progress_bar = tqdm(desc=desc, leave=False)

            try:
                result = func(*args, **kwargs)
            finally:
                progress_bar.close()

            # Log the completion of the function
            logging.info("Function completed.")

            return result

        return wrapper

    return decorator



def generate_output_name(input_name):
    return input_name.split('/')[-1].split('.')[0]
