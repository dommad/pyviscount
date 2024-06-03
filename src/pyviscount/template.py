"""Abstract classes for validation by partition analysis"""

from typing import List
from abc import ABC, abstractmethod
import pandas as pd
import numpy as np
import scipy.stats as st
from tqdm import tqdm


class FileReading:

    def __init__(self, parser) -> None:
        self.parser = parser


    def read_search_results(self, *args) -> List[pd.DataFrame]:
        """
        Read search results for target-only, target-decoy, and decoy-only searches.

        Returns:
        tuple: A tuple containing DataFrames for target-only, target-decoy, and decoy-only search results.
        """

        dataframes = []

        for idx in tqdm(range(len(args)), desc="Reading input files"):
            try:
                df = self.parser.parse(args[idx])
                dataframes.append(df)
            except Exception as e:
                # Handle specific exceptions, log the error, and continue or raise a more informative exception
                raise ValueError(f"Error reading file {args[idx]}: {e}") from e

        return dataframes




class QualityFiltering(ABC):

    @abstractmethod
    def _remove_below_threshold(self, df: pd.DataFrame, threshold: float):
        pass

    @abstractmethod
    def _update_identification_status_labels(self, threshold: float):
        pass


class ConfidenceInterval:

    def __init__(self) -> None:
        pass

    def get_mean_cis(self, bin_results, l_cutoff, u_cutoff, confidence_level):

        truncated = [x[int(l_cutoff * len(x)): int(u_cutoff * len(x))] for x in bin_results]
        stats = np.array([self.calculate_confidence_interval(x, confidence_level) for x in truncated])

        return stats


    def calculate_confidence_interval(self, sample, conf_level):

        sample_mean = np.mean(sample)
        margin_of_error = self.calculate_margin_of_error(sample, conf_level)

        return sample_mean - margin_of_error, sample_mean, sample_mean + margin_of_error


    @staticmethod
    def calculate_margin_of_error(sample, conf_level):
        sample_std = np.std(sample, ddof=1)
        alpha = 1 - conf_level
        critical_value = st.t.ppf((1 - alpha) / 2, df=len(sample))
        std_error = sample_std / np.sqrt(len(sample))
        margin_of_error = critical_value * std_error
        return margin_of_error



class AveragedFdpFdr:

    def __init__(self, fdr_bins) -> None:
        
        self.fdr_bins = fdr_bins
        self.range_bins = range(len(fdr_bins))


    def put_in_bins(self, fdrs, fdps):

        fdr_digitized = [np.digitize(fdr_arr, bins=self.fdr_bins) for fdr_arr in fdrs]
        zipped_fdr_bin_fdp = [np.array(list(zip(fdr_digitized[idx], fdps[idx]))) for idx in range(len(fdps))]
        results = [self.separate_zipped(x) for x in zipped_fdr_bin_fdp]
        extracted = [self.extract_from_each_bin(results, idx) for idx in self.range_bins]
        aggregated_for_each_bin = [self.aggregate_for_each_bin(extracted, idx)[0] for idx in self.range_bins]

        return aggregated_for_each_bin


    @staticmethod
    def aggregate_for_each_bin(extracted, bin_idx):
        return [[item for array in extracted[bin_idx] for item in array]]

    @staticmethod
    def extract_from_each_bin(results, idx):
        return [results[item][idx] for item in range(len(results))]


    def separate_zipped(self, zipped):
        return [zipped[zipped[:, 0] == idx + 1][:, 1] for idx in self.range_bins]
    
