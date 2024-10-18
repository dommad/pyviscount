"""Analysis of the search space"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List, Tuple, Dict
from collections import Counter
import numpy as np
from .utils import count_pattern_occurrences



@dataclass
class AminoAcidData:
    """Class responsible for storing amino acid data.
    aa_dict: dictionary with amino acid: mass pairs
    protein_dict: dictionary obtained from analysis of a database"""
    aa_dict: Dict
    protein_dict: Dict


@dataclass
class AminoAcidStatistics:
    """A dataclass containing useful statistics about
    amino acids from the analyzed search space"""
    aa_counts: List[Counter]
    total_aa_counts: List[Tuple[str, int]]
    aa_proportion: np.array
    all_items_mass: List[float]
    tryptic_site_counts: List[List[int]]



class AminoAcidCounter:


    def __init__(self, amino_acid_data: AminoAcidData) -> None:
        self.aa_data = amino_acid_data
        self.aa_counts = []

    def count_amino_acids_in_all_sequences(self) -> None:
        self.aa_counts = [self.count_amino_acids(x) for x in self.aa_data.protein_dict.values()]

    def count_amino_acids(self, sequence):
        return Counter(char.upper() for char in sequence if char.isalpha())
    


class AminoAcidStatisticsCalculator:

    def __init__(self, aa_counter: AminoAcidCounter) -> None:
        self.aa_counter = aa_counter


    def get_statistics(self) -> AminoAcidStatistics:

        self.aa_counter.count_amino_acids_in_all_sequences()

        aa_counts = self.aa_counter.aa_counts
        total_aa_counts = self.calculate_total_aa_counts(aa_counts)
        aa_proportion = self.calculate_aa_proportion(total_aa_counts)
        all_items_mass = self.calculate_all_items_mass()
        tryptic_site_counts = self.calculate_total_tryptic_site_counts()

        statistics = AminoAcidStatistics(
            aa_counts=aa_counts,
            total_aa_counts=total_aa_counts,
            aa_proportion=aa_proportion,
            all_items_mass=all_items_mass,
            tryptic_site_counts=tryptic_site_counts
        )

        return statistics


    def calculate_total_aa_counts(self, aa_counts):

        total_aa_counts = []
        for amino_acid in self.aa_counter.aa_data.aa_dict:
            all_counts = sum(item[amino_acid] for item in aa_counts)
            total_aa_counts.append((amino_acid, all_counts))

        return total_aa_counts
    

    def calculate_total_tryptic_site_counts(self):

        individual_counts = [self.count_tryptic_and_proline_sites(x) for x in self.aa_counter.aa_data.protein_dict.values()]
        return self.sum_tryptic_site_counts(individual_counts)
        
    @staticmethod
    def sum_tryptic_site_counts(dict_list):
        return {key: sum(cur_dict[key] for cur_dict in dict_list) for key in ('KP', 'RP', 'K', 'R')}


    def calculate_aa_proportion(self, total_aa_counts):
        amino_acids, counts = zip(*total_aa_counts)
        counts = np.array(counts)
        aa_proportion = counts / np.sum(counts)
        return aa_proportion


    def calculate_all_items_mass(self):
        return [self.calculate_item_mass(x) for x in self.aa_counter.aa_data.protein_dict.values()]


    def count_amino_acids(self, sequence):
        return Counter(char.upper() for char in sequence if char.isalpha())


    def count_tryptic_and_proline_sites(self, sequence):
        patterns = ["KP", "RP", "K", "R"]
        return count_pattern_occurrences(sequence, patterns)


    def calculate_item_mass(self, sequence):

        amino_acids = list(sequence)
        total_mass = sum(self.aa_counter.aa_data.aa_dict.get(amino_acid, 0) for amino_acid in amino_acids)
        return total_mass

 


class AAStatisticsComparison(ABC):
    """Comparison of various aspects of two sets of data related
    to two different search spaces, e.g., two protein sequence databases.

    Attributes:
        aa_stats (List[AminoAcidStatistics]): List of AminoAcidStatistics objects.
            It should contain exactly two datasets for comparison.
    """

    def __init__(self, aa_stats: List[AminoAcidStatistics]) -> None:
        self.aa_stats = aa_stats

        # we want exactly two datasets to compare
        if len(self.aa_stats) != 2:
            raise ValueError("AAStatisticsComparison requires exactly two sets of statistics.")


    @abstractmethod
    def compare(self):
        """Abstract method to be implemented by subclasses."""
        pass



class AAProportionsComparison(AAStatisticsComparison):


    def compare(self):
        """Returning the proportions of amino acids present in the datasets 
        and the difference between them"""

        aa_proportions = dict((idx, data.aa_proportion) for idx, data in enumerate(self.aa_stats))
        # difference = aa_proportions[0] - aa_proportions[1]

        return aa_proportions # , difference



class MassDistributionComparison(AAStatisticsComparison):

    
    def compare(self):
        """Returning the masses of proteins/peptides in the datasets"""
        
        return [data.all_items_mass for data in self.aa_stats]



class TrypticSitesComparison(AAStatisticsComparison):

    def compare(self):
        """Returning the number of tryptic sites in the datasets"""
        
        return [data.tryptic_site_counts for data in self.aa_stats]


class AAComparisonRunner:

    def __init__(self, comparers: List[AAStatisticsComparison]):
        self.comparers = comparers

    def run_aa_comparisons(self, aa_statistics):

        return [comparer(aa_statistics).compare() for comparer in self.comparers]



class PartitionQualityEvaluation:
    """This class contains methods for evaluating the quality of partition"""

    def __init__(self, aa_stats: List[AminoAcidStatistics], comparison_runner: AAComparisonRunner) -> None:
        self.aa_stats = aa_stats
        self.comparison_runner = comparison_runner


    def run_analysis(self):

        output = self.comparison_runner.run_aa_comparisons(self.aa_stats)
        return output
