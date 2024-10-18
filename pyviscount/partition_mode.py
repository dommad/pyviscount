"""Partition modes available for the analysis"""

from typing import List, Tuple
import logging
from abc import ABC, abstractmethod
import numpy as np
from .sequence_processing import digest_by_trypsin




class PartitionMode(ABC):
    """Blueprint for a class responsible for partition in a certain way"""

    @abstractmethod
    def generate_partitions(self, protein_dict, n_subsets) -> List[Tuple[int, Tuple]]:
        pass


class ProteinLevelPartition(PartitionMode):
    """Partition on protein level. The set of proteins from the original sequence database
    is split into n_subsets randomly."""

    @staticmethod
    def generate_partitions(protein_dict, n_subsets, *args, **kwargs) -> List[Tuple[int, Tuple]]:
        subset_indices = np.random.randint(0, n_subsets, len(protein_dict))
        return list(zip(subset_indices, protein_dict.items()))


class PeptideLevelPartition(PartitionMode):
    """Partition on peptide level. Proteins from the search space are first digested into
    peptides, and those peptides are then split into n_subsets. The origin of the peptide
    is tracked, so the info is stored as a tuple (protein, peptide) in case we want to 
    trace back which protein given peptide originated from. This also necessitates removal
    of possible duplicates that is done randomly, i.e., if multiple proteins have the same
    peptide, the peptide's protein identifier is selected randomly from the available values."""
    
    @staticmethod
    def generate_partitions(protein_dict: dict, n_subsets: int, num_missed_cleavages: int) -> List[Tuple[int, Tuple]]:
        peptide_list = []
        for protein_id, protein_seq in protein_dict.items():
            cur_peptides = digest_by_trypsin(protein_seq, num_missed_cleavages)
            peptide_list.extend([(protein_id, pep) for pep in cur_peptides])

        # remove duplicates
        peptide_set = set(peptide_list)
        logging.warning(f"{len(peptide_list) - len(peptide_set)} duplicate peptides removed! Total peptides now: {len(set(peptide_set))}")
        subset_indices = np.random.randint(0, n_subsets, len(peptide_set))

        return list(zip(subset_indices, peptide_set))


