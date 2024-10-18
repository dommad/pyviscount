"""This module contains methods for generation
of partitioned search spaces: spectral libraries and sequence databases"""
from typing import List
from abc import ABC, abstractmethod
import numpy as np
from .exporter import Exporter
from .sequence_processing import clean_sequence
from . import utils, exporter, partition_mode



class Partition(ABC):
    
    @abstractmethod
    def export_split_search_space(self) -> None:
        pass

    @abstractmethod
    def split_search_space(self) -> None:
        pass


class PartitionSeqDB(Partition):
    """Partitioning the sequence database. It can be executed using different modes specified
    by the 'mode' argument: protein level or peptide level."""

    def __init__(self, config, protein_dict: dict) -> None:

        mode = config.get('partition.general', 'partition_mode', fallback='Protein').split(" ")[0].strip()
        self.partition_mode = getattr(partition_mode, mode.capitalize() + "LevelPartition")
        self.exporter_instance = getattr(exporter, mode.capitalize()+ "Exporter")
    
        self.n_subsets = int(config.get('partition.general', 'num_partitions', fallback=2).split(" ")[0].strip())
        self.mc = int(config.get('partition.general', 'num_missed_cleavages', fallback=2).split(" ")[0].strip())
        self.output_path = config.get('partition.general', 'output_path').split(" ")[0].strip()
        self.protein_dict = protein_dict
        self.split_target_sequences: List = None

        self.core_name = config.get('partition.general', 'search_space_path', fallback='test').split(" ")[0].strip().split('/')[-1].split('.')[0]


    def split_search_space(self):
        self.split_target_sequences = self.partition_mode.generate_partitions(self.protein_dict, self.n_subsets, self.mc)
        return self


    def export_split_search_space(self):

        if not self.split_target_sequences or self.n_subsets == 0:
            return  # No need to proceed if there are no results or subsets

        for k_subset in range(self.n_subsets):
            k_sequences = (x[1] for x in self.split_target_sequences if x[0] == k_subset)
            self.exporter_instance(self.output_path+self.core_name).export_to_fasta(k_subset, k_sequences)


class PartitionSpTXT(Partition):

    def __init__(self, input_file, n_subsets, library_spectra, exporter_instance: Exporter) -> None:
        self.input_file = input_file
        self.n_subsets = n_subsets
        self.library_spectra = library_spectra
        self.out_name = utils.generate_output_name(input_file)
        self.exporter_instance = exporter_instance


    def split_search_space(self, n_split=2):

        peptides = [clean_sequence(x.fullname) for x in self.library_spectra]
        indices = np.random.randint(0, n_split, len(peptides))

        return [peptides[indices == split] for split in range(n_split)]


    def export_split_search_space(self, peptides_split):
        self.exporter_instance.export_to_tsv(peptides_split, self.out_name)

