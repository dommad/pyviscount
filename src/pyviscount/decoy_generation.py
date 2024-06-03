"""Decoy generation methods"""

from abc import ABC, abstractmethod
from typing import List
from functools import reduce
import numpy as np


class DecoyGenerator(ABC):

    @abstractmethod
    def generate(self):
        pass


class DecoyLevel(ABC):

    @abstractmethod
    def modify_sequences(self, sequence_dict: dict, cut_sites: List=None, no_cut_before: List=None):
        """Generate sequences that will be subjected to 
        decoy-generating action"""
        pass


class ProteinDecoy(DecoyLevel):


    def __init__(self, sequence_dict, cut_sites_list=None, no_cut_before_list=None) -> None:
        self.sequence_dict = sequence_dict

    @staticmethod
    def modify_sequences(sequence_dict):
        return sequence_dict.values()


class PeptideDecoy(DecoyLevel):

    def __init__(self, sequence_dict, cut_sites_list=None, no_cut_before_list=None) -> None:
        self.sequence_dict = sequence_dict
        self.cut_list = cut_sites_list
        self.no_cut_list = no_cut_before_list
    
    def modify_sequences(self):
        
        return [self.process_single_sequence(seq) for seq in self.sequence_dict.values()]
        
      
    def process_single_sequence(self, seq):

        length = len(seq)
        cut_sites = [idx + 1 for idx, aa in enumerate(seq[:-1]) if aa in self.cut_list and seq[idx + 1] not in self.no_cut_list]
        cut_sites_all = [0] + cut_sites + [length]

        return [seq[cut_sites_all[idx]:cut_sites_all[idx + 1]] for idx in range(len(cut_sites_all) - 1)]
        

class ShuffledDecoy(DecoyGenerator):
    
    def __init__(self, sequence_dict: dict, decoy_level: DecoyLevel, cut_sites_list=None, no_cut_before_list=None) -> None:
        self.sequence_dict = sequence_dict
        self.decoy_level = decoy_level(sequence_dict, cut_sites_list, no_cut_before_list)


    def generate(self):
        
        shuffled_sequences = []
        sequences_to_process = self.decoy_level.modify_sequences()

        for (seq_name, seq_list) in zip(self.sequence_dict.keys(), sequences_to_process):
            for seq in seq_list:
                shuffled_seq = self.shuffle_sequence(seq)
                
                if shuffled_seq == seq:
                    shuffled_seq = self.reshuffle_sequence(shuffled_seq, seq, 3)
                    if shuffled_seq == seq:
                        continue

                entry = (seq_name, shuffled_seq)
                shuffled_sequences.append(entry)
        
        return shuffled_sequences


    def reshuffle_sequence(self, shuffled_seq, org_seq, n_rep):

        rep_idx = 0
        while (shuffled_seq == org_seq) and rep_idx < n_rep:
            shuffled_seq = self.shuffle_sequence(org_seq)
            rep_idx += 1
        
        return shuffled_seq
                    

    @staticmethod
    def shuffle_sequence(pep_seq):

        aa_list = list(pep_seq)
        shuffled_indices = np.arange(len(pep_seq)-1)
        np.random.shuffle(shuffled_indices)

        return "".join([aa_list[i] for i in shuffled_indices]) + aa_list[-1]






class ReverseDecoy(DecoyGenerator):
    
    def __init__(self, decoy_level: DecoyLevel) -> None:
        pass

    def generate(self):
        pass


    def reverse_sequence(self):
        pass


        

