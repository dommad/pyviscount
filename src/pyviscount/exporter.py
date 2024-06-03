"""Exporters"""

from abc import ABC, abstractmethod
import logging
import pandas as pd



class Exporter(ABC):

    def __init__(self, out_name):
        self.out_name = out_name

    @abstractmethod
    def export_to_fasta(self):
        pass

    @abstractmethod
    def export_to_tsv(self):
        pass

    @staticmethod
    def write_sequences_to_fasta(file, k_sequences, k_subset, mode='target'):
        if mode == 'target':
            prefix = ""
        else:
            prefix = "decoy_"
        for idx, (protein_id, sequence) in enumerate(k_sequences):
            protein_code = protein_id.split(" ")[0]
            file.write(f">{prefix}{protein_code}_partition_{k_subset}_{idx}\n{sequence}\n")


class SpectraSTPeptideExporter(Exporter):

    def export_to_tsv(self, peptides_split, out_name):
        for split, split_peptides in enumerate(peptides_split):
            df = pd.DataFrame(split_peptides)
            df.to_csv(f"{out_name}_pep_{split}.tsv", sep='\t', header=None, index=None)



class ProteinExporter(Exporter):


    def export_to_fasta(self, k_subset, k_sequences, mode='target'):
        try:
            file_path = f"{self.out_name}_protein_{k_subset}_{mode}.fasta"
            with open(file_path, 'w') as file:
                Exporter.write_sequences_to_fasta(file, k_sequences, k_subset, mode)
        except IOError as e:
            logging.error(f"Error exporting proteins to {file_path}: {e}")


    def export_to_tsv(self):
        NotImplementedError("This function is not implemented yet.")


class PeptideExporter(Exporter):


    def export_to_fasta(self, k_subset, k_sequences, mode='target'):
        try:
            file_path = f"{self.out_name}_peptide_{k_subset}.fasta"
            with open(file_path, 'w') as file:
                Exporter.write_sequences_to_fasta(file, k_sequences, k_subset, mode)
        except IOError as e:
            logging.error(f"Error exporting peptides to {file_path}: {e}")


    def export_to_tsv(self):
        NotImplementedError("This function is not implemented yet.")

