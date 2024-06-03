"""Process the post-search partition data"""
from abc import ABC, abstractmethod
import numpy as np
import pandas as pd
from tqdm import tqdm
from .utils import get_thresholds_dictionary
from .template import QualityFiltering
from .fdr import FdpFdrCalculator, FdrLevel
from . import fdr



class PostSearchPartition:
    """Creating peptide partition mapping based on the 
    user-specified number of partitions"""

    def __init__(self, config) -> None:
        self.num_partitions = int(config.get('partition.general', 'num_partitions', fallback=2).split(" ")[0].strip())


    def create_peptide_subset_id_mapping(self, target_df: pd.DataFrame):

        if 'pep_mod_index' not in target_df.columns:
            raise ValueError("The pep_mod_index column is not present in the dataframe provided.")
        
        unique_pep_mod_idxs = target_df['pep_mod_index'].unique()
        partition_ids = np.random.randint(0, self.num_partitions, len(unique_pep_mod_idxs))

        return partition_ids


    @staticmethod
    def add_peptide_modification_index_to_target_df(target_df: pd.DataFrame):

        target_df['pep_mod'] = list(zip(target_df['sequence'], target_df['modifications']))
        target_df['pep_mod_index'], pep_mods = target_df['pep_mod'].factorize()
        pep_mod_dict = dict(zip(pep_mods, target_df['pep_mod_index'].values))
        return target_df, pep_mod_dict




class PartitionSearchSimulation(ABC):

    def __init__(self, config, peptide_id_mapping):

        self.thr_score = config.get("validation.general", "threshold_score").split(" ")[0].strip()
        self.subset_id =  int(config.get("validation.general", "subset_id").split(" ")[0].strip())
        self.peptide_id_mapping = peptide_id_mapping


    @abstractmethod
    def update_dataframe(self, threshold):
        pass


    def extract_scans_above_threshold(self, target_df, threshold):
        """
        Filter spectra based on a score threshold in a target-only search.

        Returns:
        list: Indices of spectra passing the score threshold in the target-only search.
        """
        top_hit_above_th = (target_df['hit_rank'] == 1) & (target_df[self.thr_score] >= threshold)
        return target_df[top_hit_above_th].scan


class TdcPartitionSearchSimulation(PartitionSearchSimulation):

    def __init__(self, config, peptide_id_mapping, pep_mod_dict, target_df, td_df) -> None:
        super().__init__(config, peptide_id_mapping)
        self.target_df = target_df
        self.td_df = td_df
        self.pep_mod_dict = pep_mod_dict
        self.td_filtered = None
        self.td_part_search = None


    def update_dataframe(self, threshold):
        """
        Retains PSMs from target-decoy search that are mapped to spectra present in the updated target dataframe;
        adds proxy ground truth labels to the dataframe.
        """
        quality_scans = self.extract_scans_above_threshold(self.target_df, threshold)
        self.td_filtered = self.filter_by_high_quality_scans(self.td_df, quality_scans)
        self.generate_pep_mod_idx_for_td()
        self.td_part_search = self.extract_partition_search_results()
        self.add_ground_truth_labels()

        return self.td_part_search


    def add_ground_truth_labels(self):
        gt_status = (~self.td_part_search['is_decoy']) & (self.td_part_search['hit_rank'] == 1)
        self.td_part_search['gt_status'] = gt_status


    @staticmethod
    def filter_by_high_quality_scans(df, quality_scans):
        return df[df.scan.isin(quality_scans)].copy()


    def extract_partition_search_results(self):
        mask = (self.td_filtered['is_decoy']) | (self.peptide_id_mapping[self.td_filtered['pep_mod_index']] == self.subset_id)
        return self.td_filtered.loc[mask].groupby('scan').head(1)


    def generate_pep_mod_idx_for_td(self):

        pep_mods = zip(self.td_filtered['sequence'], self.td_filtered['modifications'])
        pep_mods_idx = list(map(lambda x: self.pep_mod_dict.get(x, -1), pep_mods))
        self.td_filtered['pep_mod_index'] = pep_mods_idx


class DecoyfreePartitionSearchSimulation(PartitionSearchSimulation):
    
    def __init__(self, config, peptide_id_mapping, pep_mod_dict, target_df, subject_df) -> None:
        super().__init__(config, peptide_id_mapping)
        self.target_df = target_df
        self.pep_mod_dict = pep_mod_dict
        self.subject_filtered = None
        self.subject_part_search = None
        self.subject_df = subject_df


    def update_dataframe(self, threshold):
        """
        Retains PSMs from target-decoy search that are mapped to spectra present in the updated target dataframe;
        adds proxy ground truth labels to the dataframe.
        """
        quality_scans = self.extract_scans_above_threshold(self.target_df, threshold)
        self.subject_filtered = self.filter_by_high_quality_scans(self.subject_df, quality_scans)
        self.generate_pep_mod_idx_for_td()
        self.subject_part_search = self.extract_partition_search_results()
        self.add_ground_truth_labels()

        return self.subject_part_search


    def generate_pep_mod_idx_for_td(self):

        pep_mods = zip(self.subject_filtered['sequence'], self.subject_filtered['modifications'])
        pep_mods_idx = list(map(lambda x: self.pep_mod_dict.get(x, -1), pep_mods))
        self.subject_filtered['pep_mod_index'] = pep_mods_idx


    def add_ground_truth_labels(self):
        gt_status = (self.subject_part_search['hit_rank'] == 1)
        self.subject_part_search['gt_status'] = gt_status


    @staticmethod
    def filter_by_high_quality_scans(df, quality_scans):
        return df[df.scan.isin(quality_scans)].copy()


    def extract_partition_search_results(self):
        mask = (self.peptide_id_mapping[self.subject_filtered['pep_mod_index']] == self.subset_id)
        return self.subject_filtered.loc[mask].groupby('scan').head(1)



class EntrapmentPartitionSearchSimulation(PartitionSearchSimulation):
    
    def __init__(self, config, peptide_id_mapping, pep_mod_dict, target_df, subject_df) -> None:
        super().__init__(config, peptide_id_mapping)
        self.target_df = target_df
        self.pep_mod_dict = pep_mod_dict
        self.subject_filtered = None
        self.subject_part_search = None
        self.subject_df = subject_df


    def update_dataframe(self, threshold):
        """
        Retains PSMs from target-decoy search that are mapped to spectra present in the updated target dataframe;
        adds proxy ground truth labels to the dataframe.
        """
        quality_scans = self.extract_scans_above_threshold(self.target_df, threshold)
        self.subject_filtered = self.filter_by_high_quality_scans(self.subject_df, quality_scans)
        self.generate_pep_mod_idx_for_td()
        self.subject_part_search = self.extract_partition_search_results()
        self.add_ground_truth_labels()

        return self.subject_part_search


    def generate_pep_mod_idx_for_td(self):

        pep_mods = zip(self.subject_filtered['sequence'], self.subject_filtered['modifications'])
        pep_mods_idx = list(map(lambda x: self.pep_mod_dict.get(x, -1), pep_mods))
        self.subject_filtered['pep_mod_index'] = pep_mods_idx


    def add_ground_truth_labels(self):
        gt_status = (self.subject_part_search['hit_rank'] == 1)
        self.subject_part_search['gt_status'] = gt_status


    @staticmethod
    def filter_by_high_quality_scans(df, quality_scans):
        return df[df.scan.isin(quality_scans)].copy()


    def extract_partition_search_results(self):
        mask = (self.peptide_id_mapping[self.subject_filtered['pep_mod_index']] == self.subset_id)
        return self.subject_filtered.loc[mask].groupby('scan').head(1)

class StdPartitionSearchSimulation(PartitionSearchSimulation):
    
    def __init__(self, config, peptide_id_mapping, pep_mod_dict, target_df, subject_df) -> None:
        super().__init__(config, peptide_id_mapping)
        self.target_df = target_df
        self.pep_mod_dict = pep_mod_dict
        self.subject_filtered = None
        self.subject_part_search = None
        self.subject_df = subject_df


    def update_dataframe(self, threshold):
        """
        Retains PSMs from target-decoy search that are mapped to spectra present in the updated target dataframe;
        adds proxy ground truth labels to the dataframe.
        """
        quality_scans = self.extract_scans_above_threshold(self.target_df, threshold)
        self.subject_filtered = self.filter_by_high_quality_scans(self.subject_df, quality_scans)
        self.generate_pep_mod_idx_for_td()
        self.subject_part_search = self.extract_partition_search_results()
        self.add_ground_truth_labels()

        return self.subject_part_search


    def generate_pep_mod_idx_for_td(self):

        pep_mods = zip(self.subject_filtered['sequence'], self.subject_filtered['modifications'])
        pep_mods_idx = list(map(lambda x: self.pep_mod_dict.get(x, -1), pep_mods))
        self.subject_filtered['pep_mod_index'] = pep_mods_idx


    def add_ground_truth_labels(self):
        gt_status = (self.subject_part_search['hit_rank'] == 1)
        self.subject_part_search['gt_status'] = gt_status


    @staticmethod
    def filter_by_high_quality_scans(df, quality_scans):
        return df[df.scan.isin(quality_scans)].copy()


    def extract_partition_search_results(self):
        mask = (self.peptide_id_mapping[self.subject_filtered['pep_mod_index']] == self.subset_id)
        return self.subject_filtered.loc[mask].groupby('scan').head(1)




class PepPartitionSearchSimulation(PartitionSearchSimulation):
    pass




class ThresholdEvaluator:

    def __init__(self, config) -> None:
        self.num_thresholds = int(config.get("validation.general", "num_thresholds").split(" ")[0].strip())
        self.threshold_score = config.get('validation.general', 'threshold_score', fallback='tev').split(" ")[0].strip()

    
    def threshold_dictionary(self):
        full_dict = get_thresholds_dictionary(self.num_thresholds)
        return full_dict[self.threshold_score]



class PostSearchSimulationFactory:
    @staticmethod
    def create_simulation(config, peptide_id_mapping, pep_mod_dict, target_df, td_df, validation_mode):
        validation_mode_name = validation_mode.capitalize() + "PartitionSearchSimulation"
        search_simulation_class = globals()[validation_mode_name]
        return search_simulation_class(config, peptide_id_mapping, pep_mod_dict, target_df, td_df)



class PostSearchRunner:

    def __init__(self, config) -> None:
        self.config = config
        self.post_partition = PostSearchPartition(config)
        self.validation_mode = config.get("validation.general", "validation_mode").split(" ")[0].strip()
        self.search_simulation_factory = PostSearchSimulationFactory()


    def run_postsearch(self, threshold, target_df, td_df=None, decoy_df=None):

        # TODO: may need to move this out to initializer to make the overall execution faster
        target_df, pep_mod_dict = self.post_partition.add_peptide_modification_index_to_target_df(target_df)
        peptide_id_mapping = self.post_partition.create_peptide_subset_id_mapping(target_df)
        
        partition_search = self.search_simulation_factory.create_simulation(
            self.config, peptide_id_mapping, pep_mod_dict, target_df, td_df, self.validation_mode
        )

        return partition_search.update_dataframe(threshold)



class PostSearchOrchestrated:

    def __init__(self, config) -> None:
        
        self.threshold_evaluator = ThresholdEvaluator(config)
        self.postsearch_runner = PostSearchRunner(config)
        self.post_partition = PostSearchPartition(config)
        fdr_calc_name = config.get("validation.general", "validation_mode").split(" ")[0].strip().capitalize() + "FdpFdrCalculation"
        print(fdr_calc_name)
        self.fdp_fdr_calculation = getattr(fdr, fdr_calc_name)(config)


    def run_postsearch_validation(self, target_df, subject_df=None, decoy_df=None):
        
        threshold_dict = self.threshold_evaluator.threshold_dictionary()
        fdr_fdp_results = []

        for threshold in tqdm(threshold_dict, desc="Evaluating thresholds"):
        
            updated_df = self.postsearch_runner.run_postsearch(threshold, target_df, subject_df)
            cur_fdr_fdp_results = self.fdp_fdr_calculation.calculate_fdp_fdr_contour(updated_df, decoy_df)
            fdr_fdp_results.append(cur_fdr_fdp_results)

        return fdr_fdp_results, threshold_dict, updated_df
    

    def run_postsearch_single_threshold(self, target_df, th_val, subject_df=None, decoy_df=None):
        
        updated_df = self.postsearch_runner.run_postsearch(th_val, target_df, subject_df)
        return updated_df


    def run_postsearch_bootstrap(self, opt_threshold, num_rep, target_df, subject_df=None, decoy_df=None):
        
        fdr_fdp_results = []
        cur_rep = 0
        while cur_rep <= num_rep:
        
            updated_df = self.postsearch_runner.run_postsearch(opt_threshold, target_df, subject_df)
            cur_fdr_fdp_results = self.fdp_fdr_calculation.calculate_fdp_fdr_contour(updated_df, decoy_df)
            fdr_fdp_results.append(cur_fdr_fdp_results)
            cur_rep += 1

        return fdr_fdp_results, opt_threshold, updated_df








class PostSearchValidation(QualityFiltering):

    """Class to execute validation by post-search partition"""

    def __init__(self, config, fdp_fdr_calculator: FdpFdrCalculator, fdr_level: FdrLevel, peptide_id_mapping, pep_mod_dict, target_df, td_df, decoy_df):

        self.fdr_score = config.get('validation.general', 'fdr_score', fallback='tev').split(" ")[0].strip()
        self.threshold_score = config.get('validation.general', 'threshold_score', fallback='tev').split(" ")[0].strip()
        self.subset_id = int(config.get('validation.general', 'subset_id', fallback=0).split(" ")[0].strip())
        self.num_thresholds = int(config.get('validation.general', 'num_thresholds', fallback=10).split(" ")[0].strip())
        self.num_partitions = float(config.get('partition.general', 'num_partitions', fallback=2).split(" ")[0].strip())

        self.target_df = target_df
        self.td_df = td_df
        self.decoy_df = decoy_df

        self.decoy_factor = float(1 / self.num_partitions)
        self.peptide_id_mapping = peptide_id_mapping
        self.pep_mod_dict = pep_mod_dict

        self.fdp_fdr_calculator = fdp_fdr_calculator
        self.fdr_level = fdr_level


    def _remove_below_threshold(self, target_df, threshold):
        """
        Filter spectra based on a score threshold in a target-only search.

        Returns:
        list: Indices of spectra passing the score threshold in the target-only search.
            """

        scans_above_thr = target_df[(target_df['hit_rank'] == 1) & (target_df[self.threshold_score] >= threshold)].scan
        scan_mask = target_df.scan.isin(scans_above_thr)
        target_df = target_df[scan_mask]
        is_in_subset_mask = self.peptide_id_mapping[target_df['pep_mod_index']] == self.subset_id

        return target_df.loc[is_in_subset_mask].groupby('scan').first().index


    def generate_pep_mod_idx_for_td(self, td_quality):

        pep_mods = zip(td_quality['sequence'], td_quality['modifications'])
        pep_mods_idx = list(map(lambda x: self.pep_mod_dict.get(x, -1), pep_mods))
        return pep_mods_idx


    def _filter_update_target_decoy_df(self, td_df, peptide_id_mapping, subset_id, quality_scans):
        """
        Retains PSMs from target-decoy search that are mapped to spectra present in the updated target dataframe;
        adds proxy ground truth labels to the dataframe.
        """


        td_quality = td_df[td_df.scan.isin(quality_scans)].copy()

        td_quality['pep_mod_index'] = self.generate_pep_mod_idx_for_td(td_quality)

        mask = (td_quality['is_decoy']) | (peptide_id_mapping[td_quality['pep_mod_index']] == subset_id)

        td_quality = td_quality.loc[mask].groupby('scan').head(1)

        td_quality['gt_status'] = (~td_quality['is_decoy']) & (td_quality['hit_rank'] == 1)

        return td_quality


    def _update_identification_status_labels(self, threshold):
    
        quality_scans = self._remove_below_threshold(self.target_df, threshold)
        updated_td_df = self._filter_update_target_decoy_df(self.td_df, self.peptide_id_mapping, self.subset_id, quality_scans)

        return updated_td_df


    def calculate_fdp_fdr_contour(self):
        """
        Generate data for a contour plot of proxy false discovery proportion (FDP) vs. false discovery rate (FDR).

        Parameters:
        - n_rep (int): Number of replicates for threshold variation.
        - fdr_score (str): Score used for calculating FDR (default: 'tev').
        - thr_score (str): Score for the quality filtering (default: 'tev').

        Returns:
        tuple: A tuple containing FDP-FDR results, threshold values, and updated target-decoy DataFrame.
        """

        threshold_dict = get_thresholds_dictionary(self.num_thresholds)
        fdr_fdp_results = []

        # t_df, td_df, _ = self._read_search_results()
        # pep_dict = self._create_peptide_subset_id_mapping(t_df)
        post_partition = PostSearchPartition(self.num_partitions)
        for threshold in tqdm(threshold_dict[self.threshold_score], desc="Evaluating thresholds"):
            self.peptide_id_mapping = post_partition.create_peptide_subset_id_mapping(self.target_df)
            new_td_df = self._update_identification_status_labels(threshold)
            new_td_df = self.adjust_for_fdr_level(new_td_df, self.fdr_level)
            fdr_fdp_results.append(self.fdp_fdr_calculator.calculate_proxy_fdp_fdr(self.fdr_score, new_td_df, self.decoy_factor))

        return fdr_fdp_results, threshold_dict[self.threshold_score]


    def adjust_for_fdr_level(self, df, fdr_level: FdrLevel):
        return fdr_level.adjust_data(df, self.fdr_score)



