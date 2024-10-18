"""Process the pre-search partition data"""
import pandas as pd
from .utils import get_thresholds_dictionary
from .template import QualityFiltering
from .fdr import FdpFdrCalculator
from . import fdr


class PreSearchValidation(QualityFiltering):
    """Process the data from pepxml files to get the proxy FDP vs. FDR results"""

    def __init__(self, config, full_target_df, subset_df, subset_td_df):
        
        self.t_df = full_target_df
        self.part_a_df = subset_df
        self.part_a_td_df = subset_td_df

        self.fdr_score = config.get('validation.general', 'fdr_score', fallback='tev').split(" ")[0].strip()
        self.threshold_score = config.get('validation.general', 'threshold_score', fallback='tev').split(" ")[0].strip()
        self.subset_id = int(config.get('validation.general', 'subset_id', fallback=0).split(" ")[0].strip())
        self.n_rep = int(config.get('validation.general', 'num_thresholds', fallback=10).split(" ")[0].strip())
        self.num_partitions = float(config.get('partition.general', 'num_partitions', fallback=2).split(" ")[0].strip())

        self.decoy_factor = 1 # because target-decoy search is done on partitoned space already

        self.part_a_df.set_index('scan', inplace=True)
        self.part_a_td_df.set_index('scan', inplace=True)

        fdr_calc_name = config.get("validation.general", "validation_mode").split(" ")[0].strip().capitalize() + "FdpFdrCalculation"
        self.fdp_fdr_calculation = getattr(fdr, fdr_calc_name)(config)



    def _remove_below_threshold(self, t_df, threshold):
        """Truncate the dataframe based on selected threshold
        use either the similarity score threshold or PeptideProphet's probability
        to filter out mislabeled ground truth"""
        updated_df = t_df.copy()
        return updated_df[updated_df[self.threshold_score] > threshold]


    @staticmethod
    def index_peptide_modifications(df: pd.DataFrame):

        df['pep_mod'] = list(zip(df['sequence'], df['modifications']))
        df['pep_mod_idx'], pep_mods = df['pep_mod'].factorize()

        return df, dict(zip(pep_mods, df['pep_mod_idx'].values))
    

    @staticmethod
    def index_peptide_modifications_subset(df: pd.DataFrame, pep_mod_dict: dict):
        
        pep_mods = list(zip(df['sequence'], df['modifications']))
        df['pep_mod_idx'] = [pep_mod_dict.get(x) if pep_mod_dict.get(x) is not None else -1 for x in pep_mods]

        return df


    def _update_identification_status_labels(self, threshold):
        """Combine and compare the search results of target-only set A, set A+B,
        target-decoy set A"""

        filtered_t_df = self._remove_below_threshold(self.t_df, threshold)
        filtered_t_df.set_index('scan', inplace=True)
        updated_part_a_df, filtered_t_df, is_shared_scan = self._update_partition_df_shared_scans(self.part_a_df, filtered_t_df)
        
        positive_scans = self._find_positive_scans(filtered_t_df, updated_part_a_df)
        #return filtered_t_df, updated_part_a_df, positive_scans, pep_mod_dict
        updated_td_df = self.part_a_td_df.loc[is_shared_scan,:].copy()
        td_df = self._set_ground_truth_labels(updated_td_df, positive_scans)

        return td_df
        #return filtered_t_df, updated_part_a_df

        

    def _find_positive_scans(self, filtered_t_df, part_a_df):
        """Getting scans of positive samples, i.e., those
        whose top candidate is the same for both full and partition search"""
        
        #print(part_a_df.index.values)
        #print(filtered_t_df.index.values)
        #positive_mask = filtered_t_df.loc[part_a_df.index.values, 'pep_mod_idx'].values == part_a_df['pep_mod_idx'].values
        filtered_t_df.sort_index(inplace=True)
        part_a_df.sort_index(inplace=True)
        positive_mask = (filtered_t_df['sequence'] == part_a_df['sequence'])
        return part_a_df.index[positive_mask]


    def _set_ground_truth_labels(self, td_df, positive_scans):

        td_df['gt_status'] = False
        td_df.loc[positive_scans, 'gt_status'] = True
        return td_df


    def _update_partition_df_shared_scans(self, part_a_df, filtered_t_df):
        """Make sure partition df and full search space df share the same scans,
        avoid potentially confusing matches that were make it full"""

        is_shared_scan = list(set.intersection(set(part_a_df.index), set(filtered_t_df.index)))
        return part_a_df.loc[is_shared_scan, :].copy(), filtered_t_df.loc[is_shared_scan, :].copy(), is_shared_scan


    def _calculate_proxy_fdp_fdr(self, threshold):
        td_df_updated = self._update_identification_status_labels(threshold)
        return self.fdp_fdr_calculation.calculate_fdp_fdr_contour(td_df_updated)
        #return self.fdr_calculator.calculate_proxy_fdp_fdr(self.fdr_score, td_df_updated, self.decoy_factor)


    def calculate_fdp_fdr_contour(self):
        """Calculate proxy FDP vs. FDR values for different
        threshold values to produce an input to contour plot"""

        threshold_dict = self._get_threshold_dict()
        fdr_fdp_results = [self._calculate_proxy_fdp_fdr(threshold) for threshold in threshold_dict]
        return fdr_fdp_results, threshold_dict


    def _get_threshold_dict(self):
        return get_thresholds_dictionary(self.n_rep)[self.threshold_score]




# class StableRegionFinder:


#     @staticmethod
#     def find_nearest(array, value):
#         """find the nearest value in an array"""
#         idx = (np.abs(array - value)).argmin()
#         return array[idx], idx


#     def get_partial(self, arr, thresholds, fdr_th):
#         """Calculate approximate values of partial derivatives"""

#         derivatives = []
#         thres = []

#         # instead of single points, get a range around the
#         # given fdr threshold and calculate this for all of them, then average
#         # otherwise the data is very spikey
#         for idx in range(1, len(arr)-1):
#             cur_fdrs, cur_fdps = arr[idx]
#             next_fdrs, next_fdps = arr[idx+1]
#             prior_fdrs, prior_fdps = arr[idx-1]

#             cur_th = thresholds[idx]
#             next_th = thresholds[idx+1]
#             prior_th = thresholds[idx-1]

#             fdr_m, idx_m = self.find_nearest(cur_fdrs, fdr_th)
#             fdr_u, idx_u = self.find_nearest(next_fdrs, fdr_th)
#             fdr_l, idx_l = self.find_nearest(prior_fdrs, fdr_th)
#             fdp_m = cur_fdps[idx_m]
#             fdp_l = prior_fdps[idx_l]
#             fdp_u = next_fdps[idx_u]

#             derivative_1 = abs((fdp_u - fdp_m) / (fdr_u - fdr_m))
#             derivative_2 = abs((fdp_m - fdp_l) / (fdr_m - fdr_l))
#             derivatives.append(np.mean([derivative_1, derivative_2]))
#             thres.append((cur_th))

#         return derivatives, thres

