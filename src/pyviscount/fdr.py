from abc import ABC, abstractmethod
from typing import Tuple
import pandas as pd
import numpy as np
from scipy.interpolate import CubicSpline



class FdrLevel(ABC):

    @staticmethod
    @abstractmethod
    def adjust_data(df: pd.DataFrame, fdr_score: str):
        pass


class PsmLevelFdr(FdrLevel):

    @staticmethod
    def adjust_data(df, fdr_score):
        return df


class PeptideLevelFdr(FdrLevel):

    @staticmethod
    def adjust_data(df, fdr_score):
        peptide_level_df = df.groupby(['sequence', 'modifications']).max(fdr_score)
        peptide_level_df.reset_index(inplace=True)
        return peptide_level_df






class FdpFdrCalculator(ABC):

    @staticmethod
    @abstractmethod
    def calculate_proxy_fdp_fdr() -> Tuple:
        """Calculating proxy FDP and FDR estimates"""
        pass



class DecoyCountingCalculator(FdpFdrCalculator):


    @staticmethod
    def calculate_proxy_fdp_fdr(df: pd.DataFrame, fdr_score: str, decoy_factor: float = 1.0) -> Tuple:
        """
        Processes the data to obtain proxy FDP values and decoy-based FDR estimates.

        Parameters:
        - dataframe (pd.DataFrame): Input DataFrame with filtered target-decoy search results.
        - fdr_score_column (str): Name of the column representing the score used for FDR estimation.

        Returns:
        tuple: A tuple containing arrays of FDR estimates and corresponding proxy FDP values.
        """

        if ("is_decoy" not in df.columns) or ("gt_status" not in df.columns):
            raise ValueError("The dataframe provided to this FDR calculator\
                             must have 'is_decoy' and 'gt_status' columns.")


        df.sort_values(fdr_score, ascending=False, inplace=True)

      
        dec_cumsum = decoy_factor * (df["is_decoy"].cumsum())
        df["cum_dec"] =  dec_cumsum / (~df['is_decoy']).cumsum()

        # add column with FDP values
        df_no_decoys = df[~df["is_decoy"]].copy()
        no_decoys_index = np.arange(1, len(df_no_decoys)+1)
        df_no_decoys["cum_neg"] = abs(1-df_no_decoys['gt_status']).cumsum() / no_decoys_index
        fdr_dec = df_no_decoys["cum_dec"].to_numpy()
        fdp = df_no_decoys["cum_neg"].to_numpy()

        return (fdr_dec, fdp)
    

class EntrapmentCountingCalculator(FdpFdrCalculator):


    @staticmethod
    def calculate_proxy_fdp_fdr(df: pd.DataFrame, fdr_score: str, entrapment_df: pd.DataFrame = None) -> Tuple:
        """
        Processes the data to obtain proxy FDP values and entrapment query-based FDR estimates.

        Parameters:
        - dataframe (pd.DataFrame): Input DataFrame with filtered target-decoy search results.
        - fdr_score_column (str): Name of the column representing the score used for FDR estimation.

        Returns:
        tuple: A tuple containing arrays of FDR estimates and corresponding proxy FDP values.
        """


        def get_pi0(pvals, final_lambda=0.95):
            pi0s = []
            lambdas = np.linspace(0.001, 0.99, 100)


            for l in lambdas:
                pi0_cur = len(pvals[pvals > l]) / (len(pvals) * (1-l))
                pi0s.append(pi0_cur)
                

            cs = CubicSpline(lambdas, pi0s)
       
            return cs([final_lambda])[0]



        def new_ent_mixmax(ws, zs, pi0, df, fdr_score):

            nwz = np.sum(df[fdr_score].values[:, np.newaxis] <= zs, axis=0)
            nzz = np.sum(zs[:, np.newaxis] <= zs, axis=0)
            nzw = np.sum(zs[:, np.newaxis] > ws, axis=0)
            nww = np.sum(df[fdr_score].values[:, np.newaxis] >= ws, axis=0)

            n = len(ws)
            n_e = len(zs)

            fg = n / n_e

            e1 = 0
            j = n_e - 1

            rs = np.zeros(n)

            for m in range(n-1, 0, -1):

                while (j > 0) and (zs[j] >= ws[m]):

                    p_est = (nwz[j] / ((1-pi0) * nzz[j])) - ((fg * pi0) / (1 - pi0))
                    e1 += min(1, max(0, p_est)) * (1- pi0)
                    j -= 1
                
                rs[m] = ( pi0 * fg * nzw[m] + fg * e1) / nww[m]
                rs[m] = min(rs[m], 1)

            return rs

        df.sort_values(fdr_score, ascending=True, inplace=True)

        ws = df[fdr_score].values
        zs = np.array(entrapment_df[fdr_score].values)
        pvals = np.array([(len(zs[zs >= x])+1) / (len(zs) + 1) for x in ws[::-1]])
        pvals = pvals[pvals >= 0]
        pvals = pvals[pvals <= 1]
        pi0_est = get_pi0(pvals)

        #pi0_est = df['gt_status'].value_counts().get(0, 0) / max(1,len(df))

        fdrs = new_ent_mixmax(ws, zs, pi0_est, df, fdr_score)[1:]

        df.sort_values(fdr_score, ascending=False, inplace=True)
        df['fdp'] = (~df['gt_status']).cumsum() / np.arange(1, len(df)+1)

        fdps = df['fdp'][::-1][1:].to_numpy()

        return (fdrs, fdps)




class BhFdpFdrCalculator(FdpFdrCalculator):
    """Estimate FDR using p-values in adaptive BH framework"""

    @staticmethod
    def calculate_proxy_fdp_fdr(df_sorted: pd.DataFrame, p_value_column: str, *args):

        df_labels = df_sorted['gt_status'].to_numpy()
        sorted_pvals = df_sorted[p_value_column].to_numpy()
        length_df = len(df_labels)
        fdr_threshold_array = np.linspace(0.001, 0.1, 100)

        pi_zero = df_sorted['gt_status'].value_counts().get(0, 0) / max(1,len(df_sorted))
        critical_vals = (1/ pi_zero) * np.arange(1, length_df + 1) / length_df
        critical_array = [critical_vals * x for x in fdr_threshold_array]
       
        # len_df_correct_labels = len(df_labels[df_labels == pos_label])
        masks = [sorted_pvals <= x for x in critical_array]
        bh_gt_labels = [df_labels[mask] for mask in masks]

        fdps = np.array([len(x[x == False]) / max(1,len(x)) for x in bh_gt_labels])
        # tprs = [len(x[x == pos_label]) / len_df_correct_labels for x in bh_gt_labels]
        return fdr_threshold_array, fdps






class FdpFdrCalculation(ABC):

    def __init__(self, config) -> None:

        self.fdr_score = config.get('validation.general', 'fdr_score', fallback='tev').split(" ")[0].strip()


    @abstractmethod
    def calculate_fdp_fdr_contour(self, df):
        pass


    def adjust_for_fdr_level(self, df: pd.DataFrame, fdr_level: FdrLevel):
        return fdr_level.adjust_data(df, self.fdr_score)



class EntrapmentFdpFdrCalculation(FdpFdrCalculation):
    """FDR calculated using separate target-decoy search results"""
    
    def __init__(self, config) -> None:
        super().__init__(config)
        fdr_level_name = config.get("validation.general", "fdr_level").split(" ")[0].strip().capitalize() + "LevelFdr"
        self.fdr_level = globals()[fdr_level_name]
        # num_partitions = float(config.get("partition.general", "num_partitions", fallback=2).split(" ")[0].strip())
        # partition_mode = config.get("validation.general", "partition_mode").split(" ")[0].strip()
        self.fdr_calculator = EntrapmentCountingCalculator


    def calculate_fdp_fdr_contour(self, subject_df, decoy_df=None):


        filtered_subject_df = self.adjust_for_fdr_level(subject_df, self.fdr_level)
        filtered_ent_df = self.adjust_for_fdr_level(decoy_df, self.fdr_level)
    
        fdr_fdp_results = self.fdr_calculator.calculate_proxy_fdp_fdr(filtered_subject_df, self.fdr_score, filtered_ent_df)

        return fdr_fdp_results
    

class StdFdpFdrCalculation(FdpFdrCalculation):
    """FDR calculated using separate target-decoy search results"""
    
    def __init__(self, config) -> None:
        super().__init__(config)
        fdr_level_name = config.get("validation.general", "fdr_level").split(" ")[0].strip().capitalize() + "LevelFdr"
        self.fdr_level = globals()[fdr_level_name]
        num_partitions = float(config.get("partition.general", "num_partitions", fallback=2).split(" ")[0].strip())
        decoy_target_ratio = float(config.get("validation.general", "decoy_target_ratio", fallback=1).split(" ")[0].strip())
        self.decoy_factor = float(1 / decoy_target_ratio)
        self.decoy_mode = config.get("validation.general", "decoy_mode", fallback='Counting').split(" ")[0].strip().capitalize()
        partition_mode = config.get("validation.general", "partition_mode").split(" ")[0].strip()

        if partition_mode == 'postsearch':
            self.decoy_factor *= float(1 / num_partitions)

        if self.decoy_mode == 'P-value':
            self.fdr_calculator = BhFdpFdrCalculator
        else:
            self.fdr_calculator = DecoyCountingCalculator


    def calculate_fdp_fdr_contour(self, subject_df, decoy_df=None):

        if self.fdr_calculator == DecoyCountingCalculator:
            real_fdr_score = self.fdr_score
            processed_df = self.append_decoy_df(subject_df, decoy_df)
        else:
            real_fdr_score = 'Decoy_p_value'
            processed_df = subject_df

        #filtered_df = self.adjust_for_fdr_level(processed_df, self.fdr_level)
        fdr_fdp_results = self.fdr_calculator.calculate_proxy_fdp_fdr(processed_df, real_fdr_score, self.decoy_factor)

        return fdr_fdp_results
    
    @staticmethod
    def append_decoy_df(subject_df, decoy_df):

        top_decoys = decoy_df[decoy_df['hit_rank'] == 1].copy()
        concat_df = pd.concat([subject_df, top_decoys], join='outer', ignore_index=True)
        return concat_df



class TdcFdpFdrCalculation(FdpFdrCalculation):

    def __init__(self, config) -> None:
        super().__init__(config)
        fdr_level_name = config.get("validation.general", "fdr_level").split(" ")[0].strip().capitalize() + "LevelFdr"
        self.fdr_level = globals()[fdr_level_name]
        num_partitions = float(config.get("partition.general", "num_partitions", fallback=2).split(" ")[0].strip())
        decoy_target_ratio = float(config.get("validation.general", "decoy_target_ratio", fallback=1).split(" ")[0].strip())
        self.decoy_factor = float(1 / decoy_target_ratio)
        self.decoy_mode = config.get("validation.general", "decoy_mode", fallback='Counting').split(" ")[0].strip().capitalize()
        partition_mode = config.get("validation.general", "partition_mode").split(" ")[0].strip()

        if partition_mode == 'postsearch':
            self.decoy_factor *= float(1 / num_partitions)

        if self.decoy_mode == 'p-value':
            self.fdr_calculator = BhFdpFdrCalculator
        else:
            self.fdr_calculator = DecoyCountingCalculator


    def calculate_fdp_fdr_contour(self, td_df, *args, **kwargs):
     
        filtered_df = self.adjust_for_fdr_level(td_df, self.fdr_level)
        fdr_fdp_results = self.fdr_calculator.calculate_proxy_fdp_fdr(filtered_df, self.fdr_score, self.decoy_factor)

        return fdr_fdp_results



class DecoyfreeFdpFdrCalculation(FdpFdrCalculation):
    
    def __init__(self, config) -> None:
        super().__init__(config)
        self.fdr_calculator = BhFdpFdrCalculator
        fdr_level_name = config.get("validation.general", "fdr_level").split(" ")[0].strip().capitalize() + "LevelFdr"
        self.fdr_level = globals()[fdr_level_name]
        self.p_value_name = config.get("validation.general", "fdr_score").split(" ")[0].strip()


    def calculate_fdp_fdr_contour(self, target_df, *args, **kwargs):
     
        filtered_df = self.adjust_for_fdr_level(target_df, self.fdr_level)
        filtered_df.sort_values(self.p_value_name, ascending=True, inplace=True)
        fdr_fdp_results = self.fdr_calculator.calculate_proxy_fdp_fdr(filtered_df, self.p_value_name)

        return fdr_fdp_results



class PepFdpFdrCalculation(FdpFdrCalculation):
    pass


