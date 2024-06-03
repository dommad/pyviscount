from abc import ABC, abstractmethod
import numpy as np
from .constants import TH_BETA




class PartitionScoreCorrection(ABC):
    
    @abstractmethod
    def correct_tdc(self, *args, **kwargs):
        pass

    @abstractmethod
    def correct_decoy_free(self, *args, **kwargs):
        pass

    @abstractmethod
    def correct_std(self, *args, **kwargs):
        pass



class TEVPartitionCorrection(PartitionScoreCorrection):

    def __init__(self, decoy_target_ratio, num_partitions):
        self.decoy_target_ratio = decoy_target_ratio
        self.num_partitions = num_partitions


    def correct_tdc(self, td_df):
        """
        add decoy adjustment to account for the fact that
        targets from TD dataframe come from searching subset target + full decoy
        so we need a correction factor for TEV
        """
        adjustment_factor = 1 - 1 / (self.num_partitions * (1 + self.decoy_target_ratio))
        td_df['tev'] -= TH_BETA * np.log(adjustment_factor)

        return td_df


    def correct_decoy_free(self, subject_df):
        
        adjustment_factor = 1 - 1 / self.num_partitions
        subject_df['tev'] -= TH_BETA * np.log(adjustment_factor)

        return subject_df


    def correct_std(self, subject_df):
        """Correct the TEV score for targets obtained from separate target-decoy search"""

        adjustment_factor = 1 - 1 / self.num_partitions
        subject_df['tev'] -= TH_BETA * np.log(adjustment_factor)

        return subject_df
    


class SidakCorrectionMixin:
    """Sidak correction of user-selected p-values"""

    def __init__(self, config) -> None:
        self.num_partitions = int(config.get('partition.general', 'num_partitions', fallback=2).split(" ")[0].strip())
        self.sidak_pval_name = config.get('validation.general', 'sidak_p_value_name').split(" ")[0].strip()

    def correct_pvals(self, df):
        df['num_candidates'].replace(0, 1, inplace=True)
        df[f"Sidak_{self.sidak_pval_name}"] = 1 - pow(1 - df[self.sidak_pval_name].values, df["num_candidates"].values / self.num_partitions)
        #df[f"sidak_{p_val_column}"] = df[p_val_column].values * df["num_candidates"].values
        return df
    





