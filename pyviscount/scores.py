from abc import ABC, abstractmethod
import scipy.stats as st
import pandas as pd
from .parsers import ParamFileParser, PercolatorParser
from .correction import TEVPartitionCorrection, SidakCorrectionMixin


class ScoreProcessing:
    """Based on user settings, correct all necessary scores"""

    def __init__(self, config) -> None:
        
        self.config = config
        self.validation_mode = config.get('validation.general', 'validation_mode').split(" ")[0].strip().capitalize()
        possible_validation_modes = ['Decoyfree', 'Tdc', 'Std', 'Entrapment']
        if self.validation_mode not in possible_validation_modes:
            raise ValueError(f"The parameter: validation_mode needs to be one of the following:\
                             {possible_validation_modes}")

        self.decoy_target_ratio = float(config.get('validation.general', 'decoy_target_ratio', fallback=1).split(" ")[0].strip())
        self.num_partitions = int(config.get('partition.general', 'num_partitions',  fallback=2).split(" ")[0].strip())

        all_postprocessor_entries = [item for item in config.sections() if item.startswith('validation.postprocessor_score')]
        self.postprocessor_entries = [entry for entry in all_postprocessor_entries if config.get(entry, 'add_score') == 'true']

        self.sidak_status = config.get('validation.general', 'sidak_correction').split(" ")[0].strip()
        self.decoy_mode = config.get('validation.general', 'decoy_mode').split(" ")[0].strip()

    def add_postprocessor_scores(self, df):

        new_df = df.copy()

        for entry in self.postprocessor_entries:
            postprocessor_name = entry.split('.')[-1].capitalize()
            score_adder = globals()[f"{postprocessor_name}ScoreAdder"]
            new_df = score_adder(self.config).add_score(new_df)

        return new_df


    def correct_scores(self, df):

        df = self.correct_tev(df)

        if self.sidak_status == 'true':
            df = self.correct_sidak(df)

        return df
    

    def add_std_decoy_pval(self, target_df, decoy_df):

        new_target_df = DecoyScoreAdder(self.config).add_score(target_df, decoy_df)
        return new_target_df


    def correct_tev(self, df):
        
        tev_correction = TEVPartitionCorrection(self.decoy_target_ratio, self.num_partitions)

        if self.validation_mode == 'Tdc':
            df = tev_correction.correct_tdc(df)
        elif self.validation_mode == 'Std':
            df = tev_correction.correct_std(df)
        elif self.validation_mode == 'Decoyfree':
            df = tev_correction.correct_decoy_free(df)
        elif self.validation_mode == 'Entrapment':
            df = tev_correction.correct_decoy_free(df)
        else:
            raise ValueError("The validation_mode parameter is not valid.")
        
        return df


    def correct_sidak(self, df):

        sidak_correction = SidakCorrectionMixin(self.config)
        return sidak_correction.correct_pvals(df)



class PostProcessorScoreAdder(ABC):

    @abstractmethod
    def add_score(self, df_to_modify, *args):
        pass


class PercolatorScoreAdder(PostProcessorScoreAdder):

    def __init__(self, config) -> None:

        percolator_results_path = config.get("validation.postprocessor_score.percolator", "file_path").split(" ")[0].strip()
        self.percolator_results = PercolatorParser(percolator_results_path).parse()
        self.percolator_score_name = config.get('validation.postprocessor_score.pylord', 'score_name').split(" ")[0].strip()


    def add_score(self, df_to_modify, scan_col='scan'):

        if scan_col not in df_to_modify.columns:
            raise ValueError(f"There is no {scan_col} column in the dataframe;\
                             cannot add Percolator score.")

        percolator_scans = self.extract_scan_list()
        percolator_scores = self.extract_percolator_scores()

        df_to_modify.index = df_to_modify.scan
        df_to_modify.loc[percolator_scans, 'percolator_score'] = percolator_scores
        df_to_modify.reset_index(inplace=True, drop=True)

        # Percolator may return scores for less PSMs than originally present in the dataset
        df_to_modify.dropna(inplace=True)

        return df_to_modify

    
    def extract_percolator_scores(self):

        return [x.score for x in self.percolator_results]

    def extract_scan_list(self):

        return [float(x.psm_id.split('_')[-3]) for x in self.percolator_results]


class PylordScoreAdder(PostProcessorScoreAdder):
    
    def __init__(self, config):
        params_path = config.get('validation.postprocessor_score.pylord', 'file_path').split(" ")[0].strip()
        self.params_dict = ParamFileParser(params_path).parse()
        self.pv_column_name = config.get('validation.postprocessor_score.pylord', 'score_name').split(" ")[0].strip()

    
    def add_score(self, df_to_modify):

        group = df_to_modify.groupby('charge')['tev']
        for label, items in group:
            df_to_modify.loc[items.index, self.pv_column_name] = 1 - st.gumbel_r.cdf(items.values, *self.params_dict[label])
        df_to_modify[self.pv_column_name].replace(0, 1e-16, inplace=True)
        return df_to_modify


class CddScoreAdder(PostProcessorScoreAdder):

    def __init__(self, config):
        params_path = config.get('validation.postprocessor_score.cdd', 'file_path').split(" ")[0].strip()
        self.pv_column_name = config.get('validation.postprocessor_score.cdd', 'score_name').split(" ")[0].strip()
        self.params_dict = ParamFileParser(params_path).parse()

    def add_score(self, df_to_modify):

        #group = df_to_modify[df_to_modify['hit_rank'] == 1].groupby('charge')['tev']
        group = df_to_modify.groupby('charge')['tev']
        for label, items in group:
            items_pvals = 1 - st.gumbel_r.cdf(items.values, *self.params_dict[label])
            df_to_modify.loc[items.index, self.pv_column_name] = items_pvals
        df_to_modify[self.pv_column_name].replace(0, 1e-16, inplace=True)
        return df_to_modify


class DecoyScoreAdder(PostProcessorScoreAdder):
    """Adding decoy-derived p-values"""

    def __init__(self, config) -> None:
        self.threshold_score = config.get('validation.general', 'threshold_score').split(" ")[0].strip()


    def add_score(self, target_df, decoy_df):

        pv_column_name = 'Decoy_p_value'
        return self.calculate_decoy_pvals(target_df, decoy_df, pv_column_name)


    def calculate_decoy_pvals(self, target_df, decoy_df, pv_column_name):
        
        td_df = self.concatenate_dfs(target_df, self.select_top_hits(decoy_df))
        td_df.sort_values(self.threshold_score, ascending=False, inplace=True)

        # calculate decoy p-values as outlined in Granholm 2011s
        td_df['cum_dec'] = td_df['is_decoy'].cumsum()
        td_df[pv_column_name] = (td_df['cum_dec'] + 1) / (len(td_df[td_df['is_decoy']]) + 1)
        target_psms = td_df[~td_df['is_decoy']].copy()

        return target_psms
    
    @staticmethod
    def concatenate_dfs(*args):
        return pd.concat(args, ignore_index=True)

    @staticmethod
    def select_top_hits(df):
        return df[df['hit_rank'] == 1]


class PeptideprophetScoreAdder(PostProcessorScoreAdder):
    pass


class IprophetScoreAdder(PostProcessorScoreAdder):
    pass

