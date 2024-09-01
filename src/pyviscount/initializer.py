from abc import ABC, abstractmethod
from .plot import PlotSearchSpaceAnalysis
from .sequence_processing import get_amino_acid_mass_dict
from .partition import PartitionSeqDB
from .postsearch import PostSearchPartition, PostSearchValidation
from . import parsers
from . import plot
from . import fdr
from .correction import PartitionScoreCorrection
from .fdr import FdpFdrCalculator
from .presearch import PreSearchValidation


class PlotSearchSpaceAnalysisInitializer:

    def __init__(self, config) -> None:
        self.config = config
        self.aa_dict = get_amino_acid_mass_dict()


    def initialize(self):
        return PlotSearchSpaceAnalysis(self.config, self.aa_dict)



class SeqDBPartitionInitializer:

    def __init__(self, config, protein_dict) -> None:
        self.protein_dict = protein_dict
        self.config = config
        
    def initialize(self):
        return PartitionSeqDB(self.config, protein_dict=self.protein_dict)



class PostSearchValidationInitializer:

    def __init__(self, config, fdp_fdr_calculator: FdpFdrCalculator, target_df, td_df, decoy_df, peptide_id_mapping, pep_mod_dict) -> None:

        self.target_df = target_df
        self.td_df = td_df
        self.decoy_df = decoy_df
        self.peptide_id_mapping = peptide_id_mapping
        self.config = config
        self.fdp_fdr = fdp_fdr_calculator
        self.pep_mod_dict = pep_mod_dict

        fdr_level_name = config.get("validation.general", 'fdr_level', fallback='PSM').split(" ")[0].strip().capitalize()
        self.fdr_level = getattr(fdr, fdr_level_name+"LevelFdr")


    def initialize(self):
        return PostSearchValidation(self.config, self.fdp_fdr, self.fdr_level, self.peptide_id_mapping, self.pep_mod_dict, self.target_df, self.td_df, self.decoy_df)



class ParserInitializer:

    def __init__(self, config) -> None:
        
        self.engine = config.get('validation.general', 'engine', fallback='Tide').split(" ")[0].strip()


    def initialize(self):

        try:
            return getattr(parsers, f"{self.engine}Parser")()
        except AttributeError as exc:
            raise ValueError(f"Unsupported or invalid engine: {self.engine}") from exc



class PostSearchPartitionInitializer:

    def __init__(self, config):
        self.num_partitions = int(config.get('partition.general', 'num_partitions', fallback=2).split(" ")[0].strip())
        self.config = config

    def initialize(self):
        return PostSearchPartition(self.config)
    

class PlotPostSearchResultsInitializer:

    def __init__(self, config, fdp_fdr_results, target_df, subject_df, decoy_df):

        self.output_path = config.get('plotting', 'output_path', fallback='./').split(" ")[0].strip()
        self.fdp_fdr_plotter = plot.FdpFdrPlotter(fdp_fdr_results, self.output_path)
        hist_score = config.get('validation.general', 'threshold_score', fallback='tev').split(" ")[0].strip()
        # self.hist_plotter = plot.HistogramPlotter(target_df, subject_df, decoy_df, 'kde', hist_score)
    
    def initialize(self):
        return plot.PlotValidation(self.fdp_fdr_plotter)
    

class PartitionScoreCorrectionInitializer:

    def __init__(self, config, score_correction: PartitionScoreCorrection):

        self.decoy_target_ratio = float(config.get('validation.general', 'decoy_target_ratio', fallback=1).split(" ")[0].strip())
        self.num_partitions = int(config.get('partition.general', 'num_partitions',  fallback=2).split(" ")[0].strip())
        self.score_correction = score_correction


    def initialize(self, *args, **kwargs):
        return self.score_correction(self.decoy_target_ratio, self.num_partitions, *args, **kwargs)
        


class PreSearchValidationInitializer:

    def __init__(self, config, full_target, part_t, part_td, fdp_fdr_calculator):

        self.full_target = full_target
        self.part_t = part_t
        self.part_td = part_td
        self.config = config
        self.fdp_fdr_calculator = fdp_fdr_calculator

    def initialize(self):
        return PreSearchValidation(self.config, self.full_target, self.part_t, self.part_td, self.fdp_fdr_calculator)


