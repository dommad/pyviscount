import logging
import pandas as pd
import numpy as np
from functools import partial
from .parsers import FastaParser
from . import sequence_processing as sp
from . import searchspace_analysis as sa
from . import initializer as init
from .utils import open_config
from .template import FileReading, AveragedFdpFdr, ConfidenceInterval
from .fdr import DecoyCountingCalculator, BhFdpFdrCalculator
from .correction import TEVPartitionCorrection, SidakCorrectionMixin
from .decoy_generation import ShuffledDecoy, PeptideDecoy
from .exporter import PeptideExporter
from .postsearch import PostSearchOrchestrated, ThresholdEvaluator
from .scores import ScoreProcessing
from . import optimal_selection as osel
import matplotlib.pyplot as plt



def run_presearch_validation_new(config_file_path, full_target_file, partition_target_file, partition_td_file):

    config = open_config(config_file_path)

    logging.info("Started reading the files...")

    read_files = partial(_read_search_files, config)
    #add_percolator_score = partial(_add_percolator_score, config)
    #calculate_fdp_fdr = partial(_calculate_fdp_fdr, config)

    full_df, part_t_df, part_td_df = read_files(full_target_file, partition_target_file, partition_td_file)
    #full_df = add_percolator_score(full_df)

    fdp_fdr_calculator = DecoyCountingCalculator
    validation_object = init.PreSearchValidationInitializer(config, full_df, part_t_df, part_td_df, fdp_fdr_calculator).initialize()
    fdp_fdr_results, thresholds = validation_object.calculate_fdp_fdr_contour()


    #fdp_fdr_results, thresholds = calculate_fdp_fdr(full_df, part_t_df, part_td_df)

    logging.info("Starting calculation of proxy FDP vs. FDR contour...")
    plot_postsearch_results(config, fdp_fdr_results, full_df, part_t_df, part_td_df)
    logging.info("Finished the analysis!")
    return fdp_fdr_results, full_df


def run_presearch_validation_single(config_file_path, full_target_file, partition_target_file, partition_td_file):

    config = open_config(config_file_path)

    logging.info("Started reading the files...")

    read_files = partial(_read_search_files, config)
    #add_percolator_score = partial(_add_percolator_score, config)
    #calculate_fdp_fdr = partial(_calculate_fdp_fdr, config)

    full_df, part_t_df, part_td_df = read_files(full_target_file, partition_target_file, partition_td_file)
    #full_df = add_percolator_score(full_df)

    fdp_fdr_calculator = DecoyCountingCalculator
    validation_object = init.PreSearchValidationInitializer(config, full_df, part_t_df, part_td_df, fdp_fdr_calculator).initialize()
    fdp_fdr_results, thresholds = validation_object.calculate_fdp_fdr_contour()


    #fdp_fdr_results, thresholds = calculate_fdp_fdr(full_df, part_t_df, part_td_df)

    logging.info("Starting calculation of proxy FDP vs. FDR contour...")
    plot_postsearch_results(config, fdp_fdr_results, full_df, part_t_df, part_td_df)
    logging.info("Finished the analysis!")
    return fdp_fdr_results, full_df


def run_presearch_validation(config_file_path, full_target_file, partition_target_file, partition_td_file):

    config = open_config(config_file_path)

    logging.info("Started reading the files...")

    read_files = partial(_read_search_files, config)
    #add_percolator_score = partial(_add_percolator_score, config)
    calculate_fdp_fdr = partial(_calculate_fdp_fdr, config)

    full_df, part_t_df, part_td_df = read_files(full_target_file, partition_target_file, partition_td_file)
    # add Percolator score for quality filtering stage
    #full_df = add_percolator_score(full_df)
    #part_t_df = add_percolator_score(part_t_df)
    fdp_fdr_results, thresholds = calculate_fdp_fdr(full_df, part_t_df, part_td_df)

    logging.info("Starting calculation of proxy FDP vs. FDR contour...")
    plot_postsearch_results(config, fdp_fdr_results, full_df, part_t_df, part_td_df)
    logging.info("Finished the analysis!")

    return fdp_fdr_results


def _read_search_files(config, *args):
    parser_object = init.ParserInitializer(config).initialize()
    return FileReading(parser_object).read_search_results(*args)



def _calculate_fdp_fdr(config, full_df, part_t_df, part_td_df):
    fdp_fdr_calculator = DecoyCountingCalculator
    validation_object = init.PreSearchValidationInitializer(config, full_df, part_t_df, part_td_df, fdp_fdr_calculator).initialize()
    return validation_object.calculate_fdp_fdr_contour()


def run_postsearch_validation_single(config_file_path, target_file, subject_file=None, decoy_file=None):

    config = open_config(config_file_path)

    logging.info("Started reading the files...")
    parser_object = init.ParserInitializer(config).initialize()
    if subject_file is None:
        subject_file = target_file
    target_df, subject_df, decoy_df = FileReading(parser_object).read_search_results(target_file, subject_file, decoy_file)
    logging.info("Finished reading the files.")

    # section dedicated to postprocessor score addition
    #if decoyfree, then subject_df = target_df
    #if STD, then subject_df = target_df
    #if TDC, then subject_df = td_df
    score_processing = ScoreProcessing(config)
    #subject_df = score_processing.correct_scores(subject_df)
    #subject_df = score_processing.add_postprocessor_scores(subject_df)
    

    if score_processing.sidak_status == 'true':
        target_df = score_processing.correct_sidak(target_df)
        subject_df = score_processing.correct_sidak(subject_df)

    if (score_processing.validation_mode == 'Std') and (score_processing.decoy_mode == 'P-value'):
        target_df = score_processing.add_std_decoy_pval(target_df, decoy_df)
        subject_df = score_processing.add_std_decoy_pval(subject_df, decoy_df)


    logging.info("Corrected search space-dependent scores and added post-processor scores (if provided).")
    logging.info("Starting calculation of proxy FDP vs. FDR contour...")

    th_val = float(config.get("validation.extra", "tev_threshold").split(" ")[0].strip())
    updated_df = PostSearchOrchestrated(config).run_postsearch_single_threshold(target_df, th_val, subject_df, decoy_df)
    return updated_df



def run_postsearch_validation_new(config_file_path, target_file, subject_file=None, decoy_file=None):

    config = open_config(config_file_path)

    logging.info("Started reading the files...")
    parser_object = init.ParserInitializer(config).initialize()
    if subject_file is None:
        subject_file = target_file
    target_df, subject_df, decoy_df = FileReading(parser_object).read_search_results(target_file, subject_file, decoy_file)
    logging.info("Finished reading the files.")

    # section dedicated to postprocessor score addition
    #if decoyfree, then subject_df = target_df
    #if STD, then subject_df = target_df
    #if TDC, then subject_df = td_df
    score_processing = ScoreProcessing(config)
    subject_df = score_processing.correct_scores(subject_df)
    subject_df = score_processing.add_postprocessor_scores(subject_df)
    

    if score_processing.sidak_status == 'true':
        target_df = score_processing.correct_sidak(target_df)
        subject_df = score_processing.correct_sidak(subject_df)

    if (score_processing.validation_mode == 'Std') and (score_processing.decoy_mode == 'P-value'):
        target_df = score_processing.add_std_decoy_pval(target_df, decoy_df)
        subject_df = score_processing.add_std_decoy_pval(subject_df, decoy_df)


    logging.info("Corrected search space-dependent scores and added post-processor scores (if provided).")
    logging.info("Starting calculation of proxy FDP vs. FDR contour...")
    fdp_fdr_results, _, updated_df = PostSearchOrchestrated(config).run_postsearch_validation(target_df, subject_df, decoy_df)

    plot_postsearch_results(config, fdp_fdr_results, target_df, subject_df, decoy_df)

    opt_threshold = get_optimal_threshold(config, fdp_fdr_results)
    num_rep = 100
    bootstrap_results, _, updated_df = PostSearchOrchestrated(config).run_postsearch_bootstrap(opt_threshold, num_rep, target_df, subject_df, decoy_df)
    plot_optimal_fdp_fdr(bootstrap_results)
    logging.info("Finished the analysis!")
    return fdp_fdr_results, bootstrap_results


def get_optimal_threshold(config, fdp_fdr_results):

    thresholds = ThresholdEvaluator(config).threshold_dictionary()

    oq = osel.OptimalQualityThresholdFinder(fdp_fdr_results)
    norm_der_means, opt_threshold = oq.run(thresholds)

    # abstract the number of thresholds out by taking 10% of the number of original thresholds
    osel.Visualization().run_plots(norm_der_means, thresholds, 5, opt_threshold)
    return opt_threshold


def plot_optimal_fdp_fdr(fdp_fdr_results):

    fdrs = [i[0] for i in fdp_fdr_results]
    fdps = [i[1] for i in fdp_fdr_results]
    fdr_bins = np.linspace(0.001, 0.1, 100)

    bin_results = AveragedFdpFdr(fdr_bins).put_in_bins(fdrs, fdps)
    cis_to_plot = ConfidenceInterval().get_mean_cis(bin_results, 0, 1, 0.68)
    plot_mean_cis_fdp_fdr(fdr_bins, cis_to_plot)


def plot_mean_cis_fdp_fdr(fdr_bins, stats):

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(fdr_bins[:-1], stats[:, 1][:-1])
    ax.fill_between(fdr_bins[:-1], stats[:, 0][:-1], stats[:, 2][:-1], color='royalblue', alpha=0.2)
    ax.plot([0, 0.1], [0, 0.1], linestyle='--', color='grey')
    
    ax.set_xlabel("FDR")
    ax.set_ylabel("Proxy FDP")
    fig.savefig("./mean_fdp_fdr.png", dpi=800, bbox_inches="tight")

def run_postsearch_validation(config_file_path, target_file, td_file, decoy_file):

    config = open_config(config_file_path)

    logging.info("Started reading the files...")
 
    parser_object = init.ParserInitializer(config).initialize()
    target_df, td_df, decoy_df = FileReading(parser_object).read_search_results(target_file, td_file, decoy_file)
    
    # for now, let's focus on TDC mode, later I will extend this to separated decoy and decoy-free

    score_correction = init.PartitionScoreCorrectionInitializer(config, TEVPartitionCorrection).initialize()
    td_df = score_correction.correct_tdc(td_df)

    logging.info("Finished reading the files.")

    postsearch_partition_object = init.PostSearchPartitionInitializer(config).initialize()
    new_target_df, pep_mod_dict = postsearch_partition_object.add_peptide_modification_index_to_target_df(target_df)
    peptide_id_mapping = postsearch_partition_object.create_peptide_subset_id_mapping(new_target_df)

    logging.info("Peptide subset ID mapping created.")


    logging.info("Starting calculation of proxy FDP vs. FDR contour...")
    fdp_fdr_calculator = BhFdpFdrCalculator
    postsearch_validation = init.PostSearchValidationInitializer(config, fdp_fdr_calculator, new_target_df, td_df, decoy_df, peptide_id_mapping, pep_mod_dict).initialize()
    fdp_fdr_results, _ = postsearch_validation.calculate_fdp_fdr_contour()

    plot_postsearch_results(config, fdp_fdr_results, target_df, target_df, target_df)
    logging.info("Finished the analysis!")
    return fdp_fdr_results, td_df


def plot_postsearch_results(config, fdp_fdr_results, target_df, subject_df, decoy_df):
    """Plotting the contour plot with proxy FDP vs. FDR"""

    plotters = init.PlotPostSearchResultsInitializer(config, fdp_fdr_results, target_df, subject_df, decoy_df).initialize()
    plotters.plot()


def run_sequence_db_partition(config_file_path):

    config = open_config(config_file_path)
    input_file_path = config.get('partition.general', 'search_space_path', fallback='./').strip()

    sequence_dict = FastaParser(input_file_path).parse()

    
    partition_instance = init.SeqDBPartitionInitializer(config, sequence_dict).initialize()
    partition_instance.split_search_space().export_split_search_space()


def run_decoy_generation(input_file_path):

    sequence_dict = FastaParser(input_file_path).parse()
    shuffled_seqs = ShuffledDecoy(sequence_dict, PeptideDecoy, ['K', "R"], ["P"]).generate()
    export_object = PeptideExporter(input_file_path.split('/')[-1].split(".")[0] + "_decoy")
    export_object.export_to_fasta("", shuffled_seqs, mode='DECOY')


def run_searchspace_analysis(input_files):

    all_aa_stats = [get_stats(file) for file in input_files]
    comparison_types = [sa.AAProportionsComparison, sa.MassDistributionComparison, sa.TrypticSitesComparison]
    comparison_runner = sa.AAComparisonRunner(comparison_types)
    results = sa.PartitionQualityEvaluation(all_aa_stats, comparison_runner).run_analysis()

    return results


def get_stats(input_file):

    aa_dict = sp.get_amino_acid_mass_dict()
    protein_dict = FastaParser(input_file).parse()

    aa_data = sa.AminoAcidData(aa_dict, protein_dict)
    aa_counter = sa.AminoAcidCounter(aa_data)
    aa_stats = sa.AminoAcidStatisticsCalculator(aa_counter).get_statistics()

    return aa_stats


def plot_searchspace_analysis_results(config_file_path, results):

    config = open_config(config_file_path)

    aa_proportion_dict, mass_distributions, tryptic_sites = results

    plot_instance = init.PlotSearchSpaceAnalysisInitializer(config).initialize()

    plot_instance.barplot_aa_proportion(aa_proportion_dict)
    plot_instance.plot_mass_distributions(mass_distributions)
    plot_instance.barplot_tryptic_sites(tryptic_sites)