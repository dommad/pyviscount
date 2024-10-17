"""Execution scripts to run different functionalities of PyViscount"""
import logging
import numpy as np
import matplotlib.pyplot as plt
from . import sequence_processing as sp
from . import searchspace_analysis as sa
from . import initializer as init
from . import optimal_selection as osel
from .utils import open_config
from .template import FileReading, AveragedFdpFdr, ConfidenceInterval
from .decoy_generation import ShuffledDecoy, PeptideDecoy
from .exporter import PeptideExporter
from .postsearch import PostSearchOrchestrated, ThresholdEvaluator
from .scores import ScoreProcessing
from .parsers import FastaParser



def read_files(config, *args):

    parser_object = init.ParserInitializer(config).initialize()
    search_results = FileReading(parser_object).read_search_results(*args)

    return search_results



def run_presearch_validation(config_file_path, full_target_file, partition_target_file, partition_td_file):

    config = open_config(config_file_path)

    logging.info("Started reading the files...")

    full_df, part_t_df, part_td_df = read_files(config_file_path, full_target_file, partition_target_file, partition_td_file)
    #full_df = add_percolator_score(full_df)

    logging.info("Starting calculation of proxy FDP vs. FDR contour...")
    validation_object = init.PreSearchValidationInitializer(config, full_df, part_t_df, part_td_df).initialize()
    fdp_fdr_results, thresholds = validation_object.calculate_fdp_fdr_contour()

    plot_postsearch_results(config, fdp_fdr_results)
    logging.info("Finished the analysis!")
    return fdp_fdr_results, full_df



def run_postsearch_validation(config_file, target_file, subject_file=None, decoy_file=None):

    config = open_config(config_file)

    subject_file = subject_file or target_file

    logging.info("Started reading the files...")
    target_df, subject_df, decoy_df = read_files(config, target_file, subject_file, decoy_file)
    logging.info("Finished reading the files.")

    target_df, subject_df, decoy_df = process_scores(config, target_df, subject_df, decoy_df)

    logging.info("Starting calculation of proxy FDP vs. FDR contour...")
    fdp_fdr_results = calculate_fdp_fdr(config, target_df, subject_df, decoy_df)

    plot_postsearch_results(config, fdp_fdr_results)

    opt_threshold = get_optimal_threshold(config, fdp_fdr_results)
    bootstrap_results = bootstrap_fdp_fdr(config, opt_threshold, target_df, subject_df, decoy_df)
    
    plot_optimal_fdp_fdr(bootstrap_results)
    logging.info("Finished the analysis!")

    return fdp_fdr_results, bootstrap_results


def calculate_fdp_fdr(config, target_df, subject_df, decoy_df):
    """Run post-search validation and return FDP vs. FDR results."""
    return PostSearchOrchestrated(config).run_postsearch_validation(target_df, subject_df, decoy_df)


def bootstrap_fdp_fdr(config, opt_threshold, target_df, subject_df, decoy_df):
    """Run bootstrap analysis for optimal threshold."""
    return PostSearchOrchestrated(config).run_postsearch_bootstrap(opt_threshold, target_df, subject_df, decoy_df)


def process_scores(config, target_df, subject_df, decoy_df):
    """Correct and process the scores based on the configuration."""
    score_processing = ScoreProcessing(config)
    subject_df = score_processing.correct_scores(subject_df)
    subject_df = score_processing.add_postprocessor_scores(subject_df)

    # apply Sidak correction
    if score_processing.sidak_status == 'true':
        target_df = score_processing.correct_sidak(target_df)
        subject_df = score_processing.correct_sidak(subject_df)

    # if user wants to use decoy-based p-values from separate target-decoy search
    if score_processing.validation_mode == 'Std' and score_processing.decoy_mode == 'P-value':
        target_df = score_processing.add_std_decoy_pval(target_df, decoy_df)
        subject_df = score_processing.add_std_decoy_pval(subject_df, decoy_df)

    logging.info("Corrected search space-dependent scores and added post-processor scores (if provided).")
    return target_df, subject_df, decoy_df



def get_optimal_threshold(config, fdp_fdr_results):

    thresholds = ThresholdEvaluator(config).threshold_dictionary()

    oq = osel.OptimalQualityThresholdFinder(fdp_fdr_results)
    norm_der_means, opt_threshold = oq.run(thresholds)

    # maybe it'd be better to abstract the number of thresholds out by taking 10% of the number of original thresholds
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



def plot_postsearch_results(config, fdp_fdr_results):
    """Plotting the contour plot with proxy FDP vs. FDR"""

    plotters = init.PlotPostSearchResultsInitializer(config, fdp_fdr_results).initialize()
    plotters.plot()



# analysis of the partition process

def run_sequence_db_partition(config_file_path):
    """
    Partition the sequence database and export the split search space.
    """

    try:

        config = open_config(config_file_path)
        input_file_path = config.get('partition.general', 'search_space_path', fallback='./').strip()

        sequence_dict = FastaParser(input_file_path).parse()

        
        partition_instance = init.SeqDBPartitionInitializer(config, sequence_dict).initialize()
        partition_instance.split_search_space().export_split_search_space()
    
    except Exception as e:
        logging.error(f"Error during sequence DB partition: {e}")
        raise


def run_decoy_generation(input_file_path):
    """
    Generate decoy sequences from input file and export to FASTA.
    """
    try:
        sequence_dict = FastaParser(input_file_path).parse()
        shuffled_seqs = ShuffledDecoy(sequence_dict, PeptideDecoy, ['K', "R"], ["P"]).generate()
        export_object = PeptideExporter(input_file_path.split('/')[-1].split(".")[0] + "_decoy")
        export_object.export_to_fasta("", shuffled_seqs, mode='DECOY')
    
    except Exception as e:
        logging.error(f"Error during decoy generation: {e}")
        raise


def run_searchspace_analysis(input_files):
    """
    Analyze search space based on input files.
    """
    try:
        all_aa_stats = [get_stats(file) for file in input_files]
        comparison_types = [sa.AAProportionsComparison, sa.MassDistributionComparison, sa.TrypticSitesComparison]
        comparison_runner = sa.AAComparisonRunner(comparison_types)
        results = sa.PartitionQualityEvaluation(all_aa_stats, comparison_runner).run_analysis()

        return results

    except Exception as e:
        logging.error(f"Error during search space analysis: {e}")
        raise


def get_stats(input_file):
    """
    Retrieve statistics from a given input file.
    """

    try:
        aa_dict = sp.get_amino_acid_mass_dict()
        protein_dict = FastaParser(input_file).parse()

        aa_data = sa.AminoAcidData(aa_dict, protein_dict)
        aa_counter = sa.AminoAcidCounter(aa_data)
        aa_stats = sa.AminoAcidStatisticsCalculator(aa_counter).get_statistics()

        return aa_stats
    
    except Exception as e:
        logging.error(f"Error while getting stats for file {input_file}: {e}")
        raise


def plot_searchspace_analysis_results(config_file_path, results):
    """
    Plot results from the search space analysis.
    """
    try:
        config = open_config(config_file_path)

        aa_proportion_dict, mass_distributions, tryptic_sites = results

        plot_instance = init.PlotSearchSpaceAnalysisInitializer(config).initialize()

        plot_instance.barplot_aa_proportion(aa_proportion_dict)
        plot_instance.plot_mass_distributions(mass_distributions)
        plot_instance.barplot_tryptic_sites(tryptic_sites)
    
    except Exception as e:
        logging.error(f"Error while plotting search space analysis results: {e}")
        raise

