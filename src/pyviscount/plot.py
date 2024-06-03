from typing import List, Dict
from functools import partial
import numpy as np
from matplotlib import pyplot as plt, colormaps as cm
import seaborn as sns
#from .sequence_processing import get_aa_proportion_difference
from abc import ABC, abstractmethod



def get_rgba_colors(*args):
    """Generate matplotlib-compatible form of RGBA colors"""
    cols = []
    for color in args:
        new = [x / 255 for x in color[:-1]]
        new.append(color[-1])
        cols.append(new)
    return iter(cols)


class Plotter(ABC):

    @abstractmethod
    def plot(self):
        pass


class FdpFdrPlotter(Plotter):

    def __init__(self, fdp_fdr_results, output_name):
        self.fdp_fdr_results = fdp_fdr_results
        self.output_name = output_name

    def plot(self, ylim=0.1):
        """
        Plot proxy False Discovery Proportion (FDP) vs. False Discovery Rate (FDR).

        Parameters:
        - decs_negs (list): A list containing tuples of FDR and FDP arrays.
        - out_name (str): The name to be used for the output plot file.

        Returns:
        None
        """

        plt.style.use('ggplot')
        plt.rcParams.update({'font.size': 13, 'font.family': 'Helvetica'})

        n_rep = len(self.fdp_fdr_results)
        my_map = cm.get_cmap("Blues")
        my_colors = iter(my_map(np.linspace(0, 1, n_rep)))
        fig, ax = plt.subplots(figsize=(7,5), constrained_layout=True)

        for fdrs, fdps in self.fdp_fdr_results:
            ax.plot(fdrs, fdps, color=next(my_colors))

        ax.plot([0, 1], [0, 1], color='grey')
        ax.set_xlim(0, 0.1)
        ax.set_ylim(0, ylim)
        ax.set_xlabel("FDR")
        ax.set_ylabel("proxy FDP")

        fig.savefig(f"./{self.output_name}_fdp_fdr_contour.png", dpi=500, bbox_inches="tight")



class HistogramPlotter(Plotter):
    
    def __init__(self, org_df, td_df, decoy_df, plot_type='kde', score='tev'):
        self.org_df = org_df
        self.td_df = td_df
        self.decoy_df = decoy_df
        self.plot_type = plot_type
        self.score = score

    def plot(self):
        """Plot KDE or histogram for the results of target, target-decoy, decoy search"""

        fig, ax = plt.subplots(figsize=(5,5))
        plot_to_function = {
            'kde': partial(sns.kdeplot, fill=True),
            'hist': partial(sns.histplot, kde=False, fill=True, element="poly", stat='density'),
        }


        plot_to_function[self.plot_type](
            self.org_df[self.org_df['hit_rank'] != 1][self.score],
            label='incorrect target')

        plot_to_function[self.plot_type](
            self.td_df[self.td_df['is_decoy']][self.score],
            label='decoy (from TD)')

        plot_to_function[self.plot_type](
            self.decoy_df[self.score],
            label='decoy (from D)')

        plt.legend()
        plt.xlabel(self.score)
        fig.savefig('./hist_test.png', dpi=500)



class PlotValidation:
    """Plotting results of the validation by post-search partition"""

    def __init__(self, *args):
        self.plotters = args

    def plot(self):
        for plotter in self.plotters:
            plotter.plot()



class PlotSearchSpaceAnalysis:

    def __init__(self, config, aa_dict: dict):
        
        self.aa_dict = aa_dict
        self.output_path = config.get('partition.general', 'output_path', fallback='./').split(" ")[0].strip()
        self.file_format = config.get('plotting', 'plot_format', fallback='pdf').split(" ")[0].strip()
        self.dpi = int(config.get('plotting', 'dpi', fallback=300).split(" ")[0].strip())

        plt.style.use('ggplot')
        plt.rcParams.update({'font.size': 12, 'font.family': 'Helvetica'})


    def save_figure(self, fig, core_name):
        """general function for saving the figures"""
        fig.savefig(f"{self.output_path}{core_name}.{self.file_format}", dpi=self.dpi)



    def barplot_aa_proportion(self, aa_proportion_dict, difference=False, difference_dataset=None):
        
        # Set the width of the bars
        bar_width = 0.35
        # Set the positions of the bars on the x-axis
        bar_positions = np.arange(len(self.aa_dict))

        dataset_a, dataset_b = aa_proportion_dict.values()

        fig, axes = plt.subplots(figsize=(7, 5), constrained_layout=True)

        if difference:
            axes.bar(bar_positions, difference_dataset)
        else:
            axes.bar(bar_positions - bar_width/2, dataset_a, width=bar_width, label='subset A')
            axes.bar(bar_positions + bar_width/2, dataset_b, width=bar_width, label='subset B')

        # Add labels and title
        axes.set_xlabel('Amino Acid')
        axes.set_ylabel('Fraction in the subset')
        #plt.title('Bar Plot of Two Datasets')
        axes.set_xticks(bar_positions, self.aa_dict)
        axes.legend()

        self.save_figure(fig, "aa_proportion")



    def plot_mass_distributions(self, datasets):

        fig, axes = plt.subplots(figsize=(7, 5), constrained_layout=True)
        sns.kdeplot(data=datasets, fill=True, multiple='layer', ax=axes)
        axes.legend(['subset A', 'subset B'])
        axes.set_xlabel("Mass [Da]")

        self.save_figure(fig, "mass_distribution")


    
    def barplot_tryptic_sites(self, tryptic_sites_dict_list: List[Dict]):

        dataset_a_dict, dataset_b_dict = tryptic_sites_dict_list
        categories = dataset_a_dict.keys()

        bar_width = 0.2
        bar_positions = np.arange(len(categories))

        fig, axes = plt.subplots(figsize=(7, 5), constrained_layout=True)

        axes.bar(bar_positions - bar_width / 2, dataset_a_dict.values(), width=bar_width, label='subset A')
        axes.bar(bar_positions + bar_width / 2, dataset_b_dict.values(), width=bar_width, label='subset B')

        # Add labels and title
        axes.set_xlabel('Tryptic site')
        axes.set_ylabel('Total number in the subset')
        #plt.title('Bar Plot of Two Datasets')
        axes.set_xticks(bar_positions, categories)
        axes.legend()

        self.save_figure(fig, "tryptic_sites")



