from abc import ABC, abstractmethod
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


class DerivativeMode(ABC):

    @abstractmethod
    def calculate_derivative(self, quality_thresholds, fdp_vals):
        pass


class ForwardDifference(DerivativeMode):

    def calculate_derivative(self, quality_thresholds, fdp_vals):

        derivatives = []

        for idx in range(len(quality_thresholds)-1):
            fdp_idx = fdp_vals[idx]
            fdp_idx_plus_one = fdp_vals[idx+1]

            thr_idx = quality_thresholds[idx]
            thr_idx_plus_one = quality_thresholds[idx+1]

            cur_derivatives = abs((thr_idx_plus_one - thr_idx) / (fdp_idx_plus_one - fdp_idx))
            #cur_derivatives[np.isnan(cur_derivatives)] = -10
            derivatives.append(cur_derivatives)

        return derivatives

class CentralDifference(DerivativeMode):

    def calculate_derivative(self, quality_thresholds, fdp_vals):
        
        derivatives = []

        for idx in range(1, len(quality_thresholds)-1):
            fdp_idx_minus_one = fdp_vals[idx-1]
            # fdp_idx = fdp_vals[idx]
            fdp_idx_plus_one = fdp_vals[idx+1]

            thr_idx_minus_one = quality_thresholds[idx-1]
            # thr_idx = quality_thresholds[idx]
            thr_idx_plus_one = quality_thresholds[idx+1]

            cur_derivatives = abs((thr_idx_plus_one - thr_idx_minus_one) / (fdp_idx_plus_one - fdp_idx_minus_one))
            #cur_derivatives[np.isnan(cur_derivatives)] = -10
            derivatives.append(cur_derivatives)

        return derivatives


class DerivativeCalculator:

    @staticmethod
    def calculate_partial_derivatives(quality_thresholds, fdp_digitized_means, derivative_mode: DerivativeMode):
        return derivative_mode.calculate_derivative(quality_thresholds, fdp_digitized_means)
    

class MovingAverageCalculator:

    @staticmethod
    def calculate_moving_average(data, win_size=5):
    
        numbers_series = pd.Series(data)
        windows = numbers_series.rolling(win_size)
        moving_averages = windows.mean()
        return np.array(moving_averages.tolist()[win_size - 1:])


class ThresholdFinder:

    def run(self, norm_der_means, thresholds):
        min_th, max_th = self.find_optimal_threshold_range(norm_der_means, thresholds)
        return min_th, max_th

    @staticmethod
    def find_optimal_threshold_range(norm_der_means, thresholds):

        window = int(0.1 * len(thresholds))
        ave_derivatives = MovingAverageCalculator.calculate_moving_average(norm_der_means, window)
        ave_thresholds = MovingAverageCalculator.calculate_moving_average(thresholds, window)
        return ave_thresholds[ave_derivatives.argmin()]


    # @staticmethod
    # def get_optimal_threshold_idxs(thresholds, opt_min, opt_max):
    #     min_idx = list(thresholds < opt_min).index(False) + 1
    #     max_idx = list(thresholds < opt_max).index(False) + 1
    #     return min_idx, max_idx




class OptimalQualityThresholdFinder:

    def __init__(self, fdp_fdr_results):
        self.fdr_bins = None
        self.fdr_bins_idxs = None
        self.fdp_fdr_results = fdp_fdr_results

    
    def run(self, thresholds):
        fdr_digitized = self.discretize_fdr_bins()
        fdp_bin_means = self.take_all_bins_fdp_means(fdr_digitized)
        derivatives = DerivativeCalculator.calculate_partial_derivatives(thresholds, fdp_bin_means, ForwardDifference())
        der_means = self.take_derivative_means(derivatives)
        norm_der_means = der_means / max(der_means)
        opt_range = ThresholdFinder.find_optimal_threshold_range(norm_der_means, thresholds)
        # opt_idxs = ThresholdFinder.get_optimal_threshold_idxs(thresholds, *opt_range)

        return norm_der_means, opt_range


    def find_min_max_fdrs(self):
        min_val = min(min(arr[0]) for arr in self.fdp_fdr_results)
        max_val = max(max(arr[0]) for arr in self.fdp_fdr_results)
        return min_val, max_val


    def discretize_fdr_bins(self, num_bins=50):
        
        range_vals = self.find_min_max_fdrs()
        self.fdr_bins = np.linspace(*range_vals, num_bins)
        self.fdr_bins_idxs = range(1, len(self.fdr_bins)+1)
        fdr_digitized = [np.digitize(arr[0], self.fdr_bins) for arr in self.fdp_fdr_results]

        return fdr_digitized


    def take_all_bins_fdp_means(self, fdr_digitized):
        
        fdp_arrays = [arr[1] for arr in self.fdp_fdr_results]
        zipped = zip(fdr_digitized, fdp_arrays)
        fdp_digitized_means = [self.take_single_bin_means(fdr_arr, fdp_arr) for fdr_arr, fdp_arr in zipped]

        return np.array(fdp_digitized_means)


    def take_single_bin_means(self, fdr_arr, fdp_arr):
        return [np.mean(fdp_arr[fdr_arr == fdr_bin_idx]) for fdr_bin_idx in self.fdr_bins_idxs]


    def take_derivative_means(self, derivatives):
        return np.nanmean(np.array(derivatives), axis=1)


    @staticmethod
    def calculate_moving_average(norm_der_means, win_size=5):

        numbers_series = pd.Series(norm_der_means)
        windows = numbers_series.rolling(win_size)
        moving_averages = windows.mean()
        moving_averages_list = moving_averages.tolist()

        return np.array(moving_averages_list[win_size - 1:])






class Visualization:

    def run_plots(self, norm_der_means, thresholds, window_size, opt_threshold):

        mov_derivatives= MovingAverageCalculator.calculate_moving_average(norm_der_means, window_size)
        mov_thresholds = MovingAverageCalculator.calculate_moving_average(thresholds, window_size)
        self.plot_moving_average(mov_derivatives, mov_thresholds, opt_threshold)

    @staticmethod
    def plot_moving_average(mov_derivatives, mov_thresholds, opt_threshold):
        
        max_der = max(mov_derivatives)
        print(mov_derivatives)
        fig, ax = plt.subplots(figsize=(7,5))

        ax.plot(mov_thresholds[:-1], mov_derivatives, color='royalblue')
        ax.set_ylabel("derivative (moving average)")
        ax.set_xlabel("quality score threshold")
        #ax.vlines(x=opt_thresholds, ymin=2 * (min(mov_derivatives)*0.9,), ymax=2 * (max(mov_derivatives)*1.1,))
        ax.vlines(x=opt_threshold, ymin = 0, ymax=max_der*1.1)
        ax.set_ylim(0, max_der*1.1)
        fig.savefig("./moving_average_threshold.png", dpi=800, bbox_inches="tight")