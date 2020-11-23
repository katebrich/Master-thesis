import math
import operator
import os
import random
import shutil
import sys
import time

import numpy as np
from matplotlib import pyplot
from scipy import stats
from collections import Counter
import Plots
from helper import getStructuresFromDirectory
import matplotlib.pyplot as plt

import Logger

logger = Logger.get_logger(os.path.basename(__file__))

#SAMPLE_SIZE = 500
X = random.randint(0,100) #todo for debug
ALPHA = 0.05

class AnalysisComputer():
    output_dir = ""
    lbs_dir = ""
    features_dir = ""
    config = ""
    p_values = []
    b_ratios = []
    p_val_perc = []
    means = []
    errors = []
    features_list = []

    def __init__(self, output_dir, lbs_dir, features_dir, features_list, config):
        self.output_dir = output_dir
        self.lbs_dir = lbs_dir
        self.config = config
        self.features_list = features_list
        self.features_dir = features_dir

    def run(self, sample_size, iterations, balance_binding_ratio, draw_plots):
        start = time.time()
        logger.info(f"Running analysis started...")

        for feature in self.features_list:
            feature_output_dir = os.path.join(self.output_dir, feature)
            if os.path.exists(feature_output_dir):
                shutil.rmtree(feature_output_dir)
            os.makedirs(feature_output_dir)
            self.process_feature(feature, sample_size, iterations, balance_binding_ratio, draw_plots)
            logger.info(f"Feature {feature}: done")

        logger.info(f"Running analysis finished in {math.ceil(time.time() - start)}s.")

        ''' #todo for debug
        #plot the results
        barWidth = 0.45
        bars = [x[1] for x in self.p_values]
        r1 = np.arange(len(bars))
        r2 = [x + barWidth/2 for x in r1]
        pyplot.bar(r2, bars, color='g', width=barWidth, edgecolor='white')'''

    def process_feature(self, feature, sample_size, iterations, balance_binding_ratio, draw_plots):
        feature_dir = os.path.join(self.features_dir, feature)
        dataset = getStructuresFromDirectory(
            feature_dir)  # compute only for structures that have the feature computed, ignore the rest
        feature_output_dir = os.path.join(self.output_dir, feature)
        feature_type = self.config.get_feature_type(feature) # binary/continuous/categorical/ordinal
        pairs = []
        statistics = []
        p_values = []
        b_ratios = []

        #pair feature values with ligand binding sites
        for structure in dataset:
            pairs = pairs + self.compute_pairs(structure, self.lbs_dir, feature_dir, feature)

        #save pairs
        file = os.path.join(feature_output_dir, f"pairs.txt")
        with open(file, 'w') as f:
            f.write('\n'.join('{} {}'.format(x[0], x[1]) for x in pairs))

        data_binding = [x[1] for x in pairs if x[0] == 1]
        data_nonbinding = [x[1] for x in pairs if x[0] == 0]

        #compute hypothesis test for all iterations
        #random.seed(42 + X)  # todo for debug
        with open(os.path.join(feature_output_dir, f"iterations.txt"), 'w') as f:
            for i in range(1, iterations + 1):
                f.write(f"********************************\n")
                f.write(f"**        ITERATION {i}       **\n")
                f.write(f"********************************\n")

                if (sample_size == 0): #no sampling
                    sample_binding = data_binding
                    sample_nonbinding = data_nonbinding
                elif not balance_binding_ratio: #do not balance the ratio, sample from the whole dataset
                    sample = random.sample(pairs, sample_size)
                    sample_binding = [x[1] for x in sample if x[0] == 1]
                    sample_nonbinding = [x[1] for x in sample if x[0] == 0]
                else: #simple sampling without replacement
                    sample_binding = random.sample(data_binding, sample_size)
                    sample_nonbinding = random.sample(data_nonbinding, sample_size)

                #run analysis and save results to file
                if (feature_type == "continuous"):
                    res = self.welchs_t_test(sample_binding, sample_nonbinding, f)
                elif (feature_type == "categorical" or feature_type == "ordinal" or feature_type == "binary"):
                    res = self.chi_squared_test(sample_binding, sample_nonbinding, f, feature)
                else:
                    logger.error(f"Unknown type of feature '{feature_type}'. Please specify the type in config.")
                    self.errors.append(feature)
                    return
                if res == None:
                    self.errors.append(feature)
                    return

                statistics.append(res[0])
                p_values.append(res[1])
                b_ratios.append(len(sample_binding)/(len(sample_binding) + len(sample_nonbinding)))

        if draw_plots:
            self.draw_plots(feature_type, data_binding, data_nonbinding, feature_output_dir, feature)

        #save values for summary
        self.p_values.append((feature, p_values))
        self.b_ratios.append((feature, b_ratios))

        #save difference of means and variance for summary
        if (feature_type == "continuous"):
            mean_binding = np.mean(data_binding)
            mean_nonbinding = np.mean(data_nonbinding)
            means_dif =  abs(mean_binding - mean_nonbinding)
            variance_binding = np.var(data_binding)
            variance_nonbinding = np.var(data_nonbinding)
            self.means.append((feature, (means_dif, mean_binding, mean_nonbinding, variance_binding, variance_nonbinding)))

        #calculate percentage of iterations for which p-value was under significance level alpha
        under_alpha = 0
        for p in p_values:
            if p <= ALPHA:
                under_alpha = under_alpha + 1
        perc = (under_alpha / iterations) * 100
        self.p_val_perc.append((feature, round(perc,2)))

        #plot statistics
        plt.clf()
        plt.scatter(range(1,iterations + 1), statistics)
        plt.xlabel('iteration')
        if (feature_type == "continuous"):
            plt.ylabel('T-statistic')
        else:
            plt.ylabel('Chi-squared')
        plt.savefig(os.path.join(feature_output_dir,f"{feature}_statistic.png"))
        #plot p-values
        plt.clf()
        plt.scatter(range(1, iterations + 1), p_values)
        axes = plt.gca()
        axes.set_ylim([0, 1])
        plt.xlabel('iteration')
        plt.ylabel('P-value')
        plt.savefig(os.path.join(feature_output_dir, f"{feature}_p_values.png"))

        #save p-values to file
        with open(os.path.join(feature_output_dir, f"p_values.txt"), 'w') as f:
            f.write('\n'.join(f"{x}" for x in p_values))

    #todo zkontrolovat logovani
    def compute_pairs(self, structure, lbs_dir, feature_dir, feature_name):
        pdb_id = structure[0]
        chain_id = structure[1]
        error = False
        pairs = []

        try:
            # get ligand binding sites values
            file = os.path.join(lbs_dir, f"{pdb_id}{chain_id}.txt")
            lbs = np.genfromtxt(file, delimiter=' ', dtype=None)
            lbs_dict = dict(lbs)

            # get feature values
            file = os.path.join(feature_dir, f"{pdb_id}{chain_id}.txt")
            feature = np.genfromtxt(file, delimiter=' ', dtype=None, encoding=None)
            feature_vals = list(feature)

            for val in feature_vals:
                res_num = val[0]
                feature_val = val[1]
                if res_num in lbs_dict:
                    lbs_val = lbs_dict[res_num]
                else:
                    continue
                pairs.append((lbs_val, feature_val))

            return pairs
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as ex:
            error = True
            logger.debug(f"{feature_name}: Error while processing {pdb_id} {chain_id}: {ex}")
        finally:
            if (error):
                logger.error(f"{feature_name}: {pdb_id} {chain_id} NOT PROCESSED ! See log for more details.")

    def draw_plots(self, feature_type, data_binding, data_nonbinding, feature_output_dir, feature_name):
        if (feature_type == "binary" or feature_type == "categorical"):
            Plots.plot_binding_nonbinding_ratios(data_binding, data_nonbinding,
                                                 os.path.join(feature_output_dir, f"{feature_name}_ratios"))
        elif feature_type == "continuous":
            Plots.plot_histogram(data_binding, data_nonbinding, 40,
                               os.path.join(feature_output_dir, f"{feature_name}_hist_bins_40"))
        elif feature_type == "ordinal":
            Plots.plot_binding_nonbinding_ratios(data_binding, data_nonbinding,
                                                 os.path.join(feature_output_dir, f"{feature_name}_ratios"))
            Plots.plot_histogram(data_binding, data_nonbinding, 40,
                                 os.path.join(feature_output_dir, f"{feature_name}_hist_bins_40"))
            Plots.plot_histogram(data_binding, data_nonbinding, 100, os.path.join(feature_output_dir, f"{feature_name}_hist_bins_100"))

    def write_summary(self):
        self.p_values.sort(key=lambda x: x[0])  # sort by name
        with open(os.path.join(self.output_dir, f"p_values.csv"), 'w') as f:
            for line in self.p_values:
                f.write(f"{str(line[0])}, ")
                f.write(', '.join('{0:.4g}'.format(x) for x in line[1]))
                f.write('\n')

        self.b_ratios.sort(key=lambda x: x[0])  # sort by name
        with open(os.path.join(self.output_dir, f"binding_ratios.csv"), 'w') as f:
            #f.write(f"Feature, # Binding, # Nonbinding, Binding/Total\n") # header
            for line in self.b_ratios:
                f.write(f"{str(line[0])}, ")
                f.write(', '.join('{0:.5f}'.format(x) for x in line[1]))
                f.write('\n')

        self.means.sort(key=lambda x: x[1][0])  # sort by means difference
        with open(os.path.join(self.output_dir, f"means_difference.csv"), 'w') as f:
            f.write(f"Feature, Means difference, Mean binding, Mean nonbinding, Variance binding, Variance nonbinding\n") # header
            for line in self.means:
                f.write(f"{str(line[0])}, ")
                f.write(', '.join('{0:.4f}'.format(x) for x in line[1]))
                f.write('\n')

        self.p_val_perc.sort(key=lambda x: x[1], reverse = True)  # sort by percentage
        with open(os.path.join(self.output_dir, f"p_vals_perc.csv"), 'w') as f:
            for line in self.p_val_perc:
                f.write(f"{str(line[0])}, {str(line[1])}%")
                f.write('\n')
        with open(os.path.join(self.output_dir, f"errors.txt"), 'w') as f:
            for line in self.errors:
                f.write(f"{line}\n")

    def welchs_t_test(self, sample_binding, sample_nonbinding, f):
        out = stats.ttest_ind(sample_binding, sample_nonbinding, axis=0, equal_var=False, nan_policy='raise') #equal_var is False to run Welch's test instead of T-test
        t_statistic = out[0]
        p_value = float(out[1])

        f.write(f"T-statistic: {t_statistic}\n")
        f.write(f"P-value: {p_value}\n")
        f.write(f"Sample binding -> mean: {np.mean(sample_binding)}; variance: {np.var(sample_binding)}\n")
        f.write(f"Sample non-binding -> mean: {np.mean(sample_nonbinding)}; variance: {np.var(sample_nonbinding)}\n")
        f.write(f"Difference of means: {np.mean(sample_binding) - np.mean(sample_nonbinding)}\n")

        return (t_statistic, p_value)

    def chi_squared_test(self, sample_binding, sample_nonbinding, f, feature):
        sample_binding = [str(i) for i in sample_binding] #fix for Counter
        sample_nonbinding = [str(i) for i in sample_nonbinding]
        binding_counts = Counter(sample_binding)
        nonbinding_counts = Counter(sample_nonbinding)

        keys1 = binding_counts.keys()
        keys2 = nonbinding_counts.keys()
        categories = list(set().union(keys1, keys2))
        categories = sorted(categories)

        binding = []
        nonbinding = []
        binding_counts_sorted = []
        nonbinding_counts_sorted = []
        for cat in categories:
            binding.append(binding_counts[cat])
            nonbinding.append(nonbinding_counts[cat])
            binding_counts_sorted.append((cat, binding_counts[cat]))
            nonbinding_counts_sorted.append((cat, nonbinding_counts[cat]))

        obs = np.array([np.array(binding), np.array(nonbinding)])
        chi2, p_value, dof, expected = stats.chi2_contingency(obs)

        f.write(f"Binding counts:\n")
        f.write(f"{binding_counts_sorted}\n")
        f.write(f"Non-binding counts:\n")
        f.write(f"{nonbinding_counts_sorted}\n")
        f.write(f"Expected values: \n")
        f.write(f"{expected}\n")

        error = False
        invalid_count = 0
        #check assumptions
        for cell in np.nditer(expected):
            if (cell <= 1): #no cell can have expected value less than 1
                error = True
            if (cell <=5):
                invalid_count = invalid_count + 1
        if (invalid_count >= 0.2 * len(categories) * 2):
            error = True
        if error:
            logger.error(f"Feature {feature}: Chi-squared test assumptions not met. More than 20% of cells have expected value of 5 or less, and/or at least one cell has expected value of 1 or less.")
            return None

        p_value = float(p_value)

        f.write(f"Chi squared: {chi2}\n")
        f.write(f"p-value: {p_value}\n")

        return (chi2, p_value)

    '''def fischers_exact_test(self, data, results_file, feature):
            with open(results_file, 'w') as f:
                counts = Counter(data)

                # todo opravit..kdyz neni ani jedna 1, pada
                binding_positive = counts[(1, 1)]
                binding_negative = counts[(1, 0)]
                non_binding_positive = counts[(0, 1)]
                non_binding_negative = counts[(0, 0)]

                odds_ratio, p_value = stats.fisher_exact([[binding_positive, binding_negative],
                                                          [non_binding_positive, non_binding_negative]])
                p_value = float(p_value)

                f.write(f"Binding positive count: {binding_positive}\n")
                f.write(f"Binding negative count: {binding_negative}\n")
                f.write(f"Nonbinding positive count: {non_binding_positive}\n")
                f.write(f"Nonbinding negative count: {non_binding_negative}\n")

                f.write(f"odds ratio: {odds_ratio}\n")
                f.write(f"p-value: {p_value}\n")
                f.write(f"Binding sites -> positives ratio: {binding_positive / (binding_positive + binding_negative)}\n")
                f.write(
                    f"Non-binding sites -> positives ratio: {non_binding_positive / (non_binding_positive + non_binding_negative)}\n")

            return (feature, p_value, binding_positive + binding_negative, non_binding_positive + non_binding_negative, "Fisher's exact")'''
