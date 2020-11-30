import math
import os
import shutil
import sys
import time

import numpy as np
from scipy import stats
from collections import Counter
from AnalysisPipeline import Plots
from helper import getStructuresFromDirectory
#from unused.statistical_analysis import welchs_t_test, fischers_exact_test, chi_squared_test

import Logger

logger = Logger.get_logger(os.path.basename(__file__))

class AnalysisComputer():
    output_dir = ""
    lbs_dir = ""
    config = ""
    feature_dir = ""
    feature_name = ""
    errors = []
    pairs = []
    total = ""
    counter = 1
    summary = []

    def __init__(self, output_dir, lbs_dir, config):
        self.output_dir = output_dir
        self.lbs_dir = lbs_dir
        self.config = config

    def run(self, feature_name, feature_dir):
        self.feature_dir = feature_dir
        self.feature_name = feature_name

        feature_output_dir = os.path.join(self.output_dir, feature_name)

        if os.path.exists(feature_output_dir):
            shutil.rmtree(feature_output_dir)
        os.makedirs(feature_output_dir)

        dataset = getStructuresFromDirectory(self.feature_dir) #compute only for structures that have the feature computed, ignore the rest

        start = time.time()
        logger.info(f"Running analysis for feature '{self.feature_name}' started...")

        self.pairs = []
        self.errors = []
        self.total = len(dataset)
        self.counter = 1

        #pair feature values with ligand binding sites
        for structure in dataset:
            self.compute_pairs(structure)

        feature_type = self.config.get_feature_type(self.feature_name)

        #save pairs
        file = os.path.join(feature_output_dir, f"pairs.txt")
        with open(file, 'w') as f:
            f.write('\n'.join('{} {}'.format(x[0], x[1]) for x in self.pairs))

        binding_data = [x[1] for x in self.pairs if x[0] == 1]
        nonBinding_data = [x[1] for x in self.pairs if x[0] == 0]

        #run analysis and save results
        file = os.path.join(feature_output_dir, f"results.txt")
        #if (feature_type == "binary"):
        #    #logger.info("Running Fischer's exact test")
        #    self.summary.append(self.fischers_exact_test(self.pairs, file, feature_name))
        #    Plots.plot_binding_nonbinding_ratios(self.pairs, feature_output_dir)
        if (feature_type == "continuous"):
            #logger.info("Running Welch's T-test")
            self.summary.append(self.welchs_t_test(self.pairs, file, feature_name))
            Plots.plot_histogram(binding_data, nonBinding_data, 50, os.path.join(feature_output_dir, f"{feature_name}_bins_40"))
            Plots.plot_histogram(binding_data, nonBinding_data, 75, os.path.join(feature_output_dir, f"{feature_name}_bins_75"))
            Plots.plot_histogram(binding_data, nonBinding_data, 100, os.path.join(feature_output_dir, f"{feature_name}_bins_100"))
        elif (feature_type == "categorical" or feature_type == "ordinal" or feature_type == "binary"):
            #logger.info("Running Chi-squared test")
            self.summary.append(self.chi_squared_test(self.pairs, file, feature_name))
            Plots.plot_binding_nonbinding_ratios(binding_data, nonBinding_data, os.path.join(feature_output_dir, f"{feature_name}_ratios"))
        else:
            #todo
            logger.error(f"Unknown type of feature '{feature_type}'. Please specify the type in config.")
            sys.exit(1)

        if (len(self.errors) == 0):
            logger.info(f"Running analysis finished in {math.ceil(time.time() - start)}s. Results saved to {file}. All structures processed successfully.")
        else:
            errors_format = '\n'.join('%s %s' % x for x in self.errors)
            logger.warning(f"Running analysis finished in {math.ceil(time.time() - start)}s. Results saved to {file}. {len(self.errors)}/{self.total} structures were not processed successfully - they were skipped in the analysis! \n{errors_format}")

    def compute_pairs(self, structure):
        pdb_id = structure[0]
        chain_id = structure[1]
        error = False
        missing_vals = []

        try:
            # get ligand binding sites values
            file = os.path.join(self.lbs_dir, f"{pdb_id}{chain_id}.txt")
            lbs = np.genfromtxt(file, delimiter=' ', dtype=None)
            lbs_dict = dict(lbs)

            # get feature values
            file = os.path.join(self.feature_dir, f"{pdb_id}{chain_id}.txt")
            feature = np.genfromtxt(file, delimiter=' ', dtype=None, encoding=None)
            feature_vals = list(feature)

            for val in feature_vals:
                res_num = val[0]
                feature_val = val[1]
                if res_num in lbs_dict:
                    lbs_val = lbs_dict[res_num]
                else:
                #    missing_vals.append(res_num)
                    continue
                self.pairs.append((lbs_val, feature_val)) #todo aby to ta metoda vracela misto rovnou prirazovala
            #if (len(missing_vals) > 0):
                #logger.debug(f"{pdb_id} {chain_id}: Missing feature values for residues: {missing_vals}")
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as ex:
            error = True
            logger.debug(f"Error while processing {pdb_id} {chain_id}: {ex}")
        finally:
            self.counter  += 1
            if (error):
                self.errors.append((structure[0], structure[1]))
                logger.error(f"{self.counter}/{self.total}: {pdb_id} {chain_id} NOT PROCESSED ! See log for more details.")
            #else:
            #    logger.debug(f"{self.counter}/{self.total}: {pdb_id} {chain_id} processed")

    def write_summary(self):
        self.summary.sort(key=lambda x: x[0])  # sort by name
        with open(os.path.join(self.output_dir, f"summary.csv"), 'w') as f:
            f.write("\"feature\", \"P-value\", \"# binding sites\", \"# nonbinding sites\", \"test\"\n") #header
            for line in self.summary:
                f.write(f"{str(line[0])}, {format(line[1], '.4g')}, {str(line[2])}, {str(line[3])}, {str(line[4])}\n")
                #f.write(', '.join(str(v) for v in line) + '\n')

    def welchs_t_test(self, data, results_file, feature):
        with open(results_file, 'w') as f:
            equal_var = False  # todo examine if variance is expected to be the same
            # alpha = 0.05  # todo as script parameter

            binding_data = [x[1] for x in data if x[0] == 1]
            nonBinding_data = [x[1] for x in data if x[0] == 0]

            out = stats.ttest_ind(binding_data, nonBinding_data, axis=0, equal_var=equal_var, nan_policy='raise')
            t_statistic = out[0]
            p_value = out[1]

            p_value = float(p_value)

            f.write(f"Binding residues count: {len(binding_data)}\n")
            f.write(f"Non-binding residues count: {len(nonBinding_data)}\n")
            f.write(f"t-statistic: {t_statistic}\n")
            f.write(f"p-value: {p_value}\n")
            f.write(f"Binding sites -> mean: {np.mean(binding_data)}; variance: {np.var(binding_data)}\n")
            f.write(f"Non-binding sites -> mean: {np.mean(nonBinding_data)}; variance: {np.var(nonBinding_data)}\n")
            f.write(f"Total -> mean: {np.mean(data)}; variance: {np.var(data)}\n")

        return (feature, p_value, len(binding_data), len(nonBinding_data), "Welch's")

    def fischers_exact_test(self, data, results_file, feature):
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

        return (feature, p_value, binding_positive + binding_negative, non_binding_positive + non_binding_negative, "Fisher's exact")

    def chi_squared_test(self, data, results_file, feature):
        with open(results_file, 'w') as f:
            from itertools import groupby
            from operator import itemgetter

            binding_data = [x[1] for x in data if x[0] == 1]
            nonBinding_data = [x[1] for x in data if x[0] == 0]

            sorter = sorted(data, key=itemgetter(0))
            grouper = groupby(sorter, key=itemgetter(0))

            res = {k: list(map(itemgetter(1), v)) for k, v in grouper}

            binding = res[1]
            non_binding = res[0]

            # binding_counts = Counter(binding)
            binding_counts = Counter(binding)
            non_binding_counts = Counter(non_binding)

            keys1 = binding_counts.keys()
            keys2 = non_binding_counts.keys()
            categories = list(set().union(keys1, keys2))

            binding = []
            non_binding = []
            for cat in categories:
                binding.append(binding_counts[cat])
                non_binding.append(non_binding_counts[cat])

            obs = np.array([np.array(binding), np.array(non_binding)])
            chi2, p_value, dof, expected = stats.chi2_contingency(obs)

            p_value = float(p_value)

            f.write(f"Binding counts: {binding_counts}\n")
            f.write(f"Non-binding counts: {non_binding_counts}\n")

            f.write(f"Chi squared: {chi2}\n")
            f.write(f"p-value: {p_value}\n")
            f.write(f"Degrees of freedom: {dof}\n")
            f.write(f"Expected values: {expected}\n")
            #try:
            #    f.write(f"Median binding: {np.median(binding_data)}") #todo smazat
            #    f.write(f"Median nonbinding: {np.median(nonBinding_data)}")
            #    f.write(f"Median: {np.median(data)}")
            #except:
            #    print("error " + feature)

        return (feature, p_value, len(binding_data), len(nonBinding_data), "Chi-squared")