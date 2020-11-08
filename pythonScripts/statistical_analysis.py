import pandas as pd
from scipy import stats
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt

alpha = 0.05 #todo

def welchs_t_test(data, results_file):
    with open(results_file, 'w') as f:
        equal_var = False  # todo examine if variance is expected to be the same
        #alpha = 0.05  # todo as script parameter

        binding_data = [x[1] for x in data if x[0] == 1]
        nonBinding_data = [x[1] for x in data if x[0] == 0]

        out = stats.ttest_ind(binding_data, nonBinding_data, axis=0, equal_var=equal_var, nan_policy='raise')
        t_statistic = out[0]
        p_value = out[1]

        f.write(f"Binding residues count: {len(binding_data)}\n")
        f.write(f"Non-binding residues count: {len(nonBinding_data)}\n")
        f.write(f"t-statistic: {t_statistic}\n")
        f.write(f"p-value: {p_value}\n")
        f.write(f"Binding sites -> mean: {np.mean(binding_data)}; variance: {np.var(binding_data)}\n")
        f.write(f"Non-binding sites -> mean: {np.mean(nonBinding_data)}; variance: {np.var(nonBinding_data)}\n")



def fischers_exact_test(data, results_file):
    with open(results_file, 'w') as f:
        counts = Counter(data)

        #todo opravit..kdyz neni ani jedna 1, pada
        binding_positive = counts[(1,1)]
        binding_negative = counts[(1,0)]
        non_binding_positive = counts[(0,1)]
        non_binding_negative = counts[(0,0)]

        odds_ratio, p_value = stats.fisher_exact([[binding_positive, binding_negative],
                                                  [non_binding_positive, non_binding_negative]])

        f.write(f"Binding positive count: {binding_positive}\n")
        f.write(f"Binding negative count: {binding_negative}\n")
        f.write(f"Nonbinding positive count: {non_binding_positive}\n")
        f.write(f"Nonbinding negative count: {non_binding_negative}\n")

        f.write(f"odds ratio: {odds_ratio}\n")
        f.write(f"p-value: {p_value}\n")
        f.write(f"Binding sites -> positives ratio: {binding_positive / (binding_positive + binding_negative)}\n")
        f.write(
            f"Non-binding sites -> positives ratio: {non_binding_positive / (non_binding_positive + non_binding_negative)}\n")


def chi_squared_test(data, results_file):
    with open(results_file, 'w') as f:
        from itertools import groupby
        from operator import itemgetter

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

        f.write(f"Binding counts: {binding_counts}\n")
        f.write(f"Non-binding counts: {non_binding_counts}\n")

        f.write(f"Chi squared: {chi2}\n")
        f.write(f"p-value: {p_value}\n")
        f.write(f"Degrees of freedom: {dof}\n")
        f.write(f"Expected values: {expected}\n")


