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

        #print(res)


        #counts = Counter(data)

        binding = res[1]
        non_binding = res[0]

        print(binding)
        print(non_binding)

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


def compute_AA_frequencies(data):
    binding_data = [x[1] for x in data if x[0] == 1]
    nonBinding_data = [x[1] for x in data if x[0] == 0]

    print(f"Binding residues count: {len(binding_data)}")
    print(f"Nonbinding residues count: {len(nonBinding_data)}")

    binding_AAs_count = Counter(binding_data)
    nonBinding_AAs_count = Counter(nonBinding_data)

    binding_AAs_frequencies = {k: v / len(binding_data) for k, v in binding_AAs_count.items()}
    nonBinding_AAs_frequencies = {k: v / len(nonBinding_data) for k, v in nonBinding_AAs_count.items()}

    binding_AAs_frequencies_sorted = {}
    nonBinding_AAs_frequencies_sorted = {}

    AAs = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

    for AA in AAs:
        freq = binding_AAs_frequencies.get(AA, 0) #default is 0
        binding_AAs_frequencies_sorted[AA] = freq
        freq = nonBinding_AAs_frequencies.get(AA, 0)  # default is 0
        nonBinding_AAs_frequencies_sorted[AA] = freq

    print("Binding AA frequencies:")
    print(binding_AAs_frequencies_sorted)

    print("Non-binding AA frequencies:")
    print(nonBinding_AAs_frequencies_sorted)

    #plot the results
    barWidth = 0.33

    # set height of bar
    bars1 = binding_AAs_frequencies_sorted.values()
    bars2 = nonBinding_AAs_frequencies_sorted.values()

    # Set position of bar on X axis
    r1 = np.arange(len(bars1))
    r2 = [x + barWidth for x in r1]

    # Make the plot
    plt.bar(r1, bars1, color='g', width=barWidth, edgecolor='white', label='binding')
    plt.bar(r2, bars2, color='r', width=barWidth, edgecolor='white', label='non-binding')

    # Add xticks on the middle of the group bars
    plt.xlabel('Amino acid', fontweight='bold')
    plt.ylabel('Frequency', fontweight='bold')
    plt.xticks([r + barWidth for r in range(len(bars1))], AAs)

    # Create legend & Show graphic
    plt.legend()
    plt.show()
