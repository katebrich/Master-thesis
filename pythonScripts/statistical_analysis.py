import pandas as pd
from scipy import stats
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt

alpha = 0.05 #todo

def welchs_t_test(data):
    equal_var = False  # todo examine if variance is expected to be the same
    alpha = 0.05  # todo as script parameter

    binding_data = [x[1] for x in data if x[0] == 1]
    nonBinding_data = [x[1] for x in data if x[0] == 0]

    out = stats.ttest_ind(binding_data, nonBinding_data, axis=0, equal_var=equal_var, nan_policy='raise')
    t_statistic = out[0]
    p_value = out[1]

    print(f"Binding residues count: {len(binding_data)}")
    print(f"Non-binding residues count: {len(nonBinding_data)}")
    print(f"t-statistic: {t_statistic}")
    print(f"p-value: {p_value}")
    print(f"Binding sites -> mean: {np.mean(binding_data)}; variance: {np.var(binding_data)}")
    print(f"Non-binding sites -> mean: {np.mean(nonBinding_data)}; variance: {np.var(nonBinding_data)}")

    # The null hypothesis....2 samples have identical average.
    if (p_value < alpha):
        print(
            f"The null hypothesis was rejected at significance level of {alpha}. The samples have different means with probability {(1 - alpha)}.")
    else:
        print(f"The null hypothesis cannot be rejected at significance level of {alpha}")


def fischers_exact_test(data):
    #todo check data argument

    counts = Counter(data)

    #todo opravit..kdyz neni ani jedna 1, pada
    binding_positive = counts[(1,1)]
    binding_negative = counts[(1,0)]
    non_binding_positive = counts[(0,1)]
    non_binding_negative = counts[(0,0)]

    odds_ratio, p_value = stats.fisher_exact([[binding_positive, binding_negative],
                                              [non_binding_positive, non_binding_negative]])

    print(f"Binding positive count: {binding_positive}")
    print(f"Binding negative count: {binding_negative}")
    print(f"Nonbinding positive count: {non_binding_positive}")
    print(f"Nonbinding negative count: {non_binding_negative}")

    print(f"odds ratio: {odds_ratio}")
    print(f"p-value: {p_value}")
    print(f"Binding sites -> positives ratio: {binding_positive / (binding_positive + binding_negative)}")
    print(
        f"Non-binding sites -> positives ratio: {non_binding_positive / (non_binding_positive + non_binding_negative)}")

    # The null hypothesis....proportions of the feature are the same for binding and non-binding sites.
    if (p_value < alpha):
        print(
            f"The null hypothesis was rejected at significance level of {alpha}. The samples have different proportions with probability {(1 - alpha)}.")
    else:
        print(f"The null hypothesis cannot be rejected at significance level of {alpha}")

def chi_squared_test():
    stats.chisquare()


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
