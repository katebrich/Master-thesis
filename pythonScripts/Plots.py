import operator
import os

import numpy as np
from collections import Counter
from matplotlib import pyplot
from scipy import stats


def plot_frequencies(pairs):
    binding_data = [x[1] for x in pairs if x[0] == 1]
    nonBinding_data = [x[1] for x in pairs if x[0] == 0]

    binding_counts = Counter(binding_data)
    nonBinding_counts = Counter(nonBinding_data)

    binding_frequencies = {k: v / len(binding_data) for k, v in binding_counts.items()}
    nonBinding_frequencies = {k: v / len(nonBinding_data) for k, v in nonBinding_counts.items()}

    binding_AAs_frequencies_sorted = {}
    nonBinding_AAs_frequencies_sorted = {}

    AAs = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

    for AA in AAs:
        freq = binding_frequencies.get(AA, 0) #default is 0
        binding_AAs_frequencies_sorted[AA] = freq
        freq = nonBinding_frequencies.get(AA, 0)  # default is 0
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
    pyplot.bar(r1, bars1, color='g', width=barWidth, edgecolor='white', label='binding')
    pyplot.bar(r2, bars2, color='r', width=barWidth, edgecolor='white', label='non-binding')

    # Add xticks on the middle of the group bars
    pyplot.xlabel('Amino acid', fontweight='bold')
    pyplot.ylabel('Frequency', fontweight='bold')
    pyplot.xticks([r + barWidth for r in range(len(bars1))], AAs)

    # Create legend & Show graphic
    pyplot.legend()
    pyplot.show()

def plot_binding_nonbinding_ratios(pairs, output_path):
    pyplot.clf()  # todo smazat i kdyz chyba!

    binding_data = [x[1] for x in pairs if x[0] == 1]
    nonBinding_data = [x[1] for x in pairs if x[0] == 0]

    binding_counts = Counter(binding_data)
    nonBinding_counts = Counter(nonBinding_data)

    #print(binding_counts)
    #print(nonBinding_counts)

    ratios={}

    for k in nonBinding_counts.keys():
        ratio = binding_counts[k] / nonBinding_counts[k]
        ratios[k] = ratio

    #binding_AAs_frequencies_sorted = {}
    #nonBinding_AAs_frequencies_sorted = {}

    #AAs = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

    #for AA in AAs:
    #    freq = binding_frequencies.get(AA, 0) #default is 0
    #    binding_AAs_frequencies_sorted[AA] = freq
    #    freq = nonBinding_frequencies.get(AA, 0)  # default is 0
    #    nonBinding_AAs_frequencies_sorted[AA] = freq
    sorted_by_ratios = sorted(ratios.items(), key=operator.itemgetter(1), reverse=True)

    #print("Binding AA frequencies:")
    #print(binding_AAs_frequencies_sorted)

    #print("Non-binding AA frequencies:")
    #print(nonBinding_AAs_frequencies_sorted)

    #print(ratios.items())

    #plot the results
    barWidth = 0.45

    # set height of bar
    bars = [x[1] for x in sorted_by_ratios]
    #bars2 = nonBinding_AAs_frequencies_sorted.values()

    # Set position of bar on X axis
    r1 = np.arange(len(bars))
    r2 = [x + barWidth/2 for x in r1]

    # Make the plot
    pyplot.bar(r2, bars, color='g', width=barWidth, edgecolor='white')
    #plt.bar(r2, bars2, color='r', width=barWidth, edgecolor='white', label='non-binding')

    # Add xticks on the middle of the group bars
    pyplot.xlabel('Amino acid', fontweight='bold')
    pyplot.ylabel('Binding/nonbinding ratio', fontweight='bold')
    pyplot.xticks([r + barWidth/2 for r in range(len(bars))], [x[0] for x in sorted_by_ratios])

    #Draw total binding/nonbinding ratio for comparison
    ratio = len(binding_data) / len(nonBinding_data)
    pyplot.hlines(ratio, xmin=0 , xmax=len(bars) - 1 + barWidth/2, colors='k', linestyles='dashed', label='total binding/non-binding') #todo xmin, xmax

    # Create legend & Show graphic
    #plt.legend()
    #pyplot.show()

    pyplot.savefig(os.path.join(output_path, "plot_ratio.png"))
    pyplot.clf() #todo smazat i kdyz chyba!


def plot_histogram(pairs, bins, output_path):

    #binding_data = stats.norm.rvs(loc=5, scale=2, size=2000, random_state=1234)
    #nonBinding_data = stats.norm.rvs(loc=7, scale=2, size=2000, random_state=1234)
    binding_data = [x[1] for x in pairs if x[0] == 1]
    nonBinding_data = [x[1] for x in pairs if x[0] == 0]

    weights_binding = np.ones_like(binding_data)/len(binding_data)
    weights_nonBinding = np.ones_like(nonBinding_data)/len(nonBinding_data)

    pyplot.clf()  # todo
    #pyplot.hist((binding_data, nonBinding_data), bins=20, alpha=0.5, label=('Binding sites', 'Non-binding sites'))
    pyplot.hist((binding_data, nonBinding_data), bins=bins, weights=(weights_binding, weights_nonBinding), alpha=0.5,
             label=('Binding sites', 'Non-binding sites'))

    pyplot.legend(loc='upper left')
    pyplot.xlabel('Value')
    pyplot.ylabel('Density')
    pyplot.savefig(os.path.join(output_path,f"histogram_bins{bins}.svg"))
    pyplot.clf()  # todo smazat i kdyz chyba!