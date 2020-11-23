import operator
import os

import numpy as np
from collections import Counter
from matplotlib import pyplot, transforms
from matplotlib.ticker import FormatStrFormatter


def plot_binding_nonbinding_ratios(binding_data, nonbinding_data, output_file, sort_item=0):
    pyplot.clf()
    binding_data = [str(i) for i in binding_data]
    nonbinding_data = [str(i) for i in nonbinding_data]
    binding_counts = Counter(binding_data)
    nonBinding_counts = Counter(nonbinding_data)

    ratios={}
    for k in nonBinding_counts.keys():
        ratio = binding_counts[k] / nonBinding_counts[k]
        ratios[k] = ratio

    if (sort_item == 0): #sort by name
        sorted_ratios = sorted(ratios.items(), key=operator.itemgetter(0))
    elif (sort_item == 1): #sort by ratios values
        sorted_ratios = sorted(ratios.items(), key=operator.itemgetter(1), reverse=True)
    else: #no sorting
        sorted_ratios = ratios.items()

    barWidth = 0.45

    # set height of bar
    bars = [x[1] for x in sorted_ratios]
    # Set position of bar on X axis
    r1 = np.arange(len(bars))
    r2 = [x + barWidth/2 for x in r1]

    fig, ax = pyplot.subplots()
    pyplot.bar(r2, bars, color='g', width=barWidth, edgecolor='white')

    pyplot.xlabel('Value', fontweight='bold')
    pyplot.ylabel('Binding/nonbinding ratio', fontweight='bold')
    # Add xticks in the middle of bars
    pyplot.xticks([r + barWidth/2 for r in range(len(bars))], [x[0] for x in sorted_ratios])

    #Draw total binding/nonbinding ratio for comparison
    ratio = len(binding_data) / len(nonbinding_data)
    #pyplot.hlines(ratio, xmin=0 , xmax=len(bars) - 1 + barWidth, colors='r', linestyles='dashed', label='total ratio') #todo xmin, xmax
    ax.axhline(y=ratio, color="red", linestyle='--' )

    trans = transforms.blended_transform_factory(
        ax.get_yticklabels()[0].get_transform(), ax.transData)
    ax.text(0, ratio, "{:.3f}".format(ratio), color="red", transform=trans,
            ha="right", va="center")

    #round y axis labels to 3 places
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

    #save figure and clear
    pyplot.savefig(output_file)
    pyplot.clf()


def plot_histogram(binding_data, nonbinding_data, bins, output_file):
    pyplot.clf()  #
    weights_binding = np.ones_like(binding_data)/len(binding_data)
    weights_nonBinding = np.ones_like(nonbinding_data) / len(nonbinding_data)
    pyplot.hist((binding_data, nonbinding_data), bins=bins, weights=(weights_binding, weights_nonBinding), alpha=0.5,
                label=('Binding sites', 'Non-binding sites'))

    pyplot.legend(loc='upper left')
    pyplot.xlabel('Value')
    pyplot.ylabel('Density')
    pyplot.savefig(output_file)
    pyplot.clf()