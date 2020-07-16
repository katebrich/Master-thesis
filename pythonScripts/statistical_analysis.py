import pandas as pd
from scipy import stats
import numpy as np

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
    #todo check dataframe

    tab = pd.crosstab(dataframe['ligand_binding_sites'], dataframe['feature']) #todo indexovat jinak

    print(tab)

    #todo opravit..kdyz neni ani jedna 1, pada
    binding_positive = tab.loc[1].loc[1]
    binding_negative = tab.loc[1].loc[0]
    non_binding_positive = tab.loc[0].loc[1]
    non_binding_negative = tab.loc[0].loc[0]

    odds_ratio, p_value = stats.fisher_exact([[binding_positive, binding_negative],
                                              [non_binding_positive, non_binding_negative]])

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