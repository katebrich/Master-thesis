import pandas as pd
from scipy import stats

alpha = 0.05 #todo

def fischers_exact_test(dataframe):
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