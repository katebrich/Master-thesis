from scipy import stats

#todo input data

alpha = 0.05 #todo script parameter

binding_positive, binding_negative = 8, 7
non_binding_positive, non_binding_negative = 14, 1

odds_ratio, p_value = stats.fisher_exact([[binding_positive, binding_negative],
                                        [non_binding_positive, non_binding_negative]])

print(f"odds ratio: {odds_ratio}")
print(f"p-value: {p_value}")
print(f"Binding sites -> positives ratio: {binding_positive/(binding_positive+binding_negative)}")
print(f"Non-binding sites -> positives ratio: {non_binding_positive/(non_binding_positive+non_binding_negative)}")

#The null hypothesis....proportions of the feature are the same for binding and non-binding sites.
if (p_value < alpha):
    print(f"The null hypothesis was rejected at significance level of {alpha}. The samples have different proportions with probability {(1-alpha)}.")
else:
    print(f"The null hypothesis cannot be rejected at significance level of {alpha}")

