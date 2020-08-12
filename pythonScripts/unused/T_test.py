from scipy import stats
import numpy as np

#todo input data

equal_var = True #todo examine if variance is expected to be the same
alpha = 0.05 #todo as script parameter

#experimental data
#loc = mean; scale = standard deviation; size = number of samples
binding_data = stats.norm.rvs(loc=5, scale=10, size=1000)
nonBinding_data = stats.norm.rvs(loc=6, scale=10, size=2000)

out = stats.ttest_ind(binding_data, nonBinding_data, axis=0, equal_var=equal_var, nan_policy='raise')
t_statistic = out[0]
p_value = out[1]


print(f"t-statistic: {t_statistic}")
print(f"p-value: {p_value}")
print(f"Binding sites -> mean: {np.mean(binding_data)}; variance: {np.var(binding_data)}")
print(f"Non-binding sites -> mean: {np.mean(nonBinding_data)}; variance: {np.var(nonBinding_data)}")

#The null hypothesis....2 samples have identical average.
if (p_value < alpha):
    print(f"The null hypothesis was rejected at significance level of {alpha}. The samples have different means with probability {(1-alpha)}.")
else:
    print(f"The null hypothesis cannot be rejected at significance level of {alpha}")

