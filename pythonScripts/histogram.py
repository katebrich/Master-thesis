import matplotlib.pyplot as plt
from scipy import stats


binding_data = stats.norm.rvs(loc=5, scale=2, size=2000, random_state=1234)
nonBinding_data = stats.norm.rvs(loc=7, scale=2, size=2000, random_state=1234)

#weights_binding = np.ones_like(binding_data)/len(binding_data)
#weights_nonBinding = np.ones_like(nonBinding_data)/len(nonBinding_data)
#plt.hist((binding_data, nonBinding_data), bins=20, weights=(weights_binding, weights_nonBinding), alpha=0.5, label=('Binding sites', 'Non-binding sites'))

plt.hist((binding_data, nonBinding_data), bins=20, alpha=0.5, label=('Binding sites', 'Non-binding sites'))

plt.legend(loc='upper left')
plt.xlabel('Value')
plt.ylabel('Count')
plt.savefig('histogram.svg')