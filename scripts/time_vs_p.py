import matplotlib.pyplot as plt
from pandasql import sqldf
import pandas as pd

results = pd.read_csv("../results/results.csv")

# noinspection SqlNoDataSourceInspection,SqlDialectInspection
queried_data = sqldf('''
SELECT mesh, p_value, preconditioner, MIN(time) as MinTime
FROM results
GROUP BY preconditioner, p_value
ORDER BY preconditioner
'''
                     )

plt.figure()
for preconditioner in queried_data['preconditioner'].unique():
    subset = queried_data[queried_data["preconditioner"] == preconditioner]
    plt.plot(subset['p_value'], subset['MinTime'], label=preconditioner)

plt.xlabel('p')
plt.ylabel('time [ms]')
plt.yscale('log')
plt.title('Execution time in ms vs sphere exponent for different preconditioners')
plt.grid()
plt.legend()
plt.savefig('../results/plot/time_vs_p')
plt.close()
