import matplotlib.pyplot as plt
from pandasql import sqldf
import pandas as pd

# Read the results from the CSV file
results = pd.read_csv("../results/results.csv")

# Use pandasql to group by preconditioner and p_value, finding the minimum number of iterations
queried_data = sqldf('''
    SELECT mesh, p_value, preconditioner, MIN(cond_number) as CN
    FROM results
    GROUP BY preconditioner, p_value
    ORDER BY preconditioner
''')

# Plotting
plt.figure()
for preconditioner in queried_data['preconditioner'].unique():
    subset = queried_data[queried_data["preconditioner"] == preconditioner]
    plt.plot(subset['p_value'], subset['CN'], label=preconditioner)

plt.xlabel('p')
plt.ylabel('Condition number')
plt.yscale('log')  # Adjust the scale if needed
plt.title('Number of Iterations vs Sphere Exponent for Different Preconditioners')
plt.grid()
plt.legend()
plt.savefig('../results/plot/cond_numb_vs_p')
plt.close()
