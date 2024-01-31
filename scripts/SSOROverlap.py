import matplotlib.pyplot as plt
from pandasql import sqldf
import pandas as pd

results = pd.read_csv("../results/results.csv")
results = results.dropna()



# noinspection SqlNoDataSourceInspection,SqlDialectInspection
queried_data = sqldf('''
SELECT p_value, MIN(cond_number) as MinTime, overlap
FROM results
WHERE mesh = 'mesh-cube-20' AND preconditioner = 'ssor'  AND p_value < 8 
GROUP BY p_value, overlap
ORDER BY p_value
'''
                     )

print(queried_data)

plt.figure()
for overlap in queried_data['overlap'].unique():
    print(queried_data[queried_data["overlap"] == overlap])
    plt.plot(queried_data[queried_data["overlap"] == overlap]['p_value'], queried_data[queried_data["overlap"] == overlap]['MinTime'], label=f'overlap = {overlap}')

plt.xlabel('p')
plt.ylabel('condition number')
plt.yscale('log')
plt.title('p_value vs condition number')
plt.grid()
plt.legend()
plt.savefig('../results/plot/SSORoverlapCondNum')
plt.close()
