# https://www.statology.org/seaborn-stacked-bar-plot/

import pandas as pd

#create DataFrame
df = pd.DataFrame({'Day': ['Mon', 'Tue', 'Wed', 'Thur', 'Fri'],
                   'Morning': [44, 46, 49, 59, 54],
                   'Evening': [33, 46, 50, 49, 60]})

print(df)
# exit()

import matplotlib.pyplot as plt
import seaborn as sns

#set seaborn plotting aesthetics
sns.set(style='white')

#create stacked bar chart
df.set_index('Day').plot(kind='bar', stacked=True, color=['steelblue', 'red'])

#add overall title
plt.title('Customers by Time & Day of Week', fontsize=16)

#add axis titles
plt.xlabel('Day of Week')
plt.ylabel('Number of Customers')

#rotate x-axis labels
plt.xticks(rotation=45)

plt.show()