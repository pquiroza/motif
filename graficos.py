import pandas as pd
from pandas import Series
from matplotlib import pyplot
series = pd.read_csv('results/Book4data.csv',header=0)
print(series)
series.plot()
pyplot.show()
