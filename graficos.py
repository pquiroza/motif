import pandas as pd
from pandas import Series
from matplotlib import pyplot
import decimal
import numpy as np


f = open("results/PS00010mejor.fts",'r')
datos = []
for l in f:
    dec = decimal.Decimal(l)
    flo = float(dec)
    datos.append(flo)
in_array = np.around(datos,decimals=3)
#series = pd.read_csv('results/Book4data.csv',header=0)
series = pd.Series(in_array)
series.plot()
pyplot.show()
