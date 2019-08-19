import pandas as pd
from pandas import Series
from matplotlib import pyplot
from matplotlib import animation
import decimal
import numpy as np
import sys


archivo = sys.argv[1]

while True:
    f = open(archivo,'r')
    datos = []
    for l in f:
        dec = decimal.Decimal(l)
        flo = float(dec)
        datos.append(flo)
    in_array = np.around(datos,decimals=3)
    series = pd.Series(in_array)
    series.plot()
    pyplot.show(block=False)
    pyplot.pause(120)
    pyplot.close()
