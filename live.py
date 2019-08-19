from pylab import * # import matplotlib before drawnow
from drawnow import drawnow, figure
from time import sleep
import numpy as np
import pandas as pd
from pandas import Series
import decimal

def draw_fig_real():
    imshow(data)

f = open("results/PS00010mejor7.fts",'r')
datos = []
for l in f:
    dec = decimal.Decimal(l)
    flo = float(dec)
    datos.append(flo)
data = np.around(datos,decimals=3)
f.close()

figure()
drawnow(draw_fig_real)
while(5>8):
    f = open("results/PS00010mejor7.fts",'r')
    datos = []
    for l in f:
        dec = decimal.Decimal(l)
        flo = float(dec)
        datos.append(flo)
    data = np.around(datos,decimals=3)
    f.close()
    sleep(0.5)
    drawnow(draw_fig_real)
