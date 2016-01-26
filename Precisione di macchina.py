# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 18:18:43 2015

@author: pelag
"""

import numpy as np
import matplotlib.pyplot as plt

"""
La funzione numpy.linspace(start, stop, num=50, endpoint=True) restituisce i punti equispaziati tra due estremi
(start e stop) nel numero indicato nel terzo parametro (num = ...). 
Il parametro booleano endpoint(di default impostato a true) consente di includere o meno nel calcolo il punto di stop.

La funzione matplotlib.pyplot.loglog(*args, **kwargs) crea un grafico con scala logaritmica in x e y.
Questa funzione supporta tutti parametri della funzone plot().
"""

#Di seguito vengono impostati i parametri che useremo per tracciare il nostro grafico
#Nell'ordine: 
#1)numero di Nepero
#2)valore di n compreso fra 10^0 e 10^16
#3)la funzione 
e = np.e 
n = np.linspace(pow(10, 0),pow(10, 16), 100000)
y = abs(e-pow((1+1/n),n)) 

#Creo e mostro il grafico
plt.loglog(n,y, color='green')
plt.show() 

def nepero(n):
    e=pow((1+1./n),n)
    return e
n=1e5
w=nepero(n)
print w, 'per n=10^5'
n=1e10
w=nepero(n)
print w, 'per n=10^10'
n=1e11
w=nepero(n)
print w, 'per n=10^11 (la migliore approssimazione del numero di nepero)'
n=1e12
w=nepero(n)
print w, 'per n=10^12'
n=1e13
w=nepero(n)
print w, 'per n=10^13'
n=1e14
w=nepero(n)
print w, 'per n=10^14'
n=1e15
w=nepero(n)
print w, 'per n=10^15'
n=1e16
w=nepero(n)
print w, 'per n=10^16'
