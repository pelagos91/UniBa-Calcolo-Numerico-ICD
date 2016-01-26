# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 09:25:54 2015

@author: pelag
"""

import numpy as np
import scipy.linalg as las
import matplotlib.pylab as plt
import scipy.interpolate as interpolate

def lagrange(xnodi,fnodi,x):
     """
     Funzione che determina in un insieme di punti
     il valore del polinomio interpolante ottenuto
     dalla formula di Lagrange e la funzione di Lebesgue
 
     f,L = lagrange(xnodi,fnodi,x)
 
     Dati di input:
 
     xnodi vettore con i nodi dell'interpolazione
     fnodi vettore con i valori nei nodi
     x   vettore con i punti in cui si vuole
          calcolare il valore del polinomio
          interpolante
     Dati di output:
     f   vettore contenente i valori assunti
          dal polinomio interpolante
     L  vettore contenente i valori assunti dalla
         funzione di lebesgue in x
 
     """
     n=xnodi.shape[0]
     m=x.shape[0]

     l=np.zeros([n,m])
     ind=np.zeros([n,n-1],dtype=int)
     
     ind[0,:]=np.arange(1,n)
    
     for i in range(1,(n)):
        ind[i,0:(n-1)]=np.floor(np.concatenate((np.arange(0,(i)),np.arange(i+1,n))))
     ind[n-1,:]=np.arange(0,n-1)
    

     for i in range(0,n):
         den = np.prod( xnodi[i]-xnodi[(ind[i,:])])
         for j in range(0,m):
             l[i,j]=np.prod( x[j] - xnodi[ind[i,:]] )/den
    


     y = np.sum(np.dot(np.diag(fnodi),l),axis=0) #Vettore di punti del olinomio interpolante
     

     L=np.sum(abs(l),axis=0) #Funzione di Lebesgue

     return  y,L

def potenze(xnodi,fnodi,xx):
     """
     Funzione che determina in un insieme di punti
     il valore del polinomio interpolante ottenuto
     dalla formula di Lagrange e la funzione di Lebesgue
 
     f = potenze(xnodi,fnodi,x)
 
     Dati di input:
 
     xnodi vettore con i nodi dell'interpolazione
     fnodi vettore con i valori nei nodi
     x   vettore con i punti in cui si vuole
          calcolare il valore del polinomio
          interpolante
     Dati di output:
     f   vettore contenente i valori assunti
          dal polinomio interpolante
 
     """
     n=xnodi.shape[0]
     A=np.zeros([n,n]) #Matrice di Vandermonde
     for j in range(0,n):
        A[:,j]=xnodi**(j)

     p = las.solve(A,fnodi)
     f = np.polyval(p[np.arange(n-1,-1,-1)],xx)
     condA = np.linalg.cond(A, np.inf)
     return f, condA
     
def  cheby(a,b,n):
     """
     Nodi di Chebyshev
     """
     c =  (a + b + (b-a)*np.cos((2*(np.arange(0,n+1))+1)*np.pi/(2*(n+1))))/2
     return c
     
def runge(x):
    """
    Funzione di Runge
    """
    y=1/(1+25*x**2)
    return y
    
def plotinterp(ftype,a,b,n,tn,baseType):
    """
    plot del polinomio interpolante di grado n la funzione
    f in [a,b] usando n+1 nodi equidistanti se tn=0
             n+1 nodi di Chebyshev se tn=1
             
    ftype indica quale delle due funzioni utilizzare:
    - 0 per la funzione coseno
    - 1 per la funzione di Runge
    
    baseType indica il tipo di base usata:
    - 0 per la base di Lagrange
    - 1 per la base delle potenze
    """
    
    if (tn==0):
       xnodi = np.linspace(a,b,n+1)
    else:
       xnodi = cheby(a,b,n) 
    
    if (ftype==0):
       fname='f=cos(x)'
       f=np.cos
       fnodi = f(xnodi)
       xx = np.linspace(a,b,500)
       ye = f(xx)
    else:
       fname='g=1/(1+25*(x**2))'
       fnodi = gfunction(xnodi)
       xx = np.linspace(a,b,500)
       ye = gfunction(xx)
    
    if(baseType==0):
       fi, L = lagrange(xnodi,fnodi,xx)
       Lc = las.norm(L, np.inf) 
    else:
       fi, condA = potenze(xnodi, fnodi, xx)
  
        
    error = np.max(np.abs(fi-ye)) 
     
    if(baseType==0):
       plt.figure(1)
       plt.cla()
       plt.title('Polinomio interpolante per la funzione %s con n= %i'%(fname, n))
       plt.plot(xx,fi,xnodi,fnodi,'o',xx,ye,'--') 
       plt.figure(2)
       plt.cla()
       plt.plot(xx,L)
       plt.show()
    else:
       plt.figure(1)
       plt.cla()
       plt.title('Polinomio interpolante per la funzione %s con n= %i'%(fname, n))
       plt.plot(xx,fi,xnodi,fnodi,'o',xx,ye,'--') 
       plt.show()
     
    if(baseType==0): 
        return error, Lc
    else:
        return error, condA
    
def splineFunction(xnodi,fnodi, xx, fType, sType):
    if(sType==0):
        s1 = interpolate.interp1d(xnodi, fnodi, 'linear')
        sname='lineare'
    else:
        s1 = interpolate.interp1d(xnodi, fnodi, 'cubic')
        sname='cubica'
        
    if(fType==0):
        fname='f=cos(x)'
    else:
        fname='g=1/(1+25*(x**2))' 
        
    ys = s1(xx)
    error.append(np.max(np.abs(ys-yy)))
    plt.figure(i)
    plt.cla()
    plt.title('Spline %s per la funzione %s' %(sname,fname))
    plt.plot(xx,ys,xnodi,fnodi,'o',xx,yy,'--')
    plt.show()

def gfunction(x_variable):
     return 1/(1+25*(x_variable**2))
 
print("________________________________")
print("|   POTENZE NODI EQUIDISTANTI  |")
print("|______________________________|")

#prima funzione per n=4
a=0
b=2
n=4

numeroCondA = np.zeros([4,2])#Vettore in cui avviene lo store dei numeri di condizione della matrice A
numeroCondA[0]= plotinterp(0,a,b,n,0,1)

#prima funzione per n=16
n=16
numeroCondA[1] = plotinterp(0,a,b,n,0,1)

#seconda funzione per n=4
a = -1
b = 1
n=4
numeroCondA[2] = plotinterp(1,a,b,n,0,1)

#seconda funzione per n=16
n=16
numeroCondA[3] = plotinterp(1,a,b,n,0,1)

for i in range(0,4):
    print numeroCondA[i]
    

print("________________________________")
print("|  LAGRANGE NODI EQUIDISTANTI  |")
print("|______________________________|")
erroreInterpE = np.zeros([4,2])#Vettore in cui avviene lo store degli errori per nodi equidistanti
#prima funzione per n=4
a = 0
b = 2
n=4
erroreInterpE[0] = plotinterp(0,a,b,n,0,0)
#prima funzione per n=16
n=16
erroreInterpE[1] = plotinterp(0,a,b,n,0,0)

#seconda funzione per n=4
a = -1
b = 1
n = 4
erroreInterpE[2] = plotinterp(1,a,b,n,0,0)
#seconda funzione per n=16
n=16
erroreInterpE[3] = plotinterp(1,a,b,n,0,0)

print("________________________________")
print("|   POTENZE NODI DI CHEBYCHEV  |")
print("|______________________________|")

#prima funzione per n=4
a=0
b=2
n=4
numeroCondAC = np.zeros([4,2])#Vettore in cui avviene lo store dei numeri di condizione della matrice A
numeroCondAC[0] = plotinterp(0,a,b,n,1,1)

#prima funzione per n=16
n=16
numeroCondAC[1] = plotinterp(0,a,b,n,1,1)

#seconda funzione per n=4
a = -1
b = 1
n=4
numeroCondAC[2] = plotinterp(1,a,b,n,1,1)

#seconda funzione per n=16
n=16
numeroCondAC[3] = plotinterp(1,a,b,n,1,1)

for i in range(0,4):
    print numeroCondA[i]



print("_________________________________")
print("|  LAGRANGE NODI DI CHEBYCHEV   |")
print("|_______________________________|")

erroreInterpC = np.zeros([4,2])#Vettore in cui avviene lo store degli errori per nodi di Chebychev

#prima funzione per n=4
a = 0
b = 2
n=4
erroreInterpC[0] = plotinterp(0,a,b,n,1,0)
#prima funzione per n=16
n=16
erroreInterpC[1] = plotinterp(0,a,b,n,1,0)

#seconda funzione per n=4
a = -1
b = 1
n = 4
erroreInterpC[2] = plotinterp(1,a,b,n,1,0)
#seconda funzione per n=16
n=16
erroreInterpC[3] = plotinterp(1,a,b,n,1,0)

# interplazione con le funzioni spline
#funzione f 
#4 nodi
f=np.cos
xx = np.linspace(a,b,200)
yy = f(xx)
error=[]
n=4
xnodi = np.linspace(a,b,n+1)
fnodi = f(xnodi)
splineFunction(xnodi, fnodi, xx, 0, 0)
splineFunction(xnodi, fnodi, xx, 0, 1)
#16 nodi
n=16
xnodi = np.linspace(a,b,n+1)
fnodi = f(xnodi) 
splineFunction(xnodi, fnodi, xx, 0, 0)
splineFunction(xnodi, fnodi, xx, 0, 1) 

#funzione g
#4 nodi

xx = np.linspace(a,b,200)
yy = gfunction(xx)
n=4
xnodi = np.linspace(a,b,n+1)
fnodi = gfunction(xnodi)
splineFunction(xnodi, fnodi, xx, 1, 0)
splineFunction(xnodi, fnodi, xx, 1, 1)
#16 nodi
n=16
xnodi = np.linspace(a,b,n+1)
fnodi = gfunction(xnodi) 
splineFunction(xnodi, fnodi, xx, 1, 0)
splineFunction(xnodi, fnodi, xx, 1, 1) 
print   
print "ERRORE DELLA SPLINE"  

for i in range(0,8):
    if (i<4):
       print 'Funzione f=cos(x)'
    else:
       print 'Funzione di Runge'
    print error[i]





print("_________________________________")
print("|       CONFRONTO ERRORI        |")
print("|          f = cos(x)           |")
print("|_______________________________|")

print("ERRORE             |   NUMERO DI CONDIZIONE VANDERMONDE/COSTANTE DI LEBESGUE")
print("-----------------------------------------------------------------------------")
print("n = 4")
print("-----------------------------------------------------------------------------")
print numeroCondA[0],"Base delle POTENZE e nodi EQUIDISTANTI"
print numeroCondAC[0],"Base delle POTENZE e nodi di CHEBYCHEV"
print erroreInterpE[0],"Base di LAGRANGE e nodi EQUIDISTANTI"
print erroreInterpC[0],"Base di LAGRANGE e nodi di CHEBYCHEV"
print("-----------------------------------------------------------------------------")
print("n = 16")
print("-----------------------------------------------------------------------------")
print numeroCondA[1],"Base delle POTENZE e nodi EQUIDISTANTI"
print numeroCondAC[1],"Base delle POTENZE e nodi di CHEBYCHEV"
print erroreInterpE[1],"Base di LAGRANGE e nodi EQUIDISTANTI"
print erroreInterpC[1],"Base di LAGRANGE e nodi di CHEBYCHEV"

print("_________________________________")
print("|       CONFRONTO ERRORI        |")
print("|       g=1/(1+25*(x**2))       |")
print("|_______________________________|")
print("ERRORE             |   NUMERO DI CONDIZIONE VANDERMONDE/COSTANTE DI LEBESGUE")
print("-----------------------------------------------------------------------------")
print("n = 4")
print("-----------------------------------------------------------------------------")
print numeroCondA[2],"Base delle POTENZE e nodi EQUIDISTANTI"
print numeroCondAC[2],"Base delle POTENZE e nodi di CHEBYCHEV"
print erroreInterpE[2],"Base di LAGRANGE e nodi EQUIDISTANTI"
print erroreInterpC[2],"Base di LAGRANGE e nodi di CHEBYCHEV"
print("-----------------------------------------------------------------------------")
print("n = 16")
print("-----------------------------------------------------------------------------")
print numeroCondA[3],"Base delle POTENZE e nodi EQUIDISTANTI"
print numeroCondAC[3],"Base delle POTENZE e nodi di CHEBYCHEV"
print erroreInterpE[3],"Base di LAGRANGE e nodi EQUIDISTANTI"
print erroreInterpC[3],"Base di LAGRANGE e nodi di CHEBYCHEV"
