# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 21:06:34 2015

@author: pelag
"""
import numpy as np

def metodoPotenze(A,z,tol):
    
    
    z=A*z
    Az=A*z
    zt=np.transpose(Az)
    aut_1=(zt*Az)/(zt*z)
    
    z=A*z
    Az=A*z  
    aut_2=aut_1
    aut_1=(zt*Az)/(zt*z)
    count=2
    
    while((abs(aut_1-aut_2)/(1+abs(aut_1)))>tol and count<100):
         error=abs(aut_1-aut_2)/(1+abs(aut_1))      
         z=A*z
         Az=A*z    
         aut_2=aut_1   
         aut_1=(zt*Az)/(zt*z)
         count=count+1
    print 'count =', count    
    print 'errore =', error   
    z=z/np.linalg.norm(z) #Normalizzo il vettore per poterlo confrontare con l'output del secondo metodo
    return aut_1, z
    
def metodoPotenzeNorm(A,z,tol): 
    
    y=z/np.linalg.norm(z)   #Inizio calcolando y0 
    z=A*y                  #Calcolo z1
    Ay=A*y                  #Calcolo il prodotto che inserirÃ² nel numeratore  
   
    aut_1=(y.T*Ay)/(y.T*y)    #Calcolo la prima approssimazione dell'autovalore dominante
    y=z/np.linalg.norm(z)   #Calcolo y1
    
  
    z=A*y                  
    Ay=A*y                 
  
    yt=np.transpose(Ay)  
    aut_2=aut_1             #aut_2 sarÃ  cosÃ¬ uguale a lambda_k-1
    aut_1=(yt*Ay)/(yt*y)   #Nuova approssimazione dell'autovalore dominante
    y=z/np.linalg.norm(z)   
    count=2    
    
    while((abs(aut_1-aut_2)/(1+abs(aut_1)))>tol and count<100):
         error=abs(aut_1-aut_2)/(1+abs(aut_1))
         z=A*y                  
         Ay=A*y                 
        
         yt=np.transpose(Ay)  
         aut_2=aut_1            #aut_2 sarÃ  cosÃ¬ uguale a lambda_k-1
         aut_1=(yt*Ay)/(yt*y)   #Nuova approssimazione dell'autovalore dominante
         y=z/np.linalg.norm(z)   
         count=count+1
    
    print 'count =', count
    print 'errore =', error
    return aut_1, y
    
A = np.matrix([
[4., 1., 0., 0., 0.,], 
[1., 4., 1., 0., 0.,],
[0., 1., 4, 1., 0.,],
[0., 0., 1., 4., 1.,],
[0., 0., 0., 1., 4.,],
])

B = np.matrix([
[-149., -50., -154.,], 
[537., 180., 546.,],
[-27., -9., 25.,],
])

C = np.matrix ([
[63.4000, 20.8000, 64.4000,],
[-180.3000, -59.1000, -184.8000,],
[-2.7000, -0.9000, -2.2000,],
])

z = np.matrix([1,1,1,1,1]).T
z3= np.matrix([1,1,1]).T


tol=1e-12

autovaloreA, autovettoreA = metodoPotenze(A,z,tol)
autovaloreB, autovettoreB = metodoPotenze(B,z3,tol)
autovaloreC, autovettoreC = metodoPotenze(C,z3,tol)

autovaloreAN, autovettoreAN = metodoPotenzeNorm(A,z,tol)
autovaloreBN, autovettoreBN = metodoPotenzeNorm(B,z3,tol)
autovaloreCN, autovettoreCN = metodoPotenzeNorm(C,z3,tol)

print 'Autovalore dominante  |  Autovettore associato'
print("-----------------------------------------------------------------------------")
print("Metodo delle potenze")
print("-----------------------------------------------------------------------------")
print autovaloreA, autovettoreA, 'Matrice A'
print autovaloreB, autovettoreB, 'Matrice B'
print autovaloreC, autovettoreC, 'Matrice C'


print("-----------------------------------------------------------------------------")
print("Metodo delle potenze normalizzato")
print("-----------------------------------------------------------------------------")
print autovaloreAN, autovettoreAN, 'Matrice A'
print autovaloreBN, autovettoreBN, 'Matrice B'
print autovaloreCN, autovettoreCN, 'Matrice C'
