#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 16:37:58 2018

@author: joaquin
"""
from numpy.linalg import inv
import numpy as np 
import scipy as sp 
from numpy import *
from matplotlib import pyplot as plt
def Edificio(d1):
    Cap = array([0, 1050, 800, 0, 750, 500, 0, 300, 0, 300, 300, 0, 250, 0, 250, 0, 250, 150, 0, 0])*0
    a = np.sum(Cap)
    b = len(Cap)
    #print a,b
    largo = [12.*3.5, 12.*3.5, 12.*3.5, 12.*3.5, 8.*3.5, 8.*3.5, 8.*3.5, 8.*3.5, 4.*3.5, 4.*3.5, 4.*3.5, 4.*3.5, 4.*3.5, 4.*3.5, 4.*3.5, 4.*3.5, 4.*3.5, 4.*3.5, 4.*3.5, 4.*3.5] 
    masa = np.zeros(20)
    
    l = 0
    while l < len(largo):
        masa[l] = 6.2*largo[l]
        l += 1 
    
    #print masa
    
    M = sp.transpose(masa)*sp.identity(20)  
    M = array(M)
    M *= 1000 
     
     
    
    #print (M)
    
    e = 23.5*1000 #GPa
    col = array([e*(600**4), e*(700**4), e*(800**4), e*(900**4), e*(1000**4)])  #12EI
    
    #print (col) 
    
    piso1 = array([6, 2, 0, 1, 4])
    piso2 = array([4, 0, 1, 0, 4]) 
    piso3 = array([0, 0, 1, 4, 0]) 
    piso4 = array([0, 1, 4, 0, 0])
    piso5 = array([5, 0, 0, 0, 0])
    
    k=[]
    for i in range(4):
        k.append(sp.sum(col*piso1))
    for i in range(4):
        k.append(sp.sum(col*piso2))
    for i in range(4): 
        k.append(sp.sum(col*piso3))
    for i in range(4): 
        k.append(sp.sum(col*piso4)) 
    for i in range(4): 
        k.append(sp.sum(col*piso5))
        
    #Dividimos por los largos al cubo
    k = array(k)
    k[0] /= 4000**3    
    k[1:] /= 2800**3
    
    
    #print (k) #[KN/m]
    
    k_total = np.zeros((20,20)) 
    
    for i in range(19): 
        k_total[i][i] = k[i] + k[i+1] 
        k_total[i][i+1] = -k[i+1] 
        k_total[i+1][i] = -k[i+1]
    
    k_total[19][19] = k[19]
    k_total = np.matrix(k_total)
    
    #print (k_total) #matriz de rigidez 
    
    
    wn = sp.sqrt(k/masa) 
    fn = wn/2/sp.pi
    Tn= 1/fn 
    dt = 0.001
    tmax = 100
    vr = 0.0001 
    d0 = d1
    v0 = 0.
    
    Id = np.matrix(sp.identity(20))
    zero = np.matrix(np.zeros((20,20)))
    M_I = inv(M)
    
    c = np.zeros((20,20)) 
    f1 = 0.2
    f2 = 2. 
    x1 = 0.025
    
    a0 = (4*sp.pi*f1*f2*x1)*(f1-f2)/(f1**2-f2**2)
    a1 = x1*(f1-f2)/(sp.pi*(f1**2-f2**2))
    
    c = a0*M + a1*k_total  
    c = np.matrix(c)
    #print (c) 
    
    
    
    #print (masa_1)
    A =  np.block([[zero,Id],[-M_I*k_total, -M_I*c]])
    #print (A)
    
    
    
    def fun(t,z):
         if t > fun.tnextreport:
            #print ("  {} at t = {}".format(fun.solver, fun.tnextreport))
            fun.tnextreport += 1
         Fr = sp.zeros(40) 
         Fr[19] = -Cap[0]*sp.tanh(z[20]/vr)
         for i in range(20): 
            Fr[i+20] = -Cap[i]*M_I[i,i]*sp.tanh(z[i+20]-z[i+19]/vr)
         
         return sp.matmul(A,z) + Fr
    
    
    t = sp.arange(0, tmax, dt)
    Nt = len(t)
    
    z_euler = sp.zeros((40,Nt+1))
    z_RK45 = sp.zeros((40,Nt+1))
    
    z0 = sp.zeros(40) 
    i = 1
    while i <= 19: 
        z0[i] = 0.5*i
        i += 1 
    z0[20:] = v0
    #print (z0)
    
    z_euler[:,0] = z0
    z_RK45[:,0] = z0 
    
    print ('Integrando con Euler')
    fun.tnextreport = 0
    fun.solver = 'Euler'
    
    i = 1
    ti = dt 
    while (ti < tmax): 
        z_euler[:,i] = dt*fun(ti, z_euler[:,i-1]) + z_euler[:,i-1]
        ti += dt 
        i += 1
    
    #print ('Integrando con RK45') 
    #fun.tnextreport = 0
    #fun.solver = 'RK45'
    #solucion_rk45 = solve_ivp(fun, [0., tmax], z0, methor = 'RK45', t_eval=t,vectorized=False ) 
    #z_RK45[:,1:] = solucion_rk45.y
    #
    #plt.figure()
     
    
    for z, lab in zip([z_euler, z_RK45], ["Euler", "RK45"]):
        u = z[0,:]
        v = z[1,:]
        #Extraer desplazamientos y velocidades
        u = z[0,:-1]
        v = z[1,:-1]
        plt.subplot(2,1,1)
        plt.plot(t, u, label=lab)
        #plt.ylim([-1.5*d0, 1.5*d0])
        plt.xlim([0, tmax])
        plt.ylabel("Despazamiento, $u = z_1$ (m)")
        plt.grid(True)
        vmax = max(abs(v))
        plt.subplot(2,1,2)
        plt.plot(t, v)
        plt.ylabel("Velocidad, $\dot{u} = z_2$ (m/s)")
        plt.xlabel("Tiempo, $t$ (s)")
        #plt.ylim([-1.5*vmax, 1.5*vmax])
        plt.xlim([0, tmax])
        plt.grid(True) 
    plt.subplot(2,1,1)
    plt.legend()
    plt.suptitle("Solucion por metodo de Euler")
    plt.show()
    np.savez('MatrixM.npz',M)
    np.savez('MatrixK.npz',k_total)
    np.savez('MatrixC.npz', c)
