import scipy as sp
import scipy.interpolate as interpol
import matplotlib.pyplot as plt
import os
from Documento import Edificio
l="/home/joaquin/Documents/Entrega 4/Sismos/" #PONER UBICACION CARPETA CONTENEDORA DE INFO SISMOS
lis = os.listdir(l)
lista=[]
for i in lis:
    ll=l+i
    datos = sp.loadtxt(ll)
    maximoa = max(abs(datos))/9.8
    if 0.35 <= maximoa and maximoa<= 0.55:
        lista.append(ll)
    else:
        os.remove(ll)
for j in lista:
    file = open(j, "r") 
    listaa = file.readlines()
    tasa_muestreo = int(listaa[1][20:25])
    nmuestras = int(listaa[3][28:])
    g=9.8
    To=0.
    D= nmuestras/tasa_muestreo
    Ia = (np.pi/(2*g))(sp.integrate((datos)**2, (x, To, To+D)))
    metadatos={"Fecha" :listaa[0][20:30] , "Hora" : listaa[0][31:39], "Estacion_Lat" :listaa[4][11:19],"Estacion_lon" :listaa[4][11:19] , "Componente" :listaa[3][29:32] , "Estacion_Nombre" : listaa[3][12:16], "Epi_Lat":listaa[6][15:20], "Epi_Lon": listaa[6][36:41], "Epi_Profundidad": listaa[7][15:16],"M":listaa[7][28:30]}
    datos=sp.loadtxt(j)
    Vt=np.arange(0, 1, D)
    Nt = datos.size
    dt = 0.01
    t = sp.arange(0,Nt*dt,dt)
    a = interpol.interp1d(t,datos,fill_value=0,kind="cubic")
    sp.savez("acc.npz", a=a, t=t, dt=dt)
    plt.plot(t,datos,"-",markersize=1)
    Edificio(max(abs(datos)))
