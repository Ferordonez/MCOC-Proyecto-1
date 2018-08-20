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
    fecha = listaa[0][20:30]
    hora = listaa[0][31:39]
    latitud =listaa[4][11:19]
    longitud =listaa[4][29:37]
    componente = listaa[3][29:32]
    estacion= listaa[3][12:16]
    metadatos={"fecha" : fecha , "Hora" : hora, "Latitud" : latitud,"longitud" : longitud , "componente" : componente , "estacion" : estacion}
    datos=sp.loadtxt(j)
    Nt = datos.size
    dt = 0.01
    t = sp.arange(0,Nt*dt,dt)
    a = interpol.interp1d(t,datos,fill_value=0,kind="cubic")
    sp.savez("acc.npz", a=a, t=t, dt=dt)
    plt.plot(t,datos,"-",markersize=1)
    plt.show()
    Edificio(max(abs(datos)))