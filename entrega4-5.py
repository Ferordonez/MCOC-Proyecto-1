
import scipy as sp
import scipy.interpolate as interpol
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.integrate import solve_ivp
l="/home/joaquin/Documents/Entrega 4/Sismos2/" 
lis=os.listdir(l)
lista=[]
lis.sort()
print (len(lis))
for i in lis:
    ll=l+i
    datos=sp.loadtxt(ll)
    maximoa=max(abs(datos))/9.8
    if 0.2<=maximoa and maximoa<=0.8:
        lista.append(ll)
    else:
        os.remove(ll)
c=0
for j in lista:
    file=open(j,"r") 
    listaa=file.readlines()
    g=9.806
    datos=sp.loadtxt(j)
    A=datos
    Nt=datos.size
    print (j)
    dt=0.01
    tasa_muestreo=float(listaa[1][-19:-14])
    nmuestras=int(listaa[2][28:])
    DT=nmuestras/tasa_muestreo
    T=sp.arange(0,Nt*DT,DT)
    t=sp.arange(0,Nt*dt,dt)
    a=interpol.interp1d(t,datos,fill_value=0,kind="cubic")
    Ia=sp.zeros(Nt)
    v=sp.zeros(Nt)
    d=sp.zeros(Nt) 
    v[1:]=sp.cumsum(datos[1:] + datos[0:-1])*dt/2
    d[1:]=sp.cumsum(v[1:] + v[0:-1])*dt/2 
    a2=datos**2 
    da2=(a2[0:-1] + a2[1:])*dt/2
    Ia[1:]=sp.cumsum(da2)*sp.pi/(2*g) 
    i_PGA=sp.argmax(abs(datos))
    t_PGA=t[i_PGA]
    Ia_inf=Ia.max()
    PGA=(datos[i_PGA])
    i_PGV=sp.argmax(abs(v))
    t_PGV=t[i_PGV]
    PGV=(v[i_PGV])
    i_PGD=sp.argmax(abs(d))
    t_PGD=t[i_PGD]
    PGD=(d[i_PGD])
    i_05=sp.argmin(abs(Ia-0.05*Ia_inf))
    i_95=sp.argmin(abs(Ia-0.95*Ia_inf))
    t_05=t[i_05]
    Ia_05=Ia[i_05]
    t_95=t[i_95]
    Ia_95=Ia[i_95]
    D_5_95=t_95-t_05
    print ("t_95 = ", t_95)
    print ("Ia_95 = ", Ia_95)
    print ("D_5_95 = ", D_5_95)
    metadatos={"Fecha":listaa[0][20:30],"Hora":listaa[0][31:39],"Estacion_Lat":listaa[4][11:19],"Estacion_lon":listaa[4][11:19],"Componente":listaa[3][29:32],"Estacion_Nombre":listaa[3][12:16],"Epi_Lat":listaa[6][15:20],"Epi_Lon": listaa[6][36:41],"Epi_Profundidad":listaa[7][15:16],"M":listaa[7][28:30],"PGA":PGA,"PGV":PGV,"PGD":PGD,"Duracion":D_5_95} 
    plt.figure().set_size_inches([9,6])
    plt.subplot(3,1,1)
    plt.plot(t,datos/g)
    plt.axvline(t_05,color="k",linestyle="--")
    plt.axvline(t_95,color="k",linestyle="--")
    plt.text(t_PGA,PGA/g,"PGA={0:0.3f}g".format(abs(PGA)/g))
    plt.plot(t_PGA,PGA/g,"ob")
    plt.ylim([-0.6,0.6])
    plt.grid(True)
    plt.ylabel("Acc, $a$ (g)")
    plt.subplot(3,1,2)
    plt.plot(t,v*100)
    plt.axvline(t_05,color="k",linestyle="--")
    plt.axvline(t_95,color="k",linestyle="--")
    plt.text(t_PGV,PGV*100,"PGV={0:0.3f}cm/s".format(abs(PGV)*100))
    plt.plot(t_PGV,PGV*100,"ob")
    plt.ylim([-15,15])
    plt.grid(True)
    plt.ylabel("Vel, $v$ (cm/s)")
    plt.subplot(3,1,3)
    plt.plot(t,d*100)
    plt.axvline(t_05,color="k",linestyle="--")
    plt.axvline(t_95,color="k",linestyle="--")
    plt.text(t_PGD,PGD*100,"PGD={0:0.3f}cm".format(abs(PGD)*100))
    plt.plot(t_PGD,PGD*100,"ob")
    plt.ylim([-15,15])
    plt.grid(True)
    plt.xlabel("Tiempo, $t$ (s)")
    plt.ylabel("Dis, $d$ (cm)")
    plt.subplot(3,1,1)
    plt.title(j[-25:-4]+"    $D_{{5-95}} = {0:5.2f}$s".format(D_5_95))
    plt.tight_layout()
    c+=1
    print (c)
    nombre="grafico"+"0"+str(c)+".png"
    nombre2="Registro_"+"0"+str(c)+".npz"
    plt.savefig(nombre)
    plt.show()
    sp.savez(nombre2,metadatos,A,T)
    
    
    
    
    
    
    
    
    
    
    
    