import numpy as np
import matplotlib.pyplot as plt

#P1, P0 in same dimension 
def diff(P1,P0,Order):
  if np.shape(P1) == np.shape(P0):
    if Order == 1:
      for i in range(np.shape(P0)-1):
        dP0 = P0[i+1]-P0[i]
        dP1 = P1[i+1]-P0[i]
        if dP0>1:
          print('differential error')
          return 0
        else:
          return dP1/dP0
    else:
      while Order > 1:
        for i in range(np.shape(P0)-1):
          dP0 = P0[i+1]-P0[i]
          dP1 = P1[i+1]-P0[i]
          if dP0>1:
            print('differential error')
            return 0
          else:       
            P0 = dP0
            P1 = dP1
            Order -= 1 
            diff(P1,P0,Order)
  else:
    print('dimension error')
    return 0
#in the Hamiltonian of our OQS master equation in which the dynamics of the proton transfer are described following the Caldeira and Leggett quantum Brownian motion mode（WM–CL）
#parameters at t0
#T = 298.15 K
hbar = 1.0546e-34  # Planck's constant over 2*pi
e = 1.6022e-19     # Electron charge
L = 5.2918e-11     # Bohr radius
N = 1000

kb = 8.617333262*10**(-5) #eV/K
T =  np.linspace(100, 1000,20)
h = hbar*2*np.pi
#W0 = ()
#Wt = p/m*Wq+Vq*Wp-

#for biologically relevant bath temperatures(T ≈ 300 K), model requires a low-temperature condition, whereby the thermal energy is replaced with the zero-point energy
wp = 0.00367 #unit a.u. the frequency of the lowest energy eigenstate in the double well
gamma = 39 #m(1)
temp = np.array(hbar*wp/2/kb/np.array(T,float)*10**34,float)
print(temp)
print(np.exp(temp))
coth1 = (np.exp(temp)+np.exp(-temp))/(np.exp(temp)-np.exp(-temp))
print(coth1)
Ezero = hbar*wp*10**34*coth1/2
print(Ezero)
print(kb*T)
a0 = plt.figure()
plt.plot(T,kb*T,color = 'red')
plt.plot(T,Ezero,'o',color = 'blue',)
plt.show()
# Note that these corsrelation functions are explicit for baths consisting of simple harmonic oscillators within this model
# power spectral density function Jrr, can be thought of as describing the fundamental properties of the heat bath and these are typically altered to fit experimental data.
a = 4*2**(1/2)
c = np.array(h/kb/T,float)
Css = a/4*(1/wp**2*(1/(np.exp(gamma)*np.exp(-gamma*1j)-1)-1/(np.exp(-gamma)*np.exp(-gamma*1j)-1)-c**2))
