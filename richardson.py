import numpy as np
import matplotlib as plt

# U0: la solution initiale, M: la matrice, f:vecteur
def richardson(U0,M,f,eps):
	rk=f-np.dot(M,U0)
	u = []
	m=0
	while max(abs(rk))>eps and m<100 :
		alphak = np.dot(rk,rk)/np.dot(rk,np.dot(M,rk))
		U=U0+alphak*rk	
		u.append(U)		
		m=m+1

def conjugate_gradient(U0,M,f.eps):
	
