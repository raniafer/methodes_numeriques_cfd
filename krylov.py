import numpy as np
import matplotlib.pyplot as plt

# donnees 
n=16
alpha=2
gamma=1.
be=1.5

x=np.linspace(0,2*np.pi,n)
hij=x[1]-x[0]

# constuire les matrices
D=np.zeros((n,n))
for i in range(n):
	for j in range(n):
		if i!=j :
			if n%2==0 : D[i,j]=(-1)**(i+j)/(2*np.tan(hij))
			else : D[i,j]=(-1)**(i+j)/(2*np.sin(hij))

D2=np.dot(D,D)
beta=be*np.eye(n,n)

M=D2-beta

f=np.cos(alpha*x)

# calcul de la solution utilisant la methode de richardson

U0=np.zeros(n)

rk=f-np.dot(M,U0)

u = []
m=0
while max(abs(rk))>10**(-3) and m<100 :
	alphak = np.dot(rk,rk)/np.dot(rk,np.dot(M,rk))
	U=U0+alphak*rk	
	u.append(U)		
	m=m+1
	#U0=U
	#r0=rk


plt.plot(x,U);
plt.show()
# construire la solution exacte
