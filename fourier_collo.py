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

# calcul de la solution
U=np.dot((np.linalg.inv(M)),f)
plt.plot(x,U,'o');
plt.show()
# construire la solution exacte
