import numpy as np
import matplotlib.pyplot as plt

# donnees 
n=16
alpha=2
gamma=1.
be=1.5

x=np.linspace(0,2*np.pi,n)
t=np.linspace(0.,2*np.pi,n)
f=np.cos(alpha*x)
#___________________________________________________________________________
def fourier_collocation(x, f, be):
	# constuire les matrices
	D=np.zeros((n,n))
	for i in range(n):
		for j in range(n):
			if i!=j :
				hij = 0.5*(x[i]-x[j]);
				if n%2==0 : 
					D[i,j]=(-1)**(i+j)/(2*np.tan(hij))
				else : 
					D[i,j]=(-1)**(i+j)/(2*np.sin(hij))

	D2=np.dot(D,D)
	beta=be*np.eye(n,n)
	M=D2-beta

	# calcul de la solution
	U=np.dot((np.linalg.inv(M)),f)
	return(U)
	
U = fourier_collocation(x, f, be)

plt.plot(x,U,'o')
plt.show()
#___________________________________________________________________________	
def fourier_galerkin(x, f, be):
	# calcul des coefficients de fourier
	fk = np.fft.fft(f)
	k=np.fft.fftfreq(n)*n
	uk=np.zeros(n)
	
	for i in range(n):
		uk=fk/(be-k**2)
		
	U=np.fft.ifft(uk).real
	return(U)

U2 = fourier_galerkin(x, f, be)

plt.plot(x,U2,'o')
plt.show()	
#_________________________________________________________________________________________________________________________
alpha2=2*np.pi
xt=np.cos(np.pi*np.arange(n)/(n-1))
ft=(be-alpha2**2)*np.cos(alpha*xt)

def tchebychev_collocation(xt, ft, be):
	# construction des matrices
	D=np.zeros((n,n))
	for i in range(n):
		for j in range(n):
			if i==0 or i==n: Ci=0
			if j==0 or j==n: Cj=0
			else : Ci=1; Cj=1
			if i!=j:
				D[i, j]=Ci/Cj*(-1)**(i+j)/(2*(xt[i]-xt[j]))
			else :
				D[i, j]=-xt[i]/(2*(1-xt[i]**2))
	D2=np.dot(D,D)	
	M=D2+be*np.eye(n)
	# CL
	M[0,0]=1
	M[-1,-1]=1
	M[n-1,:] = 0
	M[n-1,n-1] = 1.	
	ft[0]=1
	ft[-1]=1	

	U=np.dot(np.linalg.inv(M),ft)
	return(U)

U3=tchebychev_collocation(xt, ft, be)

plt.plot(x,U3,'o')
plt.show()
#___________________________________________________________________________	

def tchebychev_galerkin(xt, ft, be):
	D=np.zeros((n,n))
	for i in range(n):
		if i==0: Ci=2
		else: Ci=1
		for j in range(i+1,n,2):
			D[i, j]=2/Ci*j
			
	D2=np.dot(D,D)	
	M=D2+be*np.eye(n)
	# CL
	M[n-2,:] = 1;
	M[n-1,0::2] = 1;
	M[n-1,1::2] = -1;
	ft[n-2] = 1.
	ft[n-1] = 1.
	
	uk=np.dot(np.linalg.inv(M),ft)
	N = np.size(uk) - 1;
	uk2 = uk.copy();
	uk2[0] *= 2.;
	uk2[N] *= 2.;
	u2k = np.concatenate((uk2,np.flipud(uk2[1:N])));
	u2  = np.real(np.fft.fft(u2k));
	U=0.5*u2[0:N+1]
	return(U)
	
U4=tchebychev_galerkin(xt, ft, be)

plt.plot(x,U4,'o')
plt.show()
# construire la solution exacte
