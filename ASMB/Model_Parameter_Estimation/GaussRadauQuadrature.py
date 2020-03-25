# Created by Kensuke Suzuki in Mar 2020

from math import *
import sys
import numpy as np

def lglnodes(N1):
        
    N=N1-1
    x=np.zeros(N1)
    xold=np.zeros(N1)
    for i in range(0,N1):
        x[i]=-cos(2*pi*i/(2*N+1))
        xold[i] = 2
    
    #% The Legendre Vandermonde Matrix
    #P=np.ones((N1,N1+1);
    
    # Compute P_(N) using the recursion relation
    # Compute its first and second derivatives and 
    # update x using the Newton-Raphson method.
       
    # Free abscissae
    free=[i for i in range(1,N1)]
    
    n=[i for i in range(0,N1+1)]
    P=np.zeros((N1,N1+1))
       
    while max(abs(x-xold)) > sys.float_info.epsilon:
    
        xold=x.copy()
        
        for i in n:
            P[0,i]=(-1)**i      
    
        P[free,0]=1
        P[free,1]=x[free]
        
        for k in range(1,N1):
            P[free,k+1]=( (2*(k+1)-1)*x[free]*P[free,k]-(k+1-1)*P[free,k-1] )/(k+1)
    
        x[free]=xold[free]-((1-xold[free])/N1)*(P[free,N]+P[free,N1])/(P[free,N]-P[free,N1])
    
    # The Legendre-Gauss-Radau Vandermonde
    P=P[0:N1,0:N1]
    
    # Compute the weights
    w=np.zeros(N1)
    w[0]=2/(N1**2);
    w[free]=(1-x[free])/(N1*P[free,N1-1])**2
    
    return w[::-1]
