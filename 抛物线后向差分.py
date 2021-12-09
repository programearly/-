import math 
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
from numpy.core.defchararray import array
from numpy.core.numeric import zeros_like


def heatbdn(xl,xr,yb,yt,M,N):
    def f(x):
        return np.sin(2*np.pi*x)**2
    D=1
    h=(xr-xl)/M
    k=(yt-yb)/N 
    m=M+1
    n=N 
    sigma = D*k/(h*h)
    a=np.diag([1-2*sigma]*m)+np.diag([sigma]*(m-1),1)+np.diag([sigma]*(m-1),-1)
    a[:,1]=np.concatenate(([-3,4,-1],np.zeros((m-3))),axis=0)
    a[:,m-1]=np.concatenate((np.zeros((m-3)),[-1,4,-3]),axis=0)
    w=np.ones((n,m))
    w[1]=f(xl+np.array(range(0,M+1))*h)
    for j in range(0,n-1):
        b=w[j]
        b[0]=0
        b[m-1]=0
        w[j+1]=np.linalg.solve(a,b)
    x=np.array(range(0,m))*h
    x=x.reshape(1,m)
    t=np.array(range(0,n))*k
    t=t.reshape(1,n)
    return w,x,t
w,x,t=heatbdn(0,1,0,1,20,20)
print(w.T[1])
ax = plt.subplot(111, projection='3d')
x, t = np.meshgrid(x, t)
ax.plot_surface(x,t,w,rstride=1,  # rstride（row）指定行的跨度
                       cstride=1 ) # cstride(column)指定列的跨度
plt.show()
