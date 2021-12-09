import math 
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
from numpy.core.defchararray import array
from numpy.core.numeric import zeros_like

def heatfd(xl,xr,yb,yt,M,N):
    def f(x):
        return np.sin(2*np.pi*x)**2
    def l(t):
        return 0*t 
    def r(t):
        return 0*t
    D=1
    h=(xr-xl)/M
    k=(yt-yb)/N
    m=M-1
    n=N
    sigma=D*k/(h*h)
    a=np.diag([1-2*sigma]*m)+np.diag([sigma]*(m-1),1)+np.diag([sigma]*(m-1),-1)
    lside = l(yb+np.array(range(0,n))*k)
    rside = r(yb+np.array(range(0,n))*k)
    w=np.ones((n,m))
    x1=f(xl+np.array(range(1,m+1))*h)
    w[1]=f(x1)
    for j in range(1,n-1):
        m1 = sigma*np.concatenate(([lside[j]],np.zeros(m-2),[rside[j]]),axis=0)
        wj=w[j].T
        w[j+1]=a.dot(wj) +m1
    w=np.concatenate(([lside],w.T,[rside]),axis=0)
    x=np.array(range(0,m+2))*h
    x=x.reshape(1,m+2)
    t=np.array(range(0,n))*k
    t=t.reshape(1,n)
    return w,x,t
w,x,t=heatfd(0,1,0,1,10,250)
print(w.T[1])
ax = plt.subplot(111, projection='3d')
x, t = np.meshgrid(x, t)
ax.plot_surface(x,t,w.T,rstride=1,  # rstride（row）指定行的跨度
                       cstride=1 ) # cstride(column)指定列的跨度
plt.show()
