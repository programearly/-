import numpy as np
import math
def halton(p,n):
    b=np.zeros((math.ceil(math.log(n)/math.log(p)),1))
    u=[0]*n
    for j in range(n):
        i=0
        b[0]=b[0]+1        
        while b[i]>p-1+np.spacing(1) :
            b[i]=0
            i+=1
            b[i]=b[i]+1
        u[j]=0
        for k in range(len(b)):
            u[j]=u[j]+b[k]*pow(p,-k-1)
    return u

print(halton(2,100))
