

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import LinearLocator
from matplotlib.colors import Normalize
import numpy as np 

import matplotlib.animation as animation
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.linalg import block_diag


np.set_printoptions(threshold=40)

N=10

h=1/4/N
dim =6*N*(N-1)*2
print(dim)
# i,j 

def alpha(x,y) :
    return 1
    if(x==0):
        return y-0.5
    elif(x==1):
        return -y+0.5
    elif(y==0):
        return x-0.5
    else:
        return -x+0.5

A =  [ [-1,4,1] , [4,-16,4] ,[1,4,-1] ] #[ [-1,4,1] , [4,-16,4] ,[1,4,1] ] 

def boundary(i,j):
    if(i==0 or i==4*N):
        if(j==0 or j==4*N):
            k , l ,b = [1,1,4*N-1,4*N-1],[1,4*N-1,1,4*N-1] ,0
        if(j>0 and j<4*N):
            k , l ,b = [1,4*N-1],[j,j] ,0
    if(i==N):
        if(j>=N and j<=3*N):
            k , l ,b = N-1,j ,alpha(i*h,j*h)

    if(i==3*N):
        if(j>=N and j<=3*N):
            k , l ,b = 3*N+1,j ,-alpha(i*h,j*h)
    
    
    if(j==0 or j==4*N):
        if(i>0 and i<4*N):
            k , l ,b = [i,i],[1,4*N-1] ,0

    if(j==N):
        if(i>N and i<3*N):
            k , l ,b = i-1,j-1 ,alpha(i*h,j*h)
    if(j==3*N):
        if(i>N and i<3*N):
            k , l ,b = i+1,j+1 ,-alpha(i*h,j*h)


    return k,l,b


def isBoundary(i,j):
    if(i==0 or i==4*N or j==0 or j==4*N  ):
        return True
    if((i>=N and i<=3*N) and (j>=N and j<=3*N) ):
        return True
    #if((j>=N and j<=3*N) and (i>=N and i<=3*N) ):
    #    return True
    return False

notBoundary = lambda i,j : not isBoundary(i, j)

def newLine():
    line = []
    for i in range(1,4*N):
        subline = []
        for j in range(1,4*N):
            if(not isBoundary(i, j)):
                subline = subline+[0]
        if(len(subline)>0):
            line= line+[subline]
            #print(len(subline))
    return line
    
def rindex(k,l):
    if(k<N or k>3*N):
        pass
    elif(l<N ):
        pass
    else:
        l=l-2*N-1
    return k-1 , l-1

def pde(i,j):
    line = newLine()
    b=0
    for kindex ,k in enumerate([i-1,i,i+1]):
        for lindex ,l in enumerate([j-1,j,j+1]):
            bkl=0
            if(not isBoundary(k, l)):
                newk,newl = k,l
            else:
                newk,newl , bkl = boundary(k, l)
        
            if(type(newk)==list):
                for index in range(len(newk)):
                    #print(newk,newl)
                    newk[index],newl[index]=rindex(newk[index],newl[index])
                    line[newk[index]][newl[index]] = line[newk[index]][newl[index]]+ A[kindex][lindex]/len(newk)
                    #pass
            else:
                newk,newl=rindex(newk,newl)
                line[newk][newl] = line[newk][newl]+ A[kindex][lindex]
            b=b+bkl
    return line,b

def mySystem():
    Ah, bh =[],[]
    for i in range(4*N):
        for j in range(4*N+1):
            if(not isBoundary(i, j)):
                line , b= pde(i, j)
                #print(sum(np.concatenate(line)))
                Ah=Ah+[np.concatenate(line)]
                bh=bh+[b]
    return np.array(Ah) ,np.array(bh)*h

myAh,mybh = mySystem()

print(np.shape(myAh), len(mybh))

notinvertible  =np.linalg.det(myAh) ==0
if(notinvertible):
    print("singular matrix :(")
else:
    print("invertible matrix :)")
    zeta=np.linalg.solve(myAh,np.transpose(mybh))


'''
print(np.sum(myAh,axis=1))

sumabs =  lambda l : sum([np.abs(v) for v in l])
#print([sumabs(myAh[:,i]) for i in range(dim)]  )
'''


subZeta1 = np.reshape(zeta[0:(N-1)*(4*N-1)],(N-1,4*N-1))

subZeta2 =np.reshape(zeta[(N-1)*(4*N-1):(N-1)*(4*N-1)+2*(N-1)*(2*N+1)],(2*N+1,2*N-2)) 

subZeta21 = subZeta2[:,0:N-1]
subZeta22 = subZeta2[:,N-1:2*N+1]

subZeta3 = np.reshape(zeta[(N-1)*(4*N-1)+2*(N-1)*(2*N+1):dim],(N-1,4*N-1))

zeroZeta = np.ones((2*N+1,2*N+1))*(max(zeta)+min(zeta))/2

#print(np.shape(subZeta21), np.shape(zeroZeta))

subZeta2 =np.concatenate ([ subZeta21,zeroZeta ,subZeta22 ],axis=1)

newZeta =np.concatenate ([ subZeta1,subZeta2 ,subZeta3 ],axis=0)


#boundaries 
# from outside
xBound = (newZeta[0,:]/2+newZeta[-1,:]/2).reshape((1,4*N-1))
newZeta = np.concatenate ([ xBound,newZeta ,xBound ],axis=0)

yBound = (newZeta[:,0]/2+newZeta[:,-1]/2 ).reshape((4*N+1,1))
newZeta = np.concatenate ([ yBound,newZeta ,yBound ],axis=1)

# inner 

#print(">>>>>>>>>>",newZeta[N:3*N,N],newZeta[N:3*N,3*N] , newZeta[N,N:3*N],newZeta[3*N,N:3*N])
newZeta[N:3*N,N]=newZeta[N-1:3*N-1,N-1]+h
newZeta[N:3*N,3*N]=newZeta[N+1:3*N+1,3*N+1]-h
newZeta[N,N:3*N]=newZeta[N-1,N:3*N]+h
newZeta[3*N,N:3*N]=newZeta[3*N+1,N:3*N]-h

#ploting

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
x = np.linspace(0,1 , 4*N+1)
y = np.linspace(0,1 , 4*N+1)

X,Y = np.meshgrid(x,y)
ax.plot_surface(X,Y,newZeta)

#saving
plt.rcParams['text.usetex'] = True
plt.title(r"Plot of $\gamma$ , $N=10$")

plt.savefig("gamma.png")

plt.show()
plt.close()
