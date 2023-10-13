import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit,bisect
from scipy.constants import G

#Earth Constants
M_earth=5.9722*10**24 # kg
R_earth=6371*10**3 # m
Rho_earth = 1.217 #kgm^-3

#Moon constants
M_moon = 0.0123*M_earth
R_moon=0.2727*R_earth
Rho_moon = 0

#Mercury constants
M_mars=0.107*M_earth #kg
R_mars=0.532*R_earth#m
Rho_mars=0.020 #kg/m^3

#Venus constants
M_venus = 0.815*M_earth
R_venus = 0.950*R_earth
Rho_venus = 65 #kgm^-3


x_f=20000

def central_diff(func,x,step):
    return (func(x+step/2)-func(x-step/2))/step

def grav_accl(M,R):
    return -G*M/(R**2)

def rk4(y,t,h,derivs,args):
    #function to implement rk4
    #y = [x,v] current state
    #t = current time
    #h = time step
    #derivs = derivative function that defines the problem
    k1,k2,k3,k4 = np.zeros(2),np.zeros(2),np.zeros(2),np.zeros(2)
    k1 = h*derivs(y,t,args)
    y_halfstep = y + k1/2. #Euler half step using k1
    k2 = h*derivs(y_halfstep,t+h/2,args)
    y_halfstep = y + k2/2. #Euler half step using k2
    k3 = h*derivs(y_halfstep,t+h/2,args)
    k4 = h*derivs(y + k3,t+h,args) #full step using k3
    y_next = y + (k1+2*k2+2*k3+k4)/6.
    return y_next

def Newt_Raph(func,x_0,h,tol,N):
    '''func takes in function to be evaluated
    x_0 is first guess
    h is accuracy for derivative
    tol dertermines how close to the exact solution we want
    N determines iteration limit'''
    while abs(func(x_0)) > tol:
        for i in range(N):
            dx = -(func(x_0))/(central_diff(func,x_0,h))
            x_0 += dx
            print(func(x_0))
        x_0 =x_0 +dx/2
        return x_0

def proj_mot(z,t,args):
    #z=[x_t,vx_t,y_t,vy_t]
    #F=-mg-Fd
    Rho, A, Cd, m, g = args
    zp=np.zeros(4)
    vmag = np.sqrt(z[1]**2 + z[3]**2)
    zp[0]=z[1]
    zp[2]=z[3]
    zp[1]=-0.5*(Rho/m)*A*Cd*z[1]*vmag
    zp[3]=g-0.5*(Rho/m)*A*Cd*z[3]*vmag
    return zp

planets=np.array([[M_earth,R_earth,Rho_earth],[M_moon,R_moon,Rho_moon],[M_mars,R_mars,Rho_mars],[M_venus,R_venus,Rho_venus]])
def motion(pl,v_0,theta):
    Cd=0.04 #Streamlined booger
    d = 1
    A= np.pi*d**2/4
    m = 2250*A*(2*d) #kg
    M=pl[0]
    R=pl[1]
    Rho=pl[2]
    g = grav_accl(M,R)
    args = [Rho, A, Cd, m, g]

    N=100000 #number of steps
    T=200 #Total time
    h=T/(float(N-1)) #step size
    time=np.arange(0,T+h,h)

    #Initial conditions
    y_0=0
    x_0=0

    vx_0=v_0*np.cos(np.radians(theta))
    vy_0=v_0*np.sin(np.radians(theta))

    states=np.zeros((N,4))
    #print(type(x_0),type(vx_0),type(y_0),type(vy_0))
    z=[x_0,vx_0,y_0,vy_0]
    states[0,:]=z

    for j in range(0,N-1):
        states[j+1,:]=rk4(states[j,:],time[j],h,proj_mot,args)
        if states[j+1,2] < 0:
            states = states[:j,:]
            #print(states[-1,:])
            break
    return states

def plot_trajec(planet,v1,v2,th,name,pos):
    y_f=1000
    xmax=0
    pl=planet

    def xfunc(veloc):
        return x_f - max(motion(pl,veloc,th)[:,0])
    
    '''Vel= bisection(xfunc,v1,v2,10)
    States = motion(planet,Vel,th)
    xmax=States[-1,0]
    y_f = States[-1,2]
    print(name, '1 :',xmax,y_f,Vel)'''
    
    Vel2=Newt_Raph(xfunc,v1,1,1e-4,15)
    States2 = motion(planet,Vel2,th)
    xmax2=States2[-1,0]
    y_f2 = States2[-1,2]
    print(name,'2 :',xmax2,y_f2,Vel2)

    ax=fig.add_subplot(2,2,pos)
    x,y=States2[:,0],States2[:,2]
    ax.plot(x,y,label=name)
    ax.legend()
    ax.set_xlim(0,20000)
    ax.set_ylim(0,20000)
    
fig=plt.figure()
'''plot_trajec(planets[0],450,500,35,'Earth',1)
plot_trajec(planets[1],160,210,35,'Moon',2)'''
plot_trajec(planets[2],255,305,35,'Mars',3)
#plot_trajec(planets[3],53350,53550,35,'Venus',4)
