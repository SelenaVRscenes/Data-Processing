
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.pyplot as plt 
import matplotlib.animation as animation
from IPython.display import HTML
import numpy as np


""""""
import numpy as np
import matplotlib.pyplot as plt

class Brownian():
    """
    A Brownian motion class constructor
    """
    def __init__(self,x0=0):
        """
        Init class
        """
        assert (type(x0)==float or type(x0)==int or x0 is None), "Expect a float or None for the initial value"
        
        self.x0 = float(x0)
    
    def gen_random_walk(self,n_step=100):
        """
        Generate motion by random walk
        
        Arguments:
            n_step: Number of steps
            
        Returns:
            A NumPy array with `n_steps` points
        """
        # Warning about the small number of steps
        if n_step < 30:
            print("WARNING! The number of steps is small. It may not generate a good stochastic process sequence!")
        
        w = np.ones(n_step)*self.x0
        
        for i in range(1,n_step):
            # Sampling from the Normal distribution with probability 1/2
            yi = np.random.choice([1,-1])
            # Weiner process
            w[i] = w[i-1]+(yi/np.sqrt(n_step))
        
        return w
    
    def gen_normal(self,n_step=100):
        """
        Generate motion by drawing from the Normal distribution
        
        Arguments:
            n_step: Number of steps
            
        Returns:
            A NumPy array with `n_steps` points
        """
        if n_step < 30:
            print("WARNING! The number of steps is small. It may not generate a good stochastic process sequence!")
        
        print("n_step = ", n_step)
        w = np.ones(n_step)*self.x0
        
        for i in range(1,n_step):
            # Sampling from the Normal distribution
            yi = np.random.normal()
            # Weiner process
            w[i] = w[i-1]+(yi/np.sqrt(n_step))
        
        return w
    
    def random_motion_gen(
                    self,
                    s0=10, 
                    mu=0.5099,#0.2312 ,#0.225, #0.22-0.24 ;
                    sigma=0.99,#0.68,
                    deltaT=1,
                    dt=0.001
                    ):
        
       
        
        """
        
        
        default:
            (
                            self,
                            s0=100,
                            mu=0.2,
                            sigma=0.68,
                            deltaT=52,
                            dt=0.1
                            ):
        
        
        
        Models a stock price S(t) using the Weiner process W(t) as
        `S(t) = S(0).exp{(mu-(sigma^2/2).t)+sigma.W(t)}`
        
        Arguments:
            s0: Iniital stock price, default 100
            mu: 'Drift' of the stock (upwards or downwards), default 1
            sigma: 'Volatility' of the stock, default 1
            deltaT: The time period for which the future prices are computed, default 52 (as in 52 weeks)
            dt (optional): The granularity of the time-period, default 0.1
        
        Returns:
            s: A NumPy array with the simulated stock prices over the time-period deltaT
        """
        
        
        
        print("s0 = ", s0)
        
        n_step = int(deltaT/dt)
        time_vector = np.linspace(0,deltaT,num=n_step)
        # Stock variation
        stock_var = (mu-(sigma**2/2))*time_vector
        # Forcefully set the initial value to zero for the stock price simulation
        self.x0=0
        # Weiner process (calls the `gen_normal` method)
        weiner_process = sigma*self.gen_normal(n_step)
        # Add two time series, take exponent, and multiply by the initial stock price
        s = s0*(np.exp(stock_var+weiner_process))
        
        return s, mu, sigma


nframes = 100

""""""

b1 = Brownian()
b2 = Brownian()

#xb = b1.gen_normal(100)
#yb = b2.gen_normal(100)


xb, mu, sigma = b1.random_motion_gen(1)
yb, mu, sigma = b2.random_motion_gen(1)


xmax,xmin,ymax,ymin = xb.max(),xb.min(),yb.max(),yb.min()
scale_factor_x = xmax - xmin
scale_factor_y = ymax - ymin

xb_out  = (xb - xmin) / scale_factor_x
yb_out  = (yb - ymin) / scale_factor_y

np.savetxt('entropyMotion mu = '+str(mu)+', sigma = '+str(sigma)+'.txt',xb_out, delimiter=',', fmt='%f')

plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)

plt.plot(xb,yb,c='b')

plt.savefig('entropyMotion mu = '+str(mu)+', sigma = '+str(sigma)+'.png')


"""""

# 3D visualization 

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

def update(num,x,y,z,u,plot):
    plot[0].remove()
    x_= -20 + xb[num]/4 + x * np.cos(u[num]) - y * np.sin(u[num])
    y_= -20 + yb[num]/4 + x * np.sin(u[num]) + y * np.cos(u[num])
    plot[0] = ax.plot_surface(x_, y_, z, cmap="magma")

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

x = 2 * np.outer(np.cos(u), np.sin(v))
y = 2 * np.outer(np.sin(u), np.sin(v))
#z = 2 * np.outer(np.ones(np.size(u)), np.cos(v))
z = 2 * np.outer(np.ones(np.size(u)), u)

ax.set_xlim(-18.1,18.1)
ax.set_ylim(-18.1,18.1)
ax.set_zlim(-18.1,18.1)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
"""""



# For Anim genegation
#plot= [ax.plot_surface(x, y, z)]
#ani = animation.FuncAnimation(fig, update, nframes, fargs=(x,y,z,u,plot), interval=100)
#ani.save('rotatevballZ.mp4', writer="ffmpeg",dpi=100)

#HTML(ani.to_html5_video())

#anim = animation.FuncAnimation(fig, update, nframes, fargs=(x,y,z,u,plot), interval=100)
#anim.save("circle.gif")
        

cperiodNum = 5 # number of periods for periodic signal
Nsamles = 1000 #Number of data points 

# Generate linear motion

xb_out = []
for p in range(0,cperiodNum):

    r1 = np.arange(0, 1, 2*cperiodNum/Nsamles) #half period
    r2 = np.arange(1, 0, -2*cperiodNum/Nsamles)
    r = np.concatenate((r1, r2))
    xb_out = np.concatenate((xb_out, r))
np.savetxt('linMotion_N='+str(Nsamles)+'.txt',xb_out, delimiter=',', fmt='%f')


# Generate sin motion

xb_out = []
x = np.arange(0, 2*cperiodNum* np.pi, 2* cperiodNum* np.pi/Nsamles) #half period
xb_out = np.sin(x)
np.savetxt('sinMotion_N='+str(Nsamles)+'.txt',xb_out, delimiter=',', fmt='%f')


# Normal distribution
xb = b1.gen_normal(Nsamles)
yb = b2.gen_normal(Nsamles)


xmax,xmin,ymax,ymin = xb.max(),xb.min(),yb.max(),yb.min()
scale_factor_x = xmax - xmin
scale_factor_y = ymax - ymin

xb_out  = (xb - xmin) / scale_factor_x
yb_out  = (yb - ymin) / scale_factor_y

np.savetxt('NormalDistMotion.txt',xb_out, delimiter=',', fmt='%f')

plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)

plt.plot(xb,yb,c='b')

plt.savefig('NormalDistMotion.png')








