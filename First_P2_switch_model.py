#### Phage P2 switch deterministic modelling ####

#### 2019-06-13 ####

# Importing packages #

import matplotlib.pyplot as plt

# EXPLICIT EULER SCHEME #

# Setting up initial conditions and parameters #

u=1
sigma=2
gamma=3
teta=2

dt=0.01
x=0.8
t=0.00

# Print initial conditions to screen #

print (t,x)

# Compute the values of x(t) for t<1.0 so in total 100 values #

t_list=list()
x_list=list()

while t<100.0:
    num= 1.0 + teta*u*x*x 
    den = 1.0 + u*x*x + sigma*u*u*x*x*x*x
    fx= num/den -gamma*x
    x = x + dt*fx
    t = t + dt
    t_list.append(t)
    x_list.append(x)

# Plot c protein dynamics over time #

fig = plt.figure()
plt.axis([0,50,0,1])
plt.title("C protein dynamics over time")
plt.xlabel("Time, units")
plt.ylabel("Protein concentration")
plt.scatter (t_list,x_list,s=1)
plt.show()

# Plot f(x) over x in order to find the number of steady states #
# The number of times f(x) crosses horizontal axis reflects the number of steady states #

x=0
dx=0.01
x_list=list()
fx_list=list()

while x<5:
    num= 1.0 + teta*u*x*x
    den = 1.0 + u*x*x + sigma*u*u*x*x*x*x
    fx= num/den -gamma*x
    x=x+dx
    x_list.append(x)
    fx_list.append(fx)

fig = plt.figure()
plt.axis([0,5,-5,5])
plt.scatter (x_list,fx_list,s=1)
plt.axhline(y=0, color='r', linestyle=':')
plt.title("f(x) versus x")
plt.xlabel("x")
plt.ylabel("f(x)")
plt.show()