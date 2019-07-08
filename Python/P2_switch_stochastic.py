#### Phage P2 switch stochastic modelling ####

#### 2019-06-27 ####

# Importing packages #

import matplotlib.pyplot as plt
import numpy as np

# Setting up initial conditions and parameters #

u=1
sigma=2
gamma=3
teta=2

dt=0.01
x=0.8
t=0.00
tau=0.01


# Compute the values of x(t) for t<1.0 so in total 100 values #
# Updated by the tau-leap method #

t_list=list()
x_list=list()

while t<100.0:
    num= 1.0 + teta*u*x*x 
    den = 1.0 + u*x*x + sigma*u*u*x*x*x*x
    fx= num/den -gamma*x
    number1 = np.random.poisson(num/den*tau,1)
    number2 = np.random.poisson(gamma*x*tau,1)
    fx=fx + number1 - number2
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
plt.scatter (t_list,x_list, s=1)
plt.show()
