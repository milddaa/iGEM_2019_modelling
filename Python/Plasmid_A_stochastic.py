#### Plasmid A stochastic modelling ####

#### 2019-06-27 ####

# Importing packages #

import matplotlib.pyplot as plt
import numpy as np

# EXPLICIT EULER SCHEME #

# Setting up initial conditions and parameters #

K1=1
K2=1
K3=1
K4=1
K5=1
K6=1
KD=1
KT=1
Kt1=10
Kt2=1
Kt3=1
Kt4=1
alpha=1
beta=1

dt=0.01
x=0.5
t=0.00
m=0.5
tau=0.01

# Print initial conditions to screen #

print (t,x, m)

# Compute the values of x(t) and m(t) for t<1.0 in a step of 0.01 so in total 100 values #
# Updated by the tau-leap method #

t_list=list()
x_list=list()
m_list=list()

while t<10.0:
    den = 1.0 + K1*KD*x*x + K2*KT*m*m*m*m + K3*K2*KT*KT*m*m*m*m*m*m*m*m + K1*K4*KD*KD*x*x*x*x + 2*K2*K5*KD*KT*x*x*m*m*m*m 
    fx= Kt1/den -alpha*x
    number1 = np.random.poisson(Kt1/den*tau,1)
    number2 = np.random.poisson(alpha*x*tau,1)
    fx= fx + number1 - number2
    num= Kt2 + Kt3*K2*KT*m*m*m*m + Kt4*K1*KD*x*x
    fm=num/den - beta*m
    number3 = np.random.poisson(num/den*tau,1)
    number4 = np.random.poisson(beta*m*tau,1)
    fm=fm + number3 - number4
    x = x + dt*fx
    m = m + dt*fm
    t = t + dt
    t_list.append(t)
    x_list.append(x)
    m_list.append(m)


# Plot c and cox protein dynamics over time #

fig = plt.figure()
plt.axis([0,10,0,2])
plt.scatter (t_list,x_list, s=1)
plt.scatter (t_list,m_list, s=1)
plt.legend (("C protein", "Cox protein"))
plt.title ("C and Cox protein dynamics over time in plasmid A")
plt.xlabel ("Time, units")
plt.ylabel ("Protein concentration")
plt.show()