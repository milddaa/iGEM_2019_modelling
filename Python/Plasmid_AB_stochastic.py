#### Plasmid AB stochastic modelling ####

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
K7=1
K8=1
K9=1

KD=1
KD2=1
KT=1

Kt1=1
Kt2=1
Kt3=1
Kt4=1
Kt5=1
Kt6=1

alpha=1
beta=1
gamma=1
A=1 #arabinose concentration

dt=0.01
t=0.00
x=0.5
m=0.5
r=0.5
tau=0.01

# Print initial conditions to screen #

print (t, x, m, r)

# Compute the values of x(t) and m(t) for t<1.0 in a step of 0.01 so in total 100 values #

t_list=list()
x_list=list()
m_list=list()
r_list=list()

while t<10.0:
    den = 1.0 + K3*A + K1*KD2*r*r + K1*K2*KD2*r*r*r*r + K7*KD*x*x + K9*K7*KD*KD*x*x*x*x + K4*KT*m*m*m*m + K4*K5*KT*KT*m*m*m*m*m*m*m*m + K8*K7*KD*x*x*KT*m*m*m*m
    num1=Kt1+Kt3*K1*KD2*r*r
    fx= num1/den -alpha*x
    number1 = np.random.poisson(num1/den*tau,1)
    number2 = np.random.poisson(alpha*x*tau,1)
    fx=fx+number1-number2
    num2= Kt2*K3*A + Kt5*K4*KT*m*m*m*m + Kt6*K7*KD*x*x + Kt4
    fm=num2/den - beta*m
    number3 = np.random.poisson(num2/den*tau,1)
    number4 = np.random.poisson(beta*m*tau,1)
    fm=fm+number3-number4
    num3= Kt2*K3*A
    fr=num3/den - gamma*r
    number5 = np.random.poisson(num3/den*tau,1)
    number6 = np.random.poisson(gamma*r*tau,1)
    fr=fr+number5-number6
    x = x + dt*fx
    m = m + dt*fm
    r = r + dt*fr
    t = t + dt
    t_list.append(t)
    x_list.append(x)
    m_list.append(m)
    r_list.append(r)


# Plot c, cox and tetR protein dynamics over time #

fig = plt.figure()
plt.axis([0,10,0,1])
plt.scatter (t_list,x_list, s=1)
plt.scatter (t_list,m_list, s=1)
plt.scatter (t_list,r_list, s=1)
plt.legend (("C protein", "Cox protein", "TetR protein"))
plt.title ("C, Cox and TetR protein dynamics over time in plasmid AB")
plt.xlabel ("Time")
plt.ylabel ("Protein concentration")
plt.show()