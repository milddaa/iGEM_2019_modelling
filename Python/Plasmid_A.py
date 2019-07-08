#### Plasmid A deterministic modelling ####

#### 2019-06-18 ####

# Importing packages #

import matplotlib.pyplot as plt

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

# Print initial conditions to screen #

print (t,x, m)

# Compute the values of x(t) and m(t) for t<1.0 in a step of 0.01 so in total 100 values #

t_list=list()
x_list=list()
m_list=list()

while t<10.0:
    den = 1.0 + K1*KD*x*x + K2*KT*m*m*m*m + K3*K2*KT*KT*m*m*m*m*m*m*m*m + K1*K4*KD*KD*x*x*x*x + 2*K2*K5*KD*KT*x*x*m*m*m*m 
    fx= Kt1/den -alpha*x
    num= Kt2 + Kt3*K2*KT*m*m*m*m + Kt4*K1*KD*x*x
    fm=num/den - beta*m
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

# Plot f(x) over x and f(m) over x in order to find the number of steady states #
# The number of times f(x) and f(m) cross horizontal axis reflects the number of steady states and the exact concentration of the steady state #

x=0
m=0
dx=0.001
dm=0.001
x_list=list()
m_list=list()
fx_list=list()
fm_list=list()

while x<5:
    den = 1.0 + K1*KD*x*x + K2*KT*m*m*m*m + K3*K2*KT*KT*m*m*m*m*m*m*m*m + K1*K4*KD*KD*x*x*x*x + 2*K2*K5*KD*KT*x*x*m*m*m*m
    fx= Kt1/den -alpha*x
    x=x+dx
    x_list.append(x)
    fx_list.append(fx)
    num= Kt2 + Kt3*K2*KT*m*m*m*m + Kt4*K1*KD*x*x
    fm=num/den - beta*m
    m=m+dm
    m_list.append(m)
    fm_list.append(fm)

fig = plt.figure()
plt.axis([0,5,-5,5])
plt.scatter (x_list,fx_list,s=1)
plt.scatter (m_list,fm_list,s=1)
plt.legend (("C protein", "Cox protein"))
plt.axhline(y=0, color='r', linestyle=':')
plt.title ("f(x) versus x")
plt.xlabel ("x")
plt.ylabel ("f(x)")
plt.show()