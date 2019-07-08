#### Plasmid B deterministic modelling ####

#### 2019-06-18 ####

# Importing packages #

import matplotlib.pyplot as plt

# EXPLICIT EULER SCHEME #

# Setting up initial conditions and parameters #

K1=1
K2=1
K3=1
KD=1
Kt1=10
Kt2=1
Kt3=1
alpha=1
beta=1
gamma=5
A=1 #arabinose concentration

dt=0.01
x=0.5
t=0.00
m=0.5
r=0.5

# Print initial conditions to screen #

print (t,x, m)

# Compute the values of x(t) and m(t) for t<1.0 in a step of 0.01 so in total 100 values #

t_list=list()
x_list=list()
m_list=list()
r_list=list()

while t<10.0:
    den = 1.0 + K3*A + K1*KD*r*r + KD*KD*K1*K2*r*r*r*r
    num1=Kt1+Kt3*K1*KD*r*r
    fx= num1/den -alpha*x
    num2= Kt2*K3*A
    fm=num2/den - beta*m
    fr=num2/den - gamma*r
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
plt.axis([0,10,0,5])
plt.scatter (t_list,x_list, s=1)
plt.scatter (t_list,m_list, s=1)
plt.scatter (t_list,r_list, s=1)
plt.legend (("C protein", "Cox protein", "TetR protein"))
plt.title ("C, Cox and TetR protein dynamics over time in plasmid B")
plt.xlabel ("Time")
plt.ylabel ("Protein concentration")
plt.show()

# Plot f(x) over x, f(m) over x and f(r) over x in order to find the number of steady states #
# The number of times f(x) and f(m) cross horizontal axis reflects the number of steady states and the exact concentration of the steady state #

x=0
m=0
r=0
dx=0.001
dm=0.001
dr=0.001
x_list=list()
m_list=list()
r_list=list()
fx_list=list()
fm_list=list()
fr_list=list()

while x<5:
    den = 1.0 + K3*A + K1*KD*r*r + KD*KD*K1*K2*r*r*r*r
    num1=Kt1+Kt3*K1*KD*r*r
    fx= num1/den -alpha*x
    num2= Kt2*K3*A
    fm=num2/den - beta*m
    fr=num2/den - gamma*r
    x=x+dx
    x_list.append(x)
    fx_list.append(fx)
    m=m+dm
    m_list.append(m)
    fm_list.append(fm)
    r=r+dr
    r_list.append(r)
    fr_list.append(fr)


fig = plt.figure()
plt.axis([0,5,-5,5])
plt.scatter (x_list,fx_list,s=1)
plt.scatter (m_list,fm_list,s=1)
plt.scatter (r_list,fr_list,s=1)
plt.legend (("C protein", "Cox protein", "TetR protein"))
plt.axhline(y=0, color='r', linestyle=':')
plt.title ("Steady state analysis in plasmid B")
plt.xlabel ("Protein concentration")
plt.ylabel ("dx/dt")
plt.show()