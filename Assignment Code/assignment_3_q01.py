
def Thomas_Algorithm(num, dia, upp, low, rhs, sol):     # Thomas Algorithm
    dia1    = np.zeros(num)
    rhs1    = np.zeros(num)
    dia1[0] = dia[0]
    rhs1[0] = rhs[0]
    for i in range(1, len(dia)):
        dia1[i] = dia[i] - low[i]*upp[i-1]/dia1[i-1]
        rhs1[i] = rhs[i] - rhs1[i-1]*low[i]/dia1[i-1]
    sol[num-1] = rhs[num-1]/dia[num-1]
    for i in range(len(dia)-2,-1,-1):
        sol[i] = (rhs1[i]-upp[i]*sol[i+1])/dia1[i]


import numpy as np
import math as m
import matplotlib.pyplot as plt
import sys

num_mesh = 21                       # number of mesh points
xmesh    = np.zeros(num_mesh)       # Mesh point
del_x    = 1.0/(num_mesh-1)         # mesh size (Delta_x)

# Compute the location of mesh points
for i in range(0, len(xmesh)):       
    xmesh[i] = i * del_x

N= 1000                             #Number of terms in series
Time = 0.5                          #Total Time

T_1 = np.zeros(N)                   # Individual terms of series       
T_2 = np.zeros(N)                   # Summation of the terms     
T = np.zeros(num_mesh)              #Analytical solution

pie = m.pi                      
T_i = 100                           #initial temperature in the wall
T_s = 300                           #Equal surface temperature at the two faces
a=1                                 #Thermal Diffusivity of the 1-D slab


#===========Calculation of Analytical solution============#
for j in range(0,num_mesh):
   for i in range(1,N):
      T_1[i] = (np.exp(-(a*Time*((i*pie)**2)))*((1-((-1)^i))/(i*pie))*(np.sin(i*pie*xmesh[j])))        #i'th value of the series 
      T_2[i] = T_2[i-1] + T_1[i]                                                                       #Addition of all the values till i
   T[j] = T_s + 2*(T_i-T_s)*T_2[N-1]                                                                   #Calculation of analytical solution at j'th mesh point


#===========Explicit Solution============#
N_E = 4000                             # Number of iterations for Explicit Method

T_new_E = np.zeros(num_mesh)  # Solution at new time step
T_old_E = np.zeros(num_mesh)  #Solution at previous time step

t_E = np.linspace(0, Time, N_E)    # mesh points in time(each in 0.000125 hr)
dt_E = t_E[1] - t_E[0]              # Calculation of timestep 
f_E = a*dt_E/((del_x)**2)


#=========Applying Boundary Condition for Explicit part===========#
for i in range(0,num_mesh):
  if i==0 or i==num_mesh-1:
    T_old_E[i]=300.0
    T_new_E[i]=300.0
  else:
    T_old_E[i] = 100.0



#=============Solution of Explicit part==================#
for n in range(0,N_E):
  for i in range(1,num_mesh-1):
    T_new_E[i] = T_old_E[i] + f_E*(T_old_E[i-1]-2*T_old_E[i]+T_old_E[i+1])    #Equation for calculating Temperature at next timestep                                                   
  T_old_E[:] = T_new_E[:]                                                     #Updating the previous timestep with current timestep

  if((N_E==m.trunc((0.1/dt_E)))):                                                                 #Value of Temperature at different mesh points for time = 0.1 hr
    T_01 = np.zeros(num_mesh)
    T_01[:] = T_new_E[:]
  elif((N_E==m.trunc((0.2/dt_E)))):                                                              #Value of Temperature at different mesh points for time = 0.2 hr
    T_02 = np.zeros(num_mesh)
    T_02[:] = T_new_E[:]
  elif((N_E==m.trunc((0.3/dt_E)))):                                                              #Value of Temperature at different mesh points for time = 0.3 hr
    T_03 = np.zeros(num_mesh)
    T_03[:] = T_new_E[:]
  elif((N_E==m.trunc((0.4/dt_E)))):                                                              #Value of Temperature at different mesh points for time = 0.4 hr
    T_04 = np.zeros(num_mesh)
    T_04[:] = T_new_E[:]



#==========Graph of Temperatures at different mesh points at each 0.1 hr interval=========#
plt.plot(xmesh,T_01,'-*',label = 'T =0.1 hr')
plt.plot(xmesh,T_02,'-*',label = 'T =0.2 hr')
plt.plot(xmesh,T_03,'-*',label = 'T =0.3 hr')
plt.plot(xmesh,T_04,'-*',label = 'T =0.4 hr')
plt.plot(xmesh,T_new_E,'-*',label = 'T =0.5 hr')
plt.xlabel('length')
plt.ylabel('Temperature')
plt.legend()
plt.show()

#==============================THOMAS_ALGORITHM for Crank Nicolson Scheme==========================#
N_CN = 4000                          # Number of iterations for Crank Nicolson Method

T_new_CN = np.zeros(num_mesh)  # Solution at new time step
T_old_CN = np.zeros(num_mesh)  #Solution at previous time step

# Diagonal elements of system matrix for crank Nicolson Scheme
d_CN    = np.zeros(num_mesh)        # main diagonal elements
u_CN    = np.zeros(num_mesh)        # upper diagonal
l_CN    = np.zeros(num_mesh)        # lower diagonal
rhs_CN    = np.zeros(num_mesh)

t_CN = np.linspace(0, Time,N_CN)    # mesh points in time(each in 0.000125 hr)
dt_CN = t_CN[1] - t_CN[0]
f_CN = a*dt_CN/((del_x)**2) 

#=========Applying Boundary Condition for Crank-Nicolson part and updating the diagonal,lower and upper values of the matrix===========#
for i in range(0,num_mesh):
  d_CN[i] = 1.0
  if i==0 or i==num_mesh-1:
    T_old_CN[i]=300.0
    T_new_CN[i]=300.0
    l_CN[i] = 0.0
    u_CN[i] = 0.0
    rhs_CN[i] = 300.0
  else:
    T_old_CN[i] = 100.0
    l_CN[i] = -(f_CN/(2*(1+f_CN)))
    u_CN[i] = -(f_CN/(2*(1+f_CN)))


#=============Solution of Crank-Nicolson part==================#
for n in range(0,N_CN):
    for j in range(1,num_mesh-1):                                 #Updating RHS at each timestep
        rhs_CN[j] = (f_CN/(2*(1+f_CN)))*T_old_CN[j+1] + (f_CN/(2*(1+f_CN)))*T_old_CN[j-1] + ((1-f_CN)/(1+f_CN))*T_old_CN[j]
    Thomas_Algorithm(num_mesh,d_CN,u_CN,l_CN,rhs_CN,T_new_CN)     #Thomas algorithm to evaluate solution
    T_old_CN[:] = T_new_CN[:]                                     #Updating the previous timestep with current timestep

    if((N_CN==m.trunc((0.1/dt_CN)))):                                                   #Value of Temperature at different mesh points for time = 0.1 hr
      T_01 = np.zeros(num_mesh)
      T_01[:] = T_new_CN[:]
    elif((N_CN==m.trunc((0.2/dt_CN)))):                                                #Value of Temperature at different mesh points for time = 0.2 hr
      T_02 = np.zeros(num_mesh)
      T_02[:] = T_new_CN[:]
    elif((N_CN==m.trunc((0.3/dt_CN)))):                                                #Value of Temperature at different mesh points for time = 0.3 hr
      T_03 = np.zeros(num_mesh)
      T_03[:] = T_new_CN[:]
    elif((N_CN==m.trunc((0.4/dt_CN)))):                                                #Value of Temperature at different mesh points for time = 0.4 hr
      T_04 = np.zeros(num_mesh)
      T_04[:] = T_new_CN[:]

#==========Graph of Temperatures at different mesh points at each 0.1 hr interval=========#
plt.plot(xmesh,T_01,'-*',label = 'T =0.1 hr')
plt.plot(xmesh,T_02,'-*',label = 'T =0.2 hr')
plt.plot(xmesh,T_03,'-*',label = 'T =0.3 hr')
plt.plot(xmesh,T_04,'-*',label = 'T =0.4 hr')
plt.plot(xmesh,T_new_CN,'-*',label = 'T =0.5 hr')
plt.xlabel('length')
plt.ylabel('Temperature')
plt.legend()
plt.show()


#===============Thomas Algorithm for Fully Implicit Scheme===============#
N_I = 4000                           # Number of iterations for Implicit Method

T_new_I = np.zeros(num_mesh)  # Solution at new time step
T_old_I = np.zeros(num_mesh)  #Solution at previous time step



# Diagonal elements of system matrix for Fully Implcit Scheme
d_I   = np.zeros(num_mesh)        # main diagonal elements
u_I    = np.zeros(num_mesh)        # upper diagonal
l_I    = np.zeros(num_mesh)        # lower diagonal
rhs_I    = np.zeros(num_mesh)

t_I = np.linspace(0, Time, N_I)    # mesh points in time(each in 0.000125 hr)
dt_I= t_I[1] - t_I[0]              # Calculation of timestep 
f_I = a*dt_I/((del_x)**2)

#=========Applying Boundary Condition for Crank-Nicolson part and updating the diagonal,lower and upper values of the matrix===========#
for i in range(0,num_mesh):
  d_I[i] = 1.0
  if i==0 or i==num_mesh-1:
    T_old_I[i]=300.0
    T_new_I[i]=300.0
    l_I[i] = 0.0
    u_I[i] = 0.0
    rhs_I[i] = 300.0
  else:
    T_old_I[i] = 100.0
    l_I[i] = -(f_I/(1+(2*f_I)))
    u_I[i] = -(f_I/(1+(2*f_I)))

#=============Solution of Fully Implicit part==================#   
for n in range(0,N_I):
    for j in range(1,num_mesh-1):                             #Updating RHS at each timestep
        rhs_I[j] = (1/(1+(2*f_I)))*T_old_I[j]
    Thomas_Algorithm(num_mesh,d_I,u_I,l_I,rhs_I,T_new_I)      #Thomas algorithm to evaluate solution
    T_old_I[:] = T_new_I[:]                                   #Updating the previous timestep with current timestep
    
    if((N_I==m.trunc((0.1/dt_I)))):                                               #Value of Temperature at different mesh points for time = 0.1 hr 
      T_01 = np.zeros(num_mesh)
      T_01[:] = T_new_I[:]
    elif((N_I==m.trunc((0.2/dt_I)))):                                            #Value of Temperature at different mesh points for time = 0.2 hr
      T_02 = np.zeros(num_mesh)
      T_02[:] = T_new_I[:]
    elif((N_I==m.trunc((0.3/dt_I)))):                                            #Value of Temperature at different mesh points for time = 0.3 hr
      T_03 = np.zeros(num_mesh)
      T_03[:] = T_new_I[:]
    elif((N_I==m.trunc((0.4/dt_I)))):                                            #Value of Temperature at different mesh points for time = 0.4 hr
      T_04 = np.zeros(num_mesh)
      T_04[:] = T_new_I[:]
    
#==========Graph of Temperatures at different mesh points at each 0.1 hr interval=========#
plt.plot(xmesh,T_01,'-*',label = 'T =0.1 hr')
plt.plot(xmesh,T_02,'-*',label = 'T =0.2 hr')
plt.plot(xmesh,T_03,'-*',label = 'T =0.3 hr')
plt.plot(xmesh,T_04,'-*',label = 'T =0.4 hr')
plt.plot(xmesh,T_new_I,'-*',label = 'T =0.5 hr')
plt.xlabel('length')
plt.ylabel('Temperature')
plt.legend()
plt.show()

#==========Graph of Temperatures at different mesh points at each 0.5 hr=========#
plt.plot(xmesh,T,'r-o',marker = '^')
plt.plot(xmesh,T_new_E,'g-o', marker='v')
plt.plot(xmesh,T_new_CN,'b-o', marker='*')
plt.plot(xmesh,T_new_I,'y-o', marker='+')
plt.xlabel('x')
plt.ylabel('Temperature')
plt.show()
