

#==============Calculation of sum absolute values of row elements of Error matrix============#
def norm_value(Error,x,y):
    for c in range(0,x):
      row_elements = np.zeros(y)
      row_elements_1 = np.zeros(y)
      for e in range(0,y):
        row_elements_1[e] = Error[c][e]
        row_elements[e] = row_elements[e-1] + np.absolute(row_elements_1[e])
      row_elements_2[c] = row_elements[y-1]
    return(row_elements_2)

#=================Thomas Algorithm=======================#
def Thomas_Algorithm(num, dia, upp, low, rhs):                # Thomas Algorithm
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
    return(sol)

import numpy as np
import math as m
import matplotlib.pyplot as plt
import sys

# Set colour interpolation and colour map.
# You can try set it to 10, or 100 to see the difference
# You can also try: colourMap = plt.cm.coolwarm
colorinterpolation = 50
colourMap = plt.cm.jet

L = 0.3           # Length in x direction
H = 0.4           # Height in y direction
k = 380           # Thermal Conductivity
diff = 0.00011234 # Thermal diffusivity
pie = m.pi

num_mesh_x = 31                      # number of mesh points
num_mesh_y = 41                      # number of mesh points
del_x    = L/(num_mesh_x-1.0)        # mesh size (Delta_x)
del_y    = H/(num_mesh_y-1.0)        # mesh size (Delta_y)
del_t    = 1

xmesh    = np.zeros(num_mesh_x)        # Mesh point
ymesh    = np.zeros(num_mesh_y)       # Mesh point
tmesh    = np.zeros(1165)

T = np.zeros((num_mesh_x,num_mesh_y))                    #Analytical solution
T_1 = np.zeros((num_mesh_x,num_mesh_y))
T_2 = np.zeros((num_mesh_x,num_mesh_y))
T_3 = np.zeros((num_mesh_x,num_mesh_y))
T_4 = np.zeros((num_mesh_x,num_mesh_y))

# Compute the location of mesh points
for i in range(0, len(xmesh)):
    xmesh[i] = i * del_x
for i in range(0, len(ymesh)):
    ymesh[i] = i * del_y
for i in range(0,len(tmesh)):
    tmesh[i] = i*del_t

#========Boundary Conditions===========#
Tabove = 10.0
Tbelow = 40.0
Tright = 0.0
Tleft = 0.0

#===============Analytical Solution===============#
for i in range(0,len(xmesh)-1):
  for j in range(0,len(ymesh)-1):
    for n in range(1,51):
      T_1[i,j] = ((1-((-1)**n))*(np.sinh(n*pie*(H-ymesh[j])/L))*(np.sin(n*pie*xmesh[i]/L)))/((np.sinh(n*pie*H/L))*(n*pie))
      T_2[i,j] = T_2[i][j] + T_1[i][j]
      T_3[i,j] = ((1-((-1)**n))*(np.sinh(n*pie*ymesh[j]/L))*(np.sin(n*pie*xmesh[i]/L)))/((np.sinh(n*pie*H/L))*(n*pie))
      T_4[i,j] = T_4[i][j] + T_3[i][j]
    T[i,j] = 2*Tbelow*T_2[i][j] + 2*Tabove*T_4[i][j]


#===========Discretization===============#
T_time = np.zeros((len(tmesh),num_mesh_x,num_mesh_y))
T_new = np.zeros((num_mesh_x,num_mesh_y))             # Temperature at current timestep
T_old = np.zeros((num_mesh_x,num_mesh_y))             # Temperature at previous timestep
T_plate = np.zeros((num_mesh_x,num_mesh_y))           # Temperature used for convergence of solution obtained due to Thomas Algorithm
sol = np.zeros(num_mesh_y)         
Error = np.zeros((num_mesh_x,num_mesh_y))             # Error between the previous timestep and current timestep
Error_1 = np.zeros((num_mesh_x,num_mesh_y))           # Error between T_new and T_plate at same timestep
row_elements_2 = np.zeros(num_mesh_y)
row_elements_4 = np.zeros(num_mesh_y)

#=========Matrix elements==========#
l = np.zeros(num_mesh_y)
u = np.zeros(num_mesh_y)
d = np.zeros(num_mesh_y)
rhs = np.zeros(num_mesh_y)
ap_0 = del_x*del_y/(del_t*diff)

#=============Boundary Conditions========#
for i in range(0,len(xmesh)):
  for j in range(0,len(ymesh)):
    if j == 0:
      T_new[i,j] = Tbelow
      T_old[i,j] = Tbelow
    elif j == len(ymesh)-1 :
      T_new[i,j] = Tabove
      T_old[i,j] = Tabove
    elif i == 0:
      T_new[i,j] = Tright
      T_old[i,j] = Tright
    elif i == len(xmesh)-1:
      T_new[i,j] = Tleft
      T_old[i,j] = Tleft
    else:  
      T_old[i,j] = 0
T_plate[:] = T_old[:]

#=========Initializing the matrix elements===========#

for j in range(0,num_mesh_y):
  if j == 0:
    d[j] = 1
    l[j] = 0
    u[j] = 0
    rhs[j] = 40
  elif j == num_mesh_y-1:
    d[j] = 1
    l[j] = 0
    u[j] = 0
    rhs[j] = 10
  else:
    l[j] = -1
    u[j] = -1
    d[j] = del_x*del_y/(del_t*diff) + 4
    
iter = 0
#===============Convergence + Solution Loop===============#
while True:
        while True:                                                              # This loop converges the solution at each timestep
              for i in range(1,len(xmesh)-1):                                    # Initializing RHS and updating  the RHS after each timestep
                for j in range(1,len(ymesh)-1):
                  rhs[j] = T_new[i-1][j] + T_new[i+1][j] + ap_0*T_old[i][j]
                  sol = Thomas_Algorithm(num_mesh_y,d,u,l,rhs)                   # Thomas Algorithm 
                for b in range(0,num_mesh_y):
                  T_new[i,b] = sol[b]
              Error_1[:][:] = T_new[:][:] - T_plate[:][:]                        # Error in the same time step
              row_elements_2 = norm_value(Error_1,len(xmesh),len(ymesh))         # Calculation of row norm of the matrix for same timestep
              norm_1 = np.max(row_elements_2)
              if(norm_1<0.001):                                                      # Convergence of Temperature in same timestep
                break                                              
              T_plate[:][:] = T_new[:][:]
        Error[:] = T_new[:] - T_old[:]                                           # Error between previous timestep and current timestep     
        row_elements_4 = norm_value(Error,len(xmesh),len(ymesh))                 # Calculation of row norm of the matrix between previous timestep and current timestep
        norm = np.max(row_elements_4)
        if(norm<0.000001):                                                              # Convergence of Temperature in between previous timestep and current timestep
          break 
        for i in range(0,len(xmesh)):
          for j in range(0,len(ymesh)):
            T_time[iter,i,j] = T_new[i][j]   
        iter = iter + 1
        T_old[:][:] = T_new[:][:]
print("Total number of ietrations :")
print(iter)


# Configure the contour
plt.title("Contour of Temperature(Line by Line)")
plt.contourf(xmesh, ymesh, np.transpose(T_new), colorinterpolation, cmap=colourMap)
plt.colorbar()                       # Set Colorbar
plt.show()                           # Show the result in the plot window

# Configure the contour
plt.title("Contour of Temperature(Analytical)")
plt.contourf(xmesh, ymesh, np.transpose(T), colorinterpolation, cmap=colourMap)
plt.colorbar()                       # Set Colorbar
plt.show()                           # Show the result in the plot window 
  
for k in range(1,9):
    plt.title(str(round((iter)*k/28))  + "sec")
    plt.contourf(xmesh, ymesh, np.transpose(T_time[round(iter*k/28)]), colorinterpolation, cmap=colourMap)
    plt.colorbar()                       # Set Colorbar
    plt.show()                           # Show the result in the plot window
