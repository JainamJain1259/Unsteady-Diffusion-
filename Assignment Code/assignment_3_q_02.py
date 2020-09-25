
def norm_value(Error,x,y):
    for c in range(0,len(x)):
      row_elements = np.zeros(num_mesh_y)
      row_elements_1 = np.zeros(num_mesh_y)
      for e in range(0,len(y)):
        row_elements_1[e] = Error[c][e]
        row_elements[e] = row_elements[e-1] + np.absolute(row_elements_1[e])
      row_elements_2[c] = row_elements[len(y)-1]
    return(row_elements_2)


import numpy as np
import math as m
import matplotlib.pyplot as plt
import sys

L=1.0                                 #Given length in x direction
H=2.0                                #Given length in y direction
N =1000                               #Number of terms taken in series
pie = m.pi

num_mesh_x = 21                       #Number of mesh points in x direction 
del_x    = L/(num_mesh_x-1.0)         # mesh size (Delta_x)
num_mesh_y = 41                       # number of mesh points in y direction
del_y    = H/(num_mesh_y-1.0)         # mesh size (Delta_y)

xmesh    = np.zeros(num_mesh_x)       # Mesh point
ymesh    = np.zeros(num_mesh_y)       # Mesh point

Tbelow=100
Tabove=0
Tright=0
Tleft=0
Tguess=50

# Set colour interpolation and colour map.
# You can try set it to 10, or 100 to see the difference
# You can also try: colourMap = plt.cm.coolwarm
colorinterpolation = 50
colourMap = plt.cm.jet

# Compute the location of mesh points
for i in range(0, len(xmesh)):
    xmesh[i] = i * del_x
for i in range(0, len(ymesh)):
    ymesh[i] = i * del_y

T1 = 100.0                                  #Temperature at y=0
T = np.zeros((num_mesh_x,num_mesh_y))       #Analytical solution
T_1 = np.zeros((num_mesh_x,num_mesh_y))     #2-D arrays used for summation of
T_2 = np.zeros((num_mesh_x,num_mesh_y))     #1000 terms of the series 

#===========Calculation of Analytical solution============#

for i in range(0,len(xmesh)-1):
  for j in range(0,len(ymesh)-1):
    for n in range(1,51):
      T_1[i,j] = ((1-((-1)**n))*(np.sinh((n*pie*(H-ymesh[j])/L)))*(np.sin(n*pie*xmesh[i]/L)))/((np.sinh(n*pie*H/L))*(n*pie))
      T_2[i,j] = T_2[i][j] + T_1[i][j]
    T[i,j] = 2*T1*T_2[i][j]
#print(T)

# Configure the contour
plt.title("Contour of Temperature(Analytical)")
plt.contourf(xmesh, ymesh, np.transpose(T), colorinterpolation, cmap=colourMap)

# Set Colorbar
plt.colorbar()

# Show the result in the plot window
plt.show()


#==================Steady State Solution using Gauss Siedel Algorithm====================#
a = del_x/del_y
b = del_y/del_x

# Set array size and set the interior value with Tguess
T_S = np.zeros((len(xmesh), len(ymesh)))          #Temperature of new value of the slab calculated using Seidel Method at current timestep
T_old = np.zeros((len(xmesh), len(ymesh)))        #Temperature of new value of the slab calculated using Seidel Method at current timestep
T_plate = np.zeros((len(xmesh), len(ymesh)))      #Guessed Temperature matrix
Error = np.zeros((len(xmesh), len(ymesh)))        #Error between the previous timestep and current timestep
Error_it = np.zeros((len(xmesh), len(ymesh))) 
row_elements_1 = np.zeros(num_mesh_y)
row_elements_2 = np.zeros(num_mesh_y)             

#===========Application of Boundary Condition===========#
for i in range(0,len(xmesh)):
  for j in range(0,len(ymesh)):
    if j == 0:
      T_old[i,j] = Tbelow
      T_S[i,j] = Tbelow
    elif j == len(ymesh)-1 :
      T_old[i,j] = Tabove
      T_S[i,j] = Tabove
    elif i == 0:
      T_old[i,j] = Tleft
      T_S[i,j] = Tleft
    elif i == len(xmesh)-1:
      T_old[i,j] = Tright
      T_S[i,j] = Tright
    else:  
      T_S[i,j] = Tguess

T_S[0,0] = 0
T_S[len(xmesh)-1,0]= 0
converged = False
iter = 0

while converged == False:
  iter = iter + 1

  for i in range(1,len(xmesh)-1):
    for j in range(1,len(ymesh)-1):
      T_S[i,j] = (b*T_S[i+1][j] + b*T_S[i-1][j] + a*T_S[i][j+1] + a*T_S[i][j-1])/(2*(a+b))
    
  Error[:] = T_S[:]- T_old[:]
  row_elements_2 = norm_value(Error,xmesh,ymesh)

  norm = np.max(row_elements_2) 

  if(norm<0.001):
    converged = True

  T_old[:] = T_S[:]   
print(iter)


# Configure the contour
plt.title("Contour of Temperature(Using Gauss-Seidel Method)")
plt.contourf(xmesh, ymesh, np.transpose(T_S), colorinterpolation, cmap=colourMap)

# Set Colorbar
plt.colorbar()

# Show the result in the plot window
plt.show()

