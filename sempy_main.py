# Ross Maguire

# Main program for 1D spectral element wave equation.  Reads user input
# files and calls subroutines to build and solve problem.

#------------------------------------------------------------------------------
from assemble_matrix import *
from sem_class import *
from obspy.signal.filter import bandpass
import matplotlib.pyplot as plt

print "WE'RE IN BUISNESS"

nt = 10000 #number of time steps
#dt = 0.00038971143 #time step
dt = 2e-4 #time step
amp = 1e7 #source amplitude

a = sem_1d()
i_source = int(a.n_global_nodes/2) #place source in middle
#i_source = 200

#make source--------------------------------------------------------------------
box = np.ones(10000)
time = np.linspace(0,1000,10000)
source = bandpass(box,freqmin=1/13.0,freqmax=1/1.0,df=1/0.03,corners=4)
source = source*amp
source_int1 = np.cumsum(source)
source_int2 = np.cumsum(source_int1)
fig,axes = plt.subplots(3)
axes[0].plot(source)
axes[1].plot(source_int1)
axes[2].plot(source_int2)
plt.show()

#build matrices-----------------------------------------------------------------
K_e = build_local_stiffness_matrix(a) #local stiffness matrix
M_e = build_local_mass_matrix(a) #local mass matrix
M = build_global_mass_matrix(a,M_e)

#setup field variables----------------------------------------------------------
u = np.zeros(a.n_global_nodes)
vel = np.zeros(a.n_global_nodes)
acc = np.zeros(a.n_global_nodes)

plot = True
if plot:
   plt.ion()
   fig = plt.figure()
   plt.hold(False) 
   plt_step = 30

#try newmark scheme-------------------------------------------------------------
for it in range(0,nt):

   print 'time = ', it * dt

   if it > 0:
      u += dt * vel + acc * dt**2 / 2
      vel += dt / 2 * acc
      acc.fill(0)
   count = 0

   for e in range(0,a.n_elem):
      ind_1 = count
      ind_2 = count+a.lpd+1
      #print 'index = ', ind_1,ind_2
      acc[ind_1:ind_2] -= np.dot(K_e[e,:,:], u[ind_1:ind_2])
      #print 'max acc_e = ', max(acc)
      #if e < a.n_elem:
      count += a.lpd
      #else:
      #   count += a.lpd
   #print 'source amp = ',source[it]
   #print 'maximum u = ', max(u)

   acc[i_source] += source[it]
   acc /= M
   vel += dt / 2 * acc

   if plot and it % plt_step == 0:
      plt.plot(u)
      plt.ylim([-0.010,0.010])
      plt.title('displacement')
      plt.draw()
