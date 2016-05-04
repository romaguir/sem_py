#Ross Maguire 12-4-15

import matplotlib.pylab as plt
import numpy as np
from scipy.interpolate import lagrange

'''
1d spectral element class
'''

class sem_1d:

    def __init__(self):

        self.lpd            = 4                                         # Lagrange polynomial degree
        self.length         = 3000.0                                    # Domain length
        self.n_elem         = 250                                       # Number of elements
        self.rho_c          = 2500.0                                    # density (kg/m3)
        self.mu_c           = 3.0e10                                    # shear modulus (Pa)
        self.vp_c           = 10.0                                      # vp (km/s)
        self.vs_c           = 5.0                                       # vs (km/s)
        self.elem_omega_pts = np.linspace(0,self.length,self.n_elem+1)  # x coordinates for boundaries between elems
        self.n_global_nodes = self.lpd * self.n_elem + 1                # total number of gll points (i.e., nodes)
        self.x_coor   = np.zeros(self.n_global_nodes)                   # x coordinates of global gll points
        self.rho_grid = np.zeros(self.n_global_nodes)                   # density grid (value at each gll point)
        self.vp_grid  = np.zeros(self.n_global_nodes)                   # vp grid (value at each gll point) 
        self.vs_grid  = np.zeros(self.n_global_nodes)                   # vs grid (value at each point)

        #---------------------------------------------------------------------------------------------------------
        #define coordinates of gll points within reference interval [-1,1]
        #---------------------------------------------------------------------------------------------------------

        #lpd = 4 (only option that is currently implemented)
        if (self.lpd == 4):
            self.gll_pts = np.array([-1.0, -0.6546, 0.0, 0.6546, 1.0])

        #---------------------------------------------------------------------------------------------------------
        #define integration weights 
        #---------------------------------------------------------------------------------------------------------
        
        #lpd = 4 (only option that is currently implemented) 
        if (self.lpd == 4):
            self.gll_wts = np.array([0.1, 0.544444444, 0.711111111, 0.544444444, 0.1])

        #---------------------------------------------------------------------------------------------------------
        #map reference element gll pts into global x coordinates
        #---------------------------------------------------------------------------------------------------------
        count = 0
        for i in range(0,self.n_elem):
           self.x_coor[count:count+5] = ( (( (1 - self.gll_pts) / 2.0) * (self.elem_omega_pts[i])) +
                                          (( (1 + self.gll_pts) / 2.0) * (self.elem_omega_pts[i+1])) )
           count += 4

        #---------------------------------------------------------------------------------------------------------
        #set default mesh grid values
        #---------------------------------------------------------------------------------------------------------
        self.rho_grid[:] = self.rho_c
        self.vp_grid[:]  = self.vp_c
        self.vs_grid[:]  = self.vs_c

        #---------------------------------------------------------------------------------------------------------
        #calculate spatial derivative of element mapping from omega onto the reference element (dx/dxi)
        #---------------------------------------------------------------------------------------------------------
        self.dx_dxi = np.zeros((self.n_elem, self.lpd+1))

        for i in range(0,self.n_elem):
            self.dx_dxi[i,:] = (self.elem_omega_pts[i+1]-self.elem_omega_pts[i])/2.0 #2.0 is lenght of ref interval

        #---------------------------------------------------------------------------------------------------------
        #calculate spatial derivative of element mapping from reference to omega element (dxi/dx)
        #---------------------------------------------------------------------------------------------------------
        self.dxi_dx = np.zeros((self.n_elem, self.lpd+1))

        for i in range(0,self.n_elem):
            self.dxi_dx[i,:] = 2.0/(self.elem_omega_pts[i+1]-self.elem_omega_pts[i]) #2.0 is lenght of ref interval
            

        #---------------------------------------------------------------------------------------------------------
        #make element dictionary
        #---------------------------------------------------------------------------------------------------------
        self.elem_dict    = {'number': np.zeros(self.n_elem), 'x_vals':np.zeros((self.n_elem,self.lpd+1))}
        count = 0
        for i in range(0,self.n_elem):
              self.elem_dict['number'][i] = i
              self.elem_dict['x_vals'][i] = self.x_coor[count:count+5]
              count = count + 4

        #---------------------------------------------------------------------------------------------------------
        #create array of values of lagrange polynomials evaluated at gll points
        #---------------------------------------------------------------------------------------------------------
        self.lp_derivatives = np.zeros(((self.lpd+1),(self.lpd+1)))
        derivs_array = np.zeros(((self.lpd+1),(self.lpd+1)))
        derivs_here = np.zeros((self.lpd+1))
        poly_array  = np.zeros(((self.lpd+1),(self.lpd+1)))
        x           = self.gll_pts
        w           = np.zeros((self.lpd+1))
        for i in range(0, self.lpd+1):
            w[:]    = 0
            w[i]    = 1.0
            l_here  = lagrange(x,w)                                 #numpy.poly1d instance of lagrange coeffs
            dl_here = l_here.deriv()
            for j in range(0,len(self.gll_pts)):
                self.lp_derivatives[i,j] = dl_here(self.gll_pts[j])


    def lagrange_poly(self,plot=True):
        derivs_array = np.zeros(((self.lpd+1),(self.lpd+1)))
        derivs_here = np.zeros((self.lpd+1))
        poly_array  = np.zeros(((self.lpd+1),(self.lpd+1)))
        x           = self.gll_pts
        w           = np.zeros((self.lpd+1))
        for i in range(0, self.lpd+1):
            w[:]    = 0
            w[i]    = 1.0
            l_here  = lagrange(x,w)                                 #numpy.poly1d instance of lagrange coeffs
            dl_here = l_here.deriv()
            for j in range(0,len(self.gll_pts)):
                derivs_array[i,j] = dl_here(self.gll_pts[j])
                
            poly_array[i,:]   = l_here.coeffs[:]

        if(plot==True):
           x_axis = np.linspace(-1,1,100)
           y_vals = np.zeros((len(x_axis)))
           for i in range(0,(self.lpd+1)):
               y_vals = ( poly_array[i,0]*x_axis**(self.lpd)   +
                          poly_array[i,1]*x_axis**(self.lpd-1) +
                          poly_array[i,2]*x_axis**(self.lpd-2) +
                          poly_array[i,3]*x_axis**(self.lpd-3) +
                          poly_array[i,4]                    )
               plt.plot(x_axis, y_vals,'b')
           plt.show()
   
        print derivs_array
        return poly_array

    def plot_mesh(self):
        ax1 = plt.subplot(311)
        ax1.plot(self.x_coor, self.rho_grid,'+')
        ax1.set_ylabel('density (kg/m3)')
        ax2 = plt.subplot(312)
        ax2.plot(self.x_coor, self.vp_grid,'+')
        ax2.set_ylabel('vp (km/s)')
        ax3 = plt.subplot(313)
        ax3.plot(self.x_coor, self.vs_grid,'+')
        ax3.set_ylabel('vs (km/s)')
        ax3.set_xlabel('x (km)')
        plt.show()

    def grid_1d(self,lpd=4, len=1.0, n_elem=100, plot=True):
        self.lpd     = lpd
        self.length  = len
        self.n_elem  = n_elem
        self.gll_pts = np.linspace(0,self.length,(self.n_elem*(self.lpd+1)))
        
        if(plot==True):
           plt.plot(self.gll_pts,self.rho_grid)
           plt.show()

    def inv_elem_mapping(self):
         for i in range(0,self.n_elem):
             coor_here = self.elem_dict['x_vals'][i]
             x_a       = coor_here[0]*1.0
             x_b       = coor_here[self.lpd]*1.0
             xi        = 2.0*((coor_here - (x_a + x_b) /2.0) / (x_b - x_a)) 
             print xi
