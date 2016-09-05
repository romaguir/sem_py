# Ross Maguire
# Functions to build and assemble the local and global mass and stiffness matrices

import numpy as np

def build_local_mass_matrix(grid):

    ''' Function to build a mass matrix for each element.
        The entire mass matrix is not calculated since it's
        diagonal. The function returns an array with the 
        diagonals of of the mass matrix for each element. '''

    M_local = np.zeros(((grid.n_elem),(grid.lpd+1)))

    for e in range(0, grid.n_elem):
        for i in range(0,(grid.lpd+1)):
        
            #-----------------------------------------------
            # build mass matrix for each element
            #-----------------------------------------------
            #TODO: replace constant density with variable density
            M_local[e,i] = grid.rho_c * np.inner(grid.gll_wts[i],grid.dx_dxi[i,i])
       
    print M_local
    return M_local
    
def build_global_mass_matrix(grid,local_mass_matrix):
    M_global = np.zeros(grid.n_global_nodes)
    count = 0
    for i in range(0, grid.n_elem):
         M_global[count:count+grid.lpd+1] += local_mass_matrix[i,:]
         count += grid.lpd

    return M_global

def build_local_stiffness_matrix(grid):      
    K_local = np.zeros(((grid.n_elem),(grid.lpd+1),(grid.lpd+1)))

    for e in range(0, grid.n_elem):
        for i in range(0,(grid.lpd+1)):
            for j in range(0,(grid.lpd+1)):

                #-----------------------------------------------
                # build stiffness matrix for each element
                #-----------------------------------------------
                #TODO replace constant shear modulus with variable
            
                d_ij = grid.lp_derivatives[i,:] * grid.lp_derivatives[j,:]
                K_local[e,i,j] =( np.inner(grid.gll_wts, d_ij) * grid.mu_c *
                                 grid.dxi_dx[e,i]**2 * grid.dx_dxi[e,i] )

    return K_local 
