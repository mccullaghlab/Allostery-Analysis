
# PREAMBLE:

import numpy as np
#from numpy.linalg import *
from scipy.linalg import block_diag

# ----------------------------------------
def perform_henm(cartesian_covariance_matrix,average_node_positions,output_directory,guess=None,max_iterations=100,alpha=0.1,kBT = 0.592,threshold=1e-4,distance_cutoff = 9999.,plotting_boolean=False,output_step=10):
        """
        """
        # ----------------------------------------
        # ASSIGNING MATRIX VARIABLES
        # ----------------------------------------
        nCartCoords = cartesian_covariance_matrix.shape[0]
        nNodes = nCartCoords//3
        nNodes_range = range(nNodes)
        inverse_distance_cutoff = 1.0/distance_cutoff

        # ----------------------------------------
        # CALCULATING THE NODE PAIR DISTANCE VECTORS
        # ----------------------------------------
        # code from P. T. Rex Lake
        #make an array of separation vectors
        xhat = average_node_positions[:,None] - average_node_positions[None,:]
        #make an array of distances
        x0 = np.sqrt(np.einsum('ijk,ijk->ij',xhat,xhat))
        #invert said array, avoiding the 0's along the diagonal
        x0 = np.divide(1., x0, out=np.zeros_like(x0), where=x0!=0)
        dist_cutoff_mask = (x0 > inverse_distance_cutoff).astype(float) # value of 1. if node pair distance is less than distance_cutoff; 0 if node pair distance is greater than distance_cutoff
        #convert the separation vectors into unit vectors
        xhat = np.einsum('ijk,ij->ijk', xhat, x0)
        #an array of outer products
        xhat2=np.einsum('ijk,ijl->ijkl',xhat,xhat)
        
        # ----------------------------------------
        # CALCULATING THE INITIAL RESIDUAL FROM THE ORIGINAL COVARIANCE MATRIX
        # ----------------------------------------
        targetResidual = np.zeros((nNodes,nNodes),dtype=np.float64)
        ### NOTE: create code to do this more efficiently/pythonically
        for i in nNodes_range[:-1]:
                iIndex = i*3    # i index assuming that each node has 3 values sequentially populating the covar matrix 
                iIndex3 = iIndex+3
                for j in nNodes_range[i+1:]:
                        jIndex = j*3    # i index assuming that each node has 3 values sequentially populating the covar matrix 
                        jIndex3 = jIndex+3
                        # ii_submatrix + jj_submatrix - ij_submatrix - ji_submatrix; results in a 3x3 matrix that describes the ...
                        temp = cartesian_covariance_matrix[iIndex:iIndex3,iIndex:iIndex3] + cartesian_covariance_matrix[jIndex:jIndex3,jIndex:jIndex3] - cartesian_covariance_matrix[iIndex:iIndex3,jIndex:jIndex3] - cartesian_covariance_matrix[iIndex:iIndex3,jIndex:jIndex3].T
                        # 1/(xhat * covar_submatrix * xhat)
                        targetResidual[i,j] = targetResidual[j,i] = 1.0/(np.dot(xhat[i,j].T,np.dot(temp,xhat[i,j])))
        
        np.savetxt(output_directory + 'target_residual.dat',targetResidual,fmt='%.8e')

        # ----------------------------------------
        # PREP THE INITIAL GUESS OF THE HESSIAN MATRIX; units: kcal mol^-1 \AA^-2
        # ----------------------------------------
        if type(guess) != type(None):
                hessian = np.array(guess)
        else:
                hessian = -10.0*np.ones((nNodes,nNodes),dtype=np.float64)
                hessian -= np.diag(np.diag(hessian))
                hessian -= np.diag(np.sum(hessian,axis=1))
        
        # ----------------------------------------
        # ITERATIVELY IMPROVE THE HESSIAN TO MATCH THE RESIDUAL
        # ----------------------------------------
        step = 0 
        dev = threshold + 9999.
        converged = 'False'
        while step < max_iterations and converged == 'False':
                # project NxN hessian to 3Nx3N hessian
                # code from P. T. Rex Lake
                hessian3N = np.einsum('ij,ijkl->ikjl',hessian,xhat2)
                #set the diagonals
                hessian3N_diag = block_diag(*np.sum(hessian3N,axis=2))
                hessian3N = np.reshape(hessian3N,(nCartCoords,nCartCoords))-hessian3N_diag
                # invert hessian using psuedo inverse; multiplying by thermal energy (kcal mol^-1) to remove energy units
                e,v = np.linalg.eigh(hessian3N)
                covar3N = kBT*np.dot(v[:,6:]/e[6:],v[:,6:].T)
                
                # save model covar3N
                if step%output_step == 0: 
                    np.savetxt(output_directory + '%05d.covar3N.dat'%(step),covar3N,fmt='%.8e')
                
                ### NOTE: create code to do this more efficiently/pythonically
                # take difference of tensor covars and project into separation vector
                diffMatrix = np.zeros((nNodes,nNodes),dtype=np.float64)
                for i in nNodes_range[:-1]:
                        iIndex = i*3    # i index assuming that each node has 3 values sequentially populating the covar matrix 
                        iIndex3 = iIndex+3
                        for j in nNodes_range[i+1:]:
                                jIndex = j*3    # i index assuming that each node has 3 values sequentially populating the covar matrix 
                                jIndex3 = jIndex+3
                                temp = covar3N[iIndex:iIndex3,iIndex:iIndex3] + covar3N[jIndex:jIndex3,jIndex:jIndex3] - covar3N[iIndex:iIndex3,jIndex:jIndex3] - covar3N[iIndex:iIndex3,jIndex:jIndex3].T
                                diffMatrix[i,j] = diffMatrix[j,i] = (1.0/np.dot(xhat[i,j],np.dot(temp,xhat[i,j]))) - targetResidual[i,j]
               
                # save diffMatrix 
                if step%output_step == 0: 
                        np.savetxt(output_directory + '%05d.diffMatrix.dat'%(step),diffMatrix,fmt='%.8e')

                # check to see if converged
                dev = np.linalg.norm(diffMatrix)       # do I apply the distance cutoff on the diffmatrix elements? or on the hessian after I have applied the diff matrix... Doing so would change this value...
                if dev < threshold:
                        converged = 'True'
                else: # update Hessian
                        hessian += alpha * diffMatrix
                        positive_value_mask = (hessian < 0.).astype(float) # 1.0 if less than 0., 0. else where 
                        hessian *= dist_cutoff_mask * positive_value_mask
                        hessian -= np.diag(np.diag(hessian))
                        hessian -= np.diag(np.sum(hessian,axis=1))
                
                print(step, dev)
                step += 1

        # ----------------------------------------
        # RETURNING HESSIAN MATRIX
        # ----------------------------------------
        return hessian

