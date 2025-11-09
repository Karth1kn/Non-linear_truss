import numpy as np
import matplotlib.pyplot as plt
#from scipy.optimize import fsolve, root
#import sympy as sp
import scipy
import scipy.linalg as slin
import time
import math
from scipy.linalg import solveh_banded

from D3_TrussNLN import bc, gsm, internal_forces, dofs, cord, disp, con, bc_lst, force_mat, bandwidth, color_plotter, alpha_range1, stress_matrix

print(cord)

dcord = np.copy(cord)

incre = 200
plots = 15
it_max = 5
Ss = 0.05
offset = 0

dof = 8

Fext = bc(force_mat,bc_lst)

dispPrev  = disp
dispPrev2 = disp

dcord = np.copy(cord)
Kglobal = bc(gsm(con,cord, dcord), bc_lst)
Rglobal = internal_forces(cord , dcord)[0]


loadincr = 0.1
Ds = loadincr
DsPrev = Ds
DsMax = Ds
DsMin = Ds

loadfactor = Ss
loadfactorPrev2 = 0.0
loadfactorPrev = 0.0

converged = 0
convergedPrev = 0

fig,ax = plt.subplots()

psi = 1.0
def arc_length_eqn(increment, iteration, it_max, gsM, Residual, Fext, Du, Dl, ds):
    psi = 1.0

    if increment>0:
        Asp = (Du[:,0]@Du)[0] + psi*(Dl**2)*(Fext[:,0]@Fext)[0] - ds**2
        a = 2.0*Du
        b = 2.0*psi*Dl*((Fext[:,0]@Fext)[0])
    else:
        Asp = 0.0
        a = 0*Du
        b = 1.0

    #diagonals = np.array([np.append(np.diag(gsM, -i), np.zeros([1,i])[0]) for i in range(bandwidth)  ])
    #dU__ = solveh_banded(diagonals, bc(Residual, bc_lst), lower= True)
    #dU_ = solveh_banded(diagonals, Fext, lower= True)

    dU__ = scipy.sparse.linalg.splu(gsM).solve(bc(Residual, bc_lst))    #slin.solve(gsM , bc(Residual, bc_lst))
    dU_ = scipy.sparse.linalg.splu(gsM).solve(Fext)     #slin.solve(gsM , Fext)
    dU__ = -1*dU__


    dL = ((a[:,0]@dU__)[0] - Asp)/(b + (a[:,0]@dU_)[0])
    dU = -dU__ + dL*dU_

    converged = 0
    #tol = math.sqrt(np.sum(dU)**2 / np.sum(Du)**2)
    #if tol<10e-6:
    if iteration == it_max:
        converged=1
        #print("converged")


    return dU, dL, converged


displacements = []
stresses = []
#LOAD INCREMENTS
for inc in range(incre):
    #print('inc', inc)
    if inc>0:
        DsFactor1 = Ds/DsPrev
        #DsFactor1 = 1
        disp     = (1.0+DsFactor1)*dispPrev - DsFactor1*dispPrev2
        loadfactor = (1.0+DsFactor1)*loadfactorPrev - DsFactor1*loadfactorPrev2

    Du = disp - dispPrev

    Dl = loadfactor - loadfactorPrev

    convergedPrev = converged
    converged = 0
    #plt.axis("equal")
    plt.scatter(-disp[dof], loadfactor, color = 'blue', s=0.6)
    #plt.scatter(disp[22], disp[17], color = 'blue', s=0.6)
    #ITERATIONS
    for it in range(it_max+1):
        #print('it', it)
        dcord = cord + disp.copy().reshape(len(cord), dofs)
        Kglobal = bc(gsm(con,cord, dcord), bc_lst)  
        Rglobal = internal_forces(cord , dcord)[0]

        Rglobal = Rglobal + loadfactor*Fext
        
        dU, dL, converged = arc_length_eqn(inc, it, it_max, Kglobal, Rglobal, Fext, Du, Dl, Ds)


        if converged:
            break

        disp += dU
        loadfactor += dL

        Du += dU
        Dl +=dL        


    if inc% plots == 0:
        plt.scatter(-disp[dof], loadfactor, color = 'red', s=4)
        stresses.append(stress_matrix(cord,dcord)) 
        displacements.append(dcord)

    if converged:
        #plt.text( -disp[7], loadfactor, f"{inc}")


        if inc == 0:
            Ds = ((Du[:,0]@Du)[0] + psi*(loadfactor**2)*(Fext[:,0]@Fext)[0])**0.5
            DsMax = Ds
            DsMin = Ds/1024


        loadfactorPrev2, loadfactorPrev = loadfactorPrev, loadfactor
        dispPrev2, dispPrev = dispPrev, disp
        DsPrev = Ds
        #cord += disp.copy().reshape(len(cord), 2)

        if convergedPrev:
            Ds = min(max(2.0*Ds, DsMin), DsMax)
        

    else:

        if convergedPrev:
            Ds = max(Ds*0.5, DsMin)
        else:
            Ds = max(Ds*0.25, DsMin)

#plt.show()
zscale = 3

l = len(stresses)
off = 0
for g in range(l):
    #continue
    #fig,ax = plt.subplots()
    range1 = alpha_range1(0.1, l, g+1, 1)
    color_plotter(displacements[g],stresses[g] , 1, g, zscale)
    
    off+= offset


#plt.axis("equal")
ax.set_facecolor('black')
plt.xlabel('Displacement')
plt.ylabel('Load Factor ')
plt.show()

