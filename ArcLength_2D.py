import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, root
#import sympy as sp
import scipy.linalg as slin
import time
import math

def trails():
    l_max, u_max = 10,10
    f_max = 120
    r = 0.15
    ui = 2
    lj, uj = 0,0

    f1 = lambda uj: -0.06*uj**3 + 1.2*uj**2 + 3*uj 
    f2 = lambda lj, uj: (((lj - li)**2/l_max**2) + ((uj - ui)**2/u_max**2))**0.5
    fi = f1(ui)
    li = fi/f_max

    def equations(vars):
        lj, uj = vars
        eqn1 = f2(lj,uj) - r
        eqn2 = f1(uj) - lj*f_max 
        return [eqn1, eqn2]
    solution = fsolve(equations, [lj,uj])
    x = np.linspace(0,30,200)

    for i in range(30):
        solution = root(equations, (lj,uj), method = "lm")
        lj, uj = solution.x

        li, ui = lj, uj

        fi = f1(ui)
        plt.scatter(ui, fi)
        plt.text(ui,fi, f"{i}")
        #print(ui, li, fi)
    fx = f1(x)
    plt.plot(x,fx)
    plt.show()


#trails()
area= np.array([1,1,1,1,1,1]) #m**2
youngs_modulus= 20.0 #N/m**2
yield_stress= 1 #N\m**2
const= area*youngs_modulus


cord = np.array([[0.0,0.0],[4.0,0.0],[2.0,4.0],[0.0,4.0], [4.0,4.0], [2.0,7.0]])
con = np.array([[0,2],[1,2],[2,5],[3,5],[4,5]])
#con = np.array([[0,1],[2,1]])
bc_lst = [0,1 , 2,3 ,4, 6,7 , 8,9, 10]
force= [[2,[0,-0.2]]]



cord = np.array([[0.0,0.0],[4.0,6.0],[8.0,0.0]])
con = np.array([[0,1],[1,2], [0,2]])
bc_lst = [0,1 , 5]
force= [[1,[0,-0.1]]]
#force= [[1,[-0.1,0]]]



area= np.ones([1,16])[0] #m**2
const= area*youngs_modulus

""" cord= np.array([[0,0],[1.800,3.118], [3.600,0.0], [5.400,3.118], [7.200,0.0],[9.000,3.118],[10.800,0.0],[12.600,3.118],[14.400,0.0]]) 
con= np.array([[0,1],[0,2],[1,2],[1,3],[2,3],[2,4],[3,4],[3,5],[4,5],[4,6],[5,6],[6,7],[7,8],[6,8],[5,7]])
bclst= [0,0, 1,1, 1,1, 1,1, 1,1, 1,1, 1,1, 1,1, 0,0] 
bc_lst= [0,1 , 16,17]
force= [[1,[0,-0.1]], [3,[0,-0.1]], [5,[0,-0.1]], [7,[0,-0.1]]] """



incre = 2000
plots = 10
it_max = 8
Ss = 0.15
offset = 0
ylim = 17
dof = 3

def force_matrix(node,force):
    f_matrix= np.zeros([2*len(cord),1])
    f_matrix[2*(node)] += force[0]
    f_matrix[2*(node)+1] += force[1]
    return f_matrix

force_mat=force_matrix(0,[0,0])
for ac in force:
    force_mat += force_matrix(ac[0],ac[1])

force_matc =force_matrix(0,[0,0])
for ac in force:
    force_matc += force_matrix(ac[0],ac[1])

def plotter(coordinate):
    for i,k in enumerate(con):
        k=[coordinate[k[0]],coordinate[k[1]]]
        c=[k[0][0],k[1][0]]
        d=[k[0][1],k[1][1]]
        plt.plot(c,d,color=(i/len(con), 1-i/len(con)**2, i/len(con)**3))
        #plotter(cord)
#plotter(cord)
def ln(node1,node2,corl):
    n1= node1
    n2= node2
    return abs((( corl[(n2),0] - corl[(n1),0] )**2 + ( corl[(n2),1] - corl[(n1),1] )**2 )**0.5)

#print(ln(2,3))
def cos(node1,node2,type):
    if node1>node2:
        #node1,node2 = node2,node1
        pass
    return ((type[node2,0]-type[node1,0])/ln(node1,node2,type))

#print(cos(1,2))
def sin(node1, node2,type):
    if node1>node2:
        #node1,node2 = node2,node1
        pass
    return ((type[node2,1]-type[node1,1])/ln(node1,node2,type))

dcord = np.copy(cord)


def local_stiff_mat(element, p, cordt, dcordt):
    #strain = np.array([np.log(ln(*i,dcordt)/ln(*i,cordt)) for i in con]) #np.array([((ln(*i,cordt) - ln(*i,dcordt))/(ln(*i,cordt))) for i in con])
    strain = np.array([((ln(*i,dcordt) - ln(*i,cordt))/(ln(*i,cordt))) for i in con])
    #q = const[p]*strain
    #q = np.array([const[j]*((ln(*i,dcordt) - ln(*i,cordt))/(ln(*i,cordt))) for j,i in enumerate(con)])
    q = np.array([const[j]*(ln(*i,dcordt)**2/ln(*i,cordt)**2 - 1)/2 for j,i in enumerate(con)])


    l= cos(*element, dcordt)
    m= sin(*element, dcordt)
    #print(l, "CX", m, "=CY", q, 'Q')
    matrix= (const[p]/ln(*element, cordt))*np.array([[  l**2 ,  l*m  , -l**2 , -l*m  ],
                    [  l*m  ,  m**2 , -l*m  , -m**2 ],
                    [ -l**2 , -l*m  ,  l**2 ,  l*m  ],
                    [ -l*m  , -m**2 ,  l*m  ,  m**2 ]]
                )+  (q[p]/ln(*element, dcordt))*np.array([[1,0,-1,0],
                                                        [0,1,0,-1],
                                                        [-1,0,1,0],
                                                        [0,-1,0,1]]) 

    #print(matrix, "local martrix") 
    return matrix



""" (q[p]/ln(*element, dcord))*np.array([[-m**2 , l*m , m**2 , -l*m],
                                [l*m , -l**2 , -l*m , l**2],
                                [m**2 , -l*m , -m**2 , l*m],
                                [-l*m , l**2 , l*m , -l**2]]) """

def gsm(coord_list, cordt, dcordt):
    matrix= np.zeros([len(cordt)*2,len(cordt)*2])
    global e_len
    e_len= []
    for p,i in enumerate(coord_list):
        t1= np.array([[i[0],i[0]], [i[0],i[1]], [i[1],i[0]], [i[1],i[1]]])
        t2= np.array([[0,0], [0,2], [2,0], [2,2]])
        k= local_stiff_mat(i,p,cordt, dcordt)
        for j in range(4):
            matrix[2*t1[j,0]:2*t1[j,0]+2 , 2*t1[j,1]:2*t1[j,1]+2] += k[t2[j,0]:t2[j,0]+2 , t2[j,1]:t2[j,1]+2]
        e_len.append(ln(*i,dcordt))
    np.set_printoptions(precision=2)
    return matrix



def color_plotter(coordinate, strs, alp, off, no):
    fig,ax = plt.subplots()
    ax.set_facecolor('black')
    ax.set_ylim(-ylim, ylim)
    ax.set_xlim(-1,20)
    zscale = 1
    j=0
    for k in con:
        k=[coordinate[k[0]],coordinate[k[1]]]
        c= np.array([k[0][0],k[1][0]])
        d= np.array([k[0][1],k[1][1]])
        
        fos= strs[j]/yield_stress
        if fos**2>=1:
            r,g,b,a = 0 ,1 ,0 ,1
            plt.annotate('material will fail at this load',(3000,-5000),color= 'white')
            ax.plot(c+off,d,color=(r,g,b,a), alpha = alp)
        elif fos<=0:
            r,g,b,a = 1+fos, 1+fos, 1,  1
            ax.plot(c+off,d,color=(r,g,b,a), alpha = alp)
        elif fos>0:
            r,g,b,a= 1, 1-fos, 1-fos,  1         
            ax.plot(c+off,d,color=(r,g,b,a), alpha = alp)
        j+=1
    fig.savefig(f'C:/Users/Karthikeyan/Desktop/Sim_images/Schur2D{no}.png', dpi = 100)
    plt.close()


def bc(matrix,bc_list):
    matrix[bc_list , :] =0
    if len(matrix[0])!=1:
        matrix[: , bc_list] =0
        for i in bc_list:
            matrix[i,i]=1
    np.set_printoptions(precision=2)    
    return matrix

def internal_forces(cordt, dcordt):
    #strain = np.array([np.log(ln(*i,dcordt)/ln(*i,cordt)) for i in con]) #((ln(*i,dcordt) - ln(*i,cordt))/(ln(*i,cordt)))
    #strain = np.array([((ln(*i,dcordt) - ln(*i,cordt))/(ln(*i,cordt))) for i in con])
    #strain = np.array([(ln(*i,dcordt)**2/ln(*i,cordt)**2 - 1)/2 for i in con]) # GREEN LAGRANGE STRAIN

    #q = const*strain
    #q = np.array([const[j]*((ln(*i,dcordt) - ln(*i,cordt))/(ln(*i,cordt))) for j,i in enumerate(con)])
    q = np.array([const[j]*(ln(*i,dcordt)**2/ln(*i,cordt)**2 - 1)/2 for j,i in enumerate(con)])

    trans_mats = np.array([[cos(*b,dcordt),sin(*b,dcordt),-cos(*b,dcordt),-sin(*b,dcordt)] for b in con])
    #print(trans_mats)
    #print(q)
    g_nodal_force = np.array([trans_mats[c]*q[c] for c in range(len(con))])
    #print(g_nodal_force,"gnodal")

    internal_force=force_matrix(0,[0,0])
    for f,val in enumerate(g_nodal_force):
        co = con[f] 
        if co[0]>co[1]:
            #co[0],co[1] = co[1],co[0]
            pass
        internal_force += force_matrix(co[0],val[:2])
        internal_force+=force_matrix(co[1],val[2:])

    
    return internal_force,q



def stress_matrix(cordt, dcordt):
    
    def_e_len= []
    e_len = []
    for i in con:
        def_e_len.append(ln(*i,dcordt))
        e_len.append(ln(*i,cordt))

    #print(def_e_len)
    global strain
    strain= []
    for j in range(len(e_len)):
        strain.append(((def_e_len[j]**2 / e_len[j]**2) -1)/e_len[j])
    stress= [youngs_modulus*k for k in strain]
    #np.set_printoptions(precision= 2)
    return stress



def alpha_range1(min,n,x, pow):
    a,c = np.linalg.inv(np.array([[1,1],[n**pow, 1]]))@np.array([min,1])
    y= a*x**pow + c
    if y>1:
        return 1
    return y


Fext = bc(force_mat,bc_lst)
#dU_ = slin.solve(gsM , Fext)
#print(dU_)
#ta = time.time()
#alpha = 1
#print(np.zeros([2*len(cord),1]))
#del_L = 0
#del_U = np.zeros([2*len(cord),1])
#Lo += del_L

###################################

disp = np.zeros([2*len(cord),1])

dispPrev  = disp
dispPrev2 = disp

dcord = cord + disp.copy().reshape(len(cord), 2)
Kglobal = bc(gsm(con,cord, dcord), bc_lst)
Rglobal = internal_forces(cord , dcord)[0]




loadincr = 0.01
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

    dU__ = slin.solve(gsM , bc(Residual, bc_lst))
    dU_ = slin.solve(gsM , Fext)
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
    if Du[3] >0:
        pass
        #print(inc)
    Dl = loadfactor - loadfactorPrev

    convergedPrev = converged
    converged = 0
    #plt.axis("equal")
    #plt.scatter(-disp[dof], loadfactor, color = 'blue', s=0.6)
    #ITERATIONS


    for it in range(it_max+1):
        #print('it', it)
        dcord = cord + disp.copy().reshape(len(cord), 2)
        Kglobal = bc(gsm(con,cord, dcord), bc_lst)  
        Rglobal = internal_forces(cord , dcord)[0]

        Rglobal = Rglobal + loadfactor*Fext
        
        dU, dL, converged = arc_length_eqn(inc, it, it_max, Kglobal, Rglobal, Fext, Du, Dl, Ds)


        #plt.scatter(-disp[3], loadfactor, color = 'blue')
        #plt.scatter(-disp[7], loadfactor)

        if converged:
            break

        disp += dU
        loadfactor += dL

        Du += dU
        Dl +=dL        

    if inc% plots == 0:
        stresses.append(stress_matrix(cord,dcord)) 
        displacements.append(dcord)
        plt.scatter(-disp[dof], loadfactor, color = 'red', s=0.6)


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


l = len(stresses)
off = 0
for g in range(l):
    continue
    range1 = alpha_range1(0.1, l, g+1, 1)
    color_plotter(displacements[g],stresses[g] , 1, off, g)
    off+= offset


#plt.axis("equal")

ax.set_facecolor('black')
plt.show()
