##RUN THIS CODE AS IT IS IN YOUR SYSTEM##

import numpy as np
import matplotlib.pyplot as plt
import math
import sys
fig, ax = plt.subplots()
a= 645.2e-6 #m**2
youngs_modulus= 70e9 #N/m**2
yield_stress= 15000e6 #N\m**2
const= float(a*youngs_modulus)
min_op = 0.3
tolerance = 0.0001

#INPUT COORDINATES
cord= np.array([[0,0],[1800,3118], [3600,0], [5400,3118], [7200,0],[9000,3118],[10800,0],[12600,3118],[14400,0]]) #CONNECTIONS
con=[[0,1],[0,2],[1,2],[1,3],[2,3],[2,4],[3,4],[3,5],[4,5],[4,6],[5,6],[6,7],[7,8],[6,8],[5,7]] #ELEMENT CONNECTIONS
lst= [0,0, 1,1, 1,1, 1,1, 1,1, 1,1, 1,1, 1,1, 0,0] #BOUNDARY CONDITION FOR EACH DOF 0 FOR FIX 1 FOR FREE
force= [[3,[0,-4000e3]], [6,[0,-4000e3]]] #X AND Y DIRECTIONS OF FORCE AT NECESSARY NODE INDEX WRT TO CORD


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

def ln(node1,node2,corl):
    n1= node1
    n2= node2
    return abs(( ( corl[(n2),0] - corl[(n1),0] )**2 + ( corl[(n2),1] - corl[(n1),1] )**2 )**0.5)

def cos(node1,node2,type):
    return ((type[node2,0]-type[node1,0])/ln(node1,node2,type))

def sin(node1, node2,type):
    return ((type[node2,1]-type[node1,1])/ln(node1,node2,type))

dcord = np.array([i for i in cord])


def local_stiff_mat(element, p, cordt, dcordt):
    strain = np.array([((ln(*i,cordt) - ln(*i,dcordt))/(ln(*i,cordt))) for i in con])
    q = const*strain
    l= cos(*element, dcordt)
    m= sin(*element, dcordt)
    matrix= (const/ln(*element, cord))*np.array([[  l**2 ,  l*m  , -l**2 , -l*m  ],
                    [  l*m  ,  m**2 , -l*m  , -m**2 ],
                    [ -l**2 , -l*m  ,  l**2 ,  l*m  ], #LINEAR STIFFNESS MATRIX
                    [ -l*m  , -m**2 ,  l*m  ,  m**2 ]]
                )+  (q[p]/ln(*element, dcord))*np.array([[-m**2 , l*m , m**2 , -l*m],
                                [l*m , -l**2 , -l*m , l**2],
                                [m**2 , -l*m , -m**2 , l*m], #GEOMETRIC STIFFNESS MATRIX
                                [-l*m , l**2 , l*m , -l**2]])
    return matrix


def global_stiff_mat(coord_list, cordt, dcordt):
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


fig,ax = plt.subplots()
def color_plotter(coordinate, strs, alp):
    j=0
    for k in con:
        k=[coordinate[k[0]],coordinate[k[1]]]
        c=[k[0][0],k[1][0]]
        d=[k[0][1],k[1][1]]
        
        fos= strs[j]/yield_stress
        if fos**2>=1:
            r,g,b,a = 0 ,1 ,0 ,1
            plt.annotate('material will fail at this load',(3000,-5000),color= 'white')
            ax.plot(c,d,color=(r,g,b,a), alpha = alp)
        elif fos<=0:
            r,g,b,a = 1+fos, 1+fos, 1,  1
            ax.plot(c,d,color=(r,g,b,a), alpha = alp)
        elif fos>0:
            r,g,b,a= 1, 1-fos, 1-fos,  1         
            ax.plot(c,d,color=(r,g,b,a), alpha = alp)
        j+=1
    ax.annotate(f'No.of iterations = {l}',(5500,1800),color= 'white')
    ax.annotate(f'Convergence criteria = {0.00001}',(4700,1500),color= 'white')
def boundary_condtions(matrix,bc_list):
        global del_lst
        del_lst= []
        for i in range(0,len(lst)):
            if bc_list[i]==0:
                del_lst.append(i) 
            else:
                continue
        matrix= np.delete(matrix,del_lst,0)
        if len(matrix[0])!=1: #TO NEGLECT COLUMN DELETION FOR FORCE MATRIX BC
            matrix= np.delete(matrix,del_lst,1)

        np.set_printoptions(precision=4)
        return matrix


def residual(cordt, dcordt):
    strain = np.array([((ln(*i,dcordt) - ln(*i,cordt))/(ln(*i,cordt))) for i in con])
    q = const*strain

    trans_mats = np.array([[-cos(*b,dcordt),-sin(*b,dcordt),cos(*b,dcordt),sin(*b,dcordt)] for b in con])
    g_nodal_force = np.array([trans_mats[c]*q[c] for c in range(len(con))])
    residual_force=force_matrix(0,[0,0])
    for f,val in enumerate(g_nodal_force):
        co = con[f] 
        residual_force += force_matrix(co[0],val[:2])
        residual_force+=force_matrix(co[1],val[2:])
    return residual_force,q

def stress_matrix(cordt, dcordt):
    
    def_e_len= []
    e_len = []
    for i in con:
        def_e_len.append(ln(*i,dcordt))
        e_len.append(ln(*i,cordt))
    global strain
    strain= []
    for j in range(len(e_len)):
        strain.append((def_e_len[j]-e_len[j])/e_len[j])
    stress= [youngs_modulus*k for k in strain]
    return stress



l=0
tol = 8
displacements = []
stresses = []
residual_force, fn = residual(cord,dcord)

ele_force1 = []
disps1 = []
tol_lst = []
while tol > tolerance:
    print(f"ITERATION {l}")
    global_stiff_matrix= global_stiff_mat(con,cord, dcord)
    bc_stiff_mat= boundary_condtions(global_stiff_matrix,lst)
    bc_force_matrix= boundary_condtions(force_mat,lst)

    #INVERTING STIFFNESS MATRIX
    bc_stiff_mat_inv= np.linalg.inv(bc_stiff_mat)
    #MULTIPLYING INVERSE TO FORCE MATRIX
    nodal_disp= bc_stiff_mat_inv@bc_force_matrix

    for i in del_lst:
        nodal_disp= np.insert(nodal_disp,i,0)
    nodal_disp1 = nodal_disp.copy()

    ######
    disps = np.array([(ln(*i,dcord) - ln(*i,cord)) for i in con])
    strain = np.array([((ln(*i,dcord) - ln(*i,cord))/(ln(*i,cord))) for i in con])
    el_force = strain*const
    #####

    cord, dcord = cord, dcord+nodal_disp1.reshape(len(cord), 2)
    displacements.append(dcord)
    stresses.append(stress_matrix(cord, dcord))

    #####
    disps = np.array([(ln(*i,dcord) - ln(*i,cord)) for i in con])
    print(el_force[2], disps[2])
    ele_force1.append(el_force)
    disps1.append(disps)
    #####

    residual_force, fn = residual(cord,dcord)
    bc_int_force_mat = boundary_condtions(residual_force, lst)

    force_mat = force_matc - residual_force 

    tol = math.sqrt(np.sum(nodal_disp)**2 / np.sum(dcord-cord)**2)
    print("TOLERANCE: ",tol)
    tol_lst.append(tol)
    l+=1

def convergence():
    fig,ax = plt.subplots()
    ax.scatter(np.array([i for i in range(len(tol_lst))]), np.array(tol_lst))
    ax.plot(np.array([i for i in range(len(tol_lst))]), np.array(tol_lst))
    ax.set_xlabel("Iterations")
    ax.set_ylabel("Convergence Number")
    ax.set_title("Convergence progression")
convergence()

alpha_range = lambda min, n, x: 1.0 if ((1-min)/(n-1))*x + (n*min-1)/(n-1)>1 else ((1-min)/(n-1))*x + (n*min-1)/(n-1)

def alpha_range1(min,n,x, pow):
    a,c = np.linalg.inv(np.array([[1,1],[n**pow, 1]]))@np.array([min,1])
    y= a*x**pow + c
    if y>1:
        return 1
    return y



for g in range(l):
    range = alpha_range1(min_op, l, g+1, 1)
    color_plotter(displacements[g],stresses[g] , range)
    print(range)


def sparser():
    oy = []
    for q in range(len(cord)):
        for w in range(len(cord)):
            if not np.array_equal(global_stiff_matrix[6*q:6*q+6 , 6*w:6*w+6] , np.zeros([6,6])):
                oy.append([w,q])

    print(f'Sparsity = {1- (len(oy)/len(cord)**2)}')
    plt.scatter(np.array(oy)[:,0] , np.array(oy)[:,1],s=1)

    jk = 0
    plt.xlim(0-jk,len(cord)+jk)
    plt.ylim(0-jk,len(cord)+jk)
    plt.gca().invert_yaxis()
    plt.show()




ax.set_facecolor('black')
plt.axis("equal")
plt.show()
