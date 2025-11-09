#INSTALL MAYAVI AND MATPLOTLIB
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as slin
import math 
import csv
import time
import os
import sys
from scipy.linalg import solveh_banded
#from mayavi import mlab

def timer(func):
    def wrapper(*args, **kwargs):
        t1 = time.time()
        mat = func(*args, **kwargs)
        t2 = time.time()
        print(f'{func.__name__} finished in: {t2-t1}')
        return mat
    return wrapper

script_dir = os.path.dirname(os.path.abspath(__file__))
coord_path = os.path.join(script_dir, 'truss_coord.csv')
con_path = os.path.join(script_dir, 'truss_connect.csv')
prop_path = os.path.join(script_dir, 'truss_properties.csv')
bc_path = os.path.join(script_dir, 'truss_bound_con.csv')
force_path = os.path.join(script_dir, 'truss_force_list.csv')


l1 = time.time()


with open(coord_path, 'r') as file:
    try:
        val_list = []
        for i in csv.reader(file):
            val_list.append([float(j) for j in i])
        val_array= np.array(val_list)
        with open(con_path,'r') as con:
            con_list= []

            for c in csv.reader(con):
                con_list.append([int(b) for b in c])
            con_array= np.array(con_list)

            for k in con_array:
                x_ele = [val_array[k[0],0],val_array[k[1],0]]
                y_ele = [val_array[k[0],1],val_array[k[1],1]]
                z_ele = [val_array[k[0],2],val_array[k[1],2]]
                #line1 = mlab.plot3d(x_ele, y_ele, z_ele, color=(0,0,0))

    except Exception:
        print('error')

cord=val_array
con= con_array
disp = np.zeros([3*len(cord),1])
dofs = 3

#### BOUNDARY CONDITIIONS
with open(bc_path, 'r') as bc:
        bc_lst = []
        c = 0
        for c1 in csv.reader(bc):
            pos = int(c1[0])*3
            if c1[1] == 'a':
                for r1 in range(3):
                    bc_lst.append(pos+r1)
                    c= c1[0]
                continue
            elif c1[1]=='x':
                if c1[0]==c:
                    continue
                bc_lst.append(pos)
            elif c1[1] == 'y':
                if c1[0]==c:
                    continue
                bc_lst.append(pos+1)
            elif c1[1] == 'z':
                if c1[0]==c:
                    continue
                bc_lst.append(pos+2)

with open(force_path, 'r') as force:
        force_list =[]
        lst = csv.reader(force)
        for row in lst:
            if row == []: continue
            for i in row[3:]:
                force_list.append([int(i),[float(j) for j in row[:3]]])

with open(prop_path,'r')as fl:
    props = []
    for i in csv.reader(fl):
        for j in i:
            props.append(float(j))
    area =props[1]  
    youngs_modulus = props[0] 
    yield_stress = props[2]
    tolerance = props[3]
    min_op = props[4]
    plop = props[5]



const = youngs_modulus*area
ln = lambda n1,n2,type: np.clip((( ( type[(n2),0] - type[(n1),0] )**2 + ( type[(n2),1] - type[(n1),1] )**2 + ( type[(n2),2] - type[(n1),2] )**2 )**0.5), a_max=np.inf, a_min=0)

op = lambda node1,node2, type,i: ((type[node2,i]-type[node1,i])/ln(node1,node2,type))
cosx = lambda node1,node2,type: op(node1,node2,type,0)
cosy = lambda node1,node2,type: op(node1,node2,type,1)
cosz = lambda node1,node2,type: op(node1,node2,type,2)

np.set_printoptions(linewidth=np.inf)



def force_matrix(node,force):
    f_matrix= np.zeros([3*len(cord),1])
    f_matrix[3*(node)] += force[0]
    f_matrix[3*(node)+1] += force[1]
    f_matrix[3*(node)+2] += force[2]
    return f_matrix

force_matc =force_matrix(0,[0,0,0])
for ac in force_list:
    force_matc += force_matrix(ac[0],ac[1])

force_mat=force_matrix(0,[0,0,0])
for ac in force_list:
    force_mat += force_matrix(ac[0],ac[1])

def local_stiff_mat(element, p, cordt, dcordt):
    l, m, n = cosx(*element, dcordt), cosy(*element, dcordt), cosz(*element, dcordt)
    #strain = np.array([((ln(*i,dcordt) - ln(*i,cordt))/(ln(*i,cordt))) for i in con]) #LINEAR STRAIN
    #strain = np.array([(ln(*i,dcordt)**2/ln(*i,cordt)**2 - 1)/2 for i in con]) # GREEN LAGRANGE STRAIN        
    strain = (ln(*con[p],dcordt)**2/ln(*con[p],cordt)**2 - 1)/2  # GREEN LAGRANGE STRAIN        

    kk=1
    if p ==4:
        kk=1
    q = const*strain

    matrix = ((const/kk)/ln(*element, cordt))*np.array([ [ l**2 ,  l*m  ,  l*n  , -l**2 ,  -l*m , -l*n  ],
                        [ l*m  ,  m**2 ,  m*n  ,  -l*m , -m**2 , -m*n  ],
                        [ l*n  ,   m*n ,  n**2 ,  -l*n ,  -m*n , -n**2 ],
                        [-l**2 ,  -l*m , -l*n  ,  l**2 ,  l*m  ,  l*n  ],
                        [ -l*m , -m**2 , -m*n  ,  l*m  ,  m**2 ,  m*n  ],
                        [ -l*n ,  -m*n , -n**2 ,  l*n  ,  m*n  ,  n**2 ] ]
                        )+(q/ln(*element, dcordt))*np.array([[1,0,0,-1,0,0],
                                                                [0,1,0,0,-1,0],
                                                                [0,0,1,0,0,-1],
                                                                [-1,0,0,1,0,0],
                                                                [0,-1,0,0,1,0],
                                                                [0,0,-1,0,0,1],]) 

    return matrix


def linear_local_stiff_mat(element, p, cord):
    l, m, n = cosx(*element, cord), cosy(*element, cord), cosz(*element, cord)
    #strain = np.array([((ln(*i,dcordt) - ln(*i,cordt))/(ln(*i,cordt))) for i in con]) #LINEAR STRAIN
    #strain = np.array([(ln(*i,dcordt)**2/ln(*i,cordt)**2 - 1)/2 for i in con]) # GREEN LAGRANGE STRAIN        
    kk=1
    if p ==4:
        kk=1

    matrix = ((const/kk)/ln(*element, cord))*np.array([ [ l**2 ,  l*m  ,  l*n  , -l**2 ,  -l*m , -l*n  ],
                        [ l*m  ,  m**2 ,  m*n  ,  -l*m , -m**2 , -m*n  ],
                        [ l*n  ,   m*n ,  n**2 ,  -l*n ,  -m*n , -n**2 ],
                        [-l**2 ,  -l*m , -l*n  ,  l**2 ,  l*m  ,  l*n  ],
                        [ -l*m , -m**2 , -m*n  ,  l*m  ,  m**2 ,  m*n  ],
                        [ -l*n ,  -m*n , -n**2 ,  l*n  ,  m*n  ,  n**2 ] ]
                        )
    return matrix



def gsm(coord_list, cordt, dcordt):
    matrix= np.zeros([len(cordt)*3,len(cordt)*3])
    global e_len
    e_len= []
    for p,i in enumerate(coord_list):
        t1= np.array([[i[0],i[0]], [i[0],i[1]], [i[1],i[0]], [i[1],i[1]]])
        t2= np.array([[0,0], [0,3], [3,0], [3,3]])
        k= local_stiff_mat(i,p,cordt, dcordt)
        for j in range(4):
            matrix[3*t1[j,0]:3*t1[j,0]+3 , 3*t1[j,1]:3*t1[j,1]+3] += k[t2[j,0]:t2[j,0]+3 , t2[j,1]:t2[j,1]+3]
        e_len.append(ln(*i,dcordt))

    np.set_printoptions(precision=2)
    return matrix

def Lgsm(coord_list, cordt):
    matrix= np.zeros([len(cordt)*3,len(cordt)*3])
    global e_len
    e_len= []
    for p,i in enumerate(coord_list):
        t1= np.array([[i[0],i[0]], [i[0],i[1]], [i[1],i[0]], [i[1],i[1]]])
        t2= np.array([[0,0], [0,3], [3,0], [3,3]])
        k= linear_local_stiff_mat(i,p,cordt)

        for j in range(4):
            matrix[3*t1[j,0]:3*t1[j,0]+3 , 3*t1[j,1]:3*t1[j,1]+3] += k[t2[j,0]:t2[j,0]+3 , t2[j,1]:t2[j,1]+3]
        e_len.append(ln(*i,cordt))

    np.set_printoptions(precision=2)
    return matrix


#fig = plt.figure(figsize=(10,10))

def local_mass_matrix(element, rho, cord):
    lnn = ln(*element, cord)

    consistant_mass = (rho*area*lnn/6)*np.array([[2,0,0,1,0,0],
                                          [0,2,0,0,1,0],
                                          [0,0,2,0,0,1],
                                          [1,0,0,2,0,0],
                                          [0,1,0,0,2,0],
                                          [0,0,1,0,0,2]])

    return consistant_mass


def gmm(coord_list, cordt, rho):
    matrix= np.zeros([len(cordt)*3,len(cordt)*3])
    global e_len
    e_len= []
    for p,i in enumerate(coord_list):
        t1= np.array([[i[0],i[0]], [i[0],i[1]], [i[1],i[0]], [i[1],i[1]]])
        t2= np.array([[0,0], [0,3], [3,0], [3,3]])
        k= local_mass_matrix(i, rho, cord)
        for j in range(4):
            matrix[3*t1[j,0]:3*t1[j,0]+3 , 3*t1[j,1]:3*t1[j,1]+3] += k[t2[j,0]:t2[j,0]+3 , t2[j,1]:t2[j,1]+3]


    np.set_printoptions(precision=2)
    return matrix



def color_plotter(coordinate, strs, alp, no, zscale):
    fig = plt.figure(figsize=(25,15))
    ax = fig.add_subplot(111, projection= '3d', aspect = 'equal')

    global t3
    t3 = time.time()
    #ax.set_aspect('equal')
    plt.axis("equal")
    ax.set_facecolor('black')
    ax.grid(visible=False)
    azi = 60
    for j,k in enumerate(con):
        ks = np.array([coordinate[k[0]],coordinate[k[1]]])
        c = np.array([ks[0,0],ks[1,0]])
        d = np.array([ks[0,1],ks[1,1]])
        e = np.array([ks[0,2],ks[1,2]])

        fos= (strs[j])/(yield_stress)
        #print(fos)
        if fos**2>=1:
            r,g,b,a = 0 ,1 ,0 ,alp
            if plop ==0:
                ax.plot3D(c,d,e,color=(r,g,b,a))
            else: 
                #line1 = mlab.plot3d(c,d,e,color=(r,g,b))
                #line1.actor.property.opacity = alp
                pass

        elif fos<=0:
            r,g,b,a = 1+fos, 1+fos, 1,  alp
            if plop ==0:
                ax.plot3D(c,d,e,color=(r,g,b,a))
            else: 
                #line1 = mlab.plot3d(c,d,e,color=(r,g,b))
                #line1.actor.property.opacity = alp
                pass
                
        elif fos>0:
            r,g,b,a= 1, 1-fos, 1-fos,  alp         
            if plop ==0:
                ax.plot3D(c,d,e,color=(r,g,b,a))
            else: 
                #line1 = mlab.plot3d(c,d,e,color=(r,g,b))
                #line1.actor.property.opacity = alp
                pass

    ax.set_zlim(-zscale, zscale)
    ax.view_init(elev=10., azim=azi-no)
    #fig.savefig(f'C:/Users/Karthikeyan/Desktop/Sim_images/Kadapa_{no}.png')
    plt.close()
    #main_sequence()

    #ax.set_facecolor('black')

    #ax.set_box_aspect([1, 1, 1])
    global t4
    t4 = time.time()

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
    strain = np.array([(ln(*i,dcordt)**2/ln(*i,cordt)**2 - 1)/2 for i in con]) # GREEN LAGRANGE STRAIN

    q = const*strain

    #q = np.array([const*((ln(*i,dcordt) - ln(*i,cordt))/(ln(*i,cordt))) for j,i in enumerate(con)])

    trans_mats = np.array([[cosx(*b,dcordt),cosy(*b,dcordt),cosz(*b,dcordt),-cosx(*b,dcordt),-cosy(*b,dcordt),-cosz(*b,dcordt)] for b in con])
    #print(trans_mats)
    #print(q)
    g_nodal_force = np.array([trans_mats[c]*q[c] for c in range(len(con))])
    #print(g_nodal_force,"gnodal")

    internal_force=force_matrix(0,[0,0,0])
    for f,val in enumerate(g_nodal_force):
        co = con[f] 
        if co[0]>co[1]:
            #co[0],co[1] = co[1],co[0]
            pass
        internal_force += force_matrix(co[0],val[:3])
        internal_force+=force_matrix(co[1],val[3:])

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
        strain.append((def_e_len[j]-e_len[j])/e_len[j])
    stress= [youngs_modulus*k for k in strain]
    #np.set_printoptions(precision= 2)
    return stress
def rnge():
    max_arr = [abs(i[0]-i[1]) for i in con_array ]
    return max(max_arr)*3+3

bandwidth = rnge()   
#print(bandwidth)    



@timer
def main_sequence():
    dcord = np.array([i for i in cord])
    force_mat=force_matrix(0,[0,0,0])
    for ac in force_list:
        force_mat += force_matrix(ac[0],ac[1])
    rnge = lambda:  max([abs(i[0]-i[1]) for i in con_array ])*3+3
    iters = 0
    displacements = []
    stresses = []
    tol_lst = []
    tol = 10
    print("**** GEOMETRIC NON-LINEAR FEA SOLVER ****")
    while tol > tolerance:
        print(f"ITERATION {iters}")
        global global_stiff_matrix
        global_stiff_matrix= gsm(con,cord, dcord)
        bc_stiff_mat= bc(global_stiff_matrix,bc_lst)
        bc_force_matrix= bc(force_mat,bc_lst)
        bandwidth = rnge()   
        diagonals = np.array([np.append(np.diag(bc_stiff_mat, -i), np.zeros([1,i])[0]) for i in range(bandwidth)  ])
        nodal_disp = solveh_banded(diagonals, bc_force_matrix, lower= True)
        #nodal_disp= slin.solve(bc_stiff_mat,bc_force_matrix )
        nodal_disp1 = nodal_disp.copy()
        dcord = dcord+nodal_disp1.reshape(len(cord), 3)
        displacements.append(dcord)
        stresses.append(stress_matrix(cord, dcord))
        residual_force = internal_forces(cord,dcord)
        force_mat = force_matc - residual_force
        tol = math.sqrt(np.sum(nodal_disp)**2 / np.sum(dcord-cord)**2)

        print("TOLERANCE: ",tol)
        tol_lst.append(tol)
        iters+=1

    return displacements, stresses, tol_lst

#main_sequence()
def convergence(tol_lst):
    fig,ax = plt.subplots()
    ax.scatter(np.array([i for i in range(len(tol_lst))]), np.array(tol_lst))
    ax.plot(np.array([i for i in range(len(tol_lst))]), np.array(tol_lst))
    

alpha_range = lambda min, n, x: 1.0 if ((1-min)/(n-1))*x + (n*min-1)/(n-1)>1 else ((1-min)/(n-1))*x + (n*min-1)/(n-1)

def alpha_range1(min,n,x, pow):
    a,c = np.linalg.inv(np.array([[1,1],[n**pow, 1]]))@np.array([min,1])
    y = a*x**pow + c
    if y>1:
        return 1.00
    return y




if __name__ == "__main__":
    # Check if the script is run as the main program
    if len(sys.argv) > 1:
        function_name = sys.argv[1]
        if function_name == "main_sequence":
            
            displacements, stresses, tol_lst = main_sequence()
            iterate = len(tol_lst)
            for g in range(iterate):
                rng = alpha_range1(min_op, iterate, g+1, 1)
                color_plotter(displacements[g],stresses[g] , rng)

                print(rng)
            if plop ==0:
                plt.show()
            else: 
                pass
                #mlab.show()
            print(f'Plotting time: {t4-t3}')
            convergence(tol_lst)
            #plt.show()
            

def sparser():
    fig1,ax1 = plt.subplots()
    
    oy = []
    global_stiff_matrix= gsm(con,cord, cord)
    for q in range(len(cord)):
        for w in range(len(cord)):
            if not np.array_equal(global_stiff_matrix[3*q:3*q+3 , 3*w:3*w+3] , np.zeros([3,3])):
                oy.append([w,q])
    ax1.scatter(np.array(oy)[:,0] , np.array(oy)[:,1],s=1)
    print(f'Sparsity = {1- (len(oy)/len(cord)**2)}')
    plt.annotate(f'Sparsity = {1- (len(oy)/len(cord)**2)})',(int(len(cord)/2), 0),color= 'white')
    jk=0
    plt.xlim(0-jk,len(cord)+jk)
    plt.ylim(0-jk,len(cord)+jk)
    plt.gca().invert_yaxis()
    plt.axis('equal')
    plt.show()


#sparser()
if __name__ == "__main__":
    if len(sys.argv) > 1:
        function_name = sys.argv[1]
        if function_name == "sparser":
            sparser()
