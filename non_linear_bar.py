import numpy as np
import matplotlib.pyplot as plt
import math


a= 645.2e-6 #mm**2
youngs_modulus= 70e9 #N/mm**2
yield_stress= 500 #N\mm**2
const= float(a*youngs_modulus)

#INPUT COORDINATES
#cord= np.array([[0,0],[1800,3118], [3600,0], [5400,3118], [7200,0],[9000,3118],[10800,0]]) # [[x1,y1],[x2,y2]...],     coordinate[node,abscissa/ordinate]
cord = np.array([[0.0,0.0],[4.0,3.0],[8.0,0.0]])
#disp = np.zeros(2*len(cord))
#disp1 = np.copy(disp)
#dcord = cord + disp.reshape(len(cord), 2)

#cord = np.array([[0,0],[0,1000],[1000,1000],[1000,0]])
# CONNECTING NODES TO FORM ELEMENTS
m = [[0,1],[1,2],[0,2]]
#m = [[1,0],[2,1],[2,0]]
#m=[[0,1],[0,2],[1,2],[1,3],[2,3],[2,4],[3,4],[3,5],[4,5],[4,6],[5,6]]
#m = [[0,1],[1,2],[2,3],[0,3],[1,3]]
#INPUT BOUNDARY CONDTIONS
#lst= [0,0, 1,1, 1,1, 1,1, 1,1, 1,1, 1,0] #[x1,y1 , x2,y2 , x3,y3 ...]
lst = [0,0 , 1,1 , 1,0 ]
force= [[1,[0,-2000e3]]]
def force_matrix(node,force):
    f_matrix= np.zeros([2*len(cord),1])
    f_matrix[2*(node)] += force[0]
    f_matrix[2*(node)+1] += force[1]
    return f_matrix

force_mat=force_matrix(0,[0,0])
for ac in force:
    force_mat += force_matrix(ac[0],ac[1])
##############################################
#NODAL FORCES
#fig, ax = plt.subplots()

def plotter(coordinate):
    for k in m:
        k=[coordinate[k[0]],coordinate[k[1]]]
        c=[k[0][0],k[1][0]]
        d=[k[0][1],k[1][1]]
        plt.plot(c,d,color=(0,0,0))
    plt.show()

#plotter(cord)


def ln(node1,node2,corl):
    n1= node1
    n2= node2
    return (( ( corl[(n2),0] - corl[(n1),0] )**2 + ( corl[(n2),1] - corl[(n1),1] )**2 )**0.5)

#print(ln(2,3))
def cos(node1,node2,type):
    if node1>node2:
        node1,node2 = node2,node1

    return ((type[node2,0]-type[node1,0])/ln(node1,node2,type))

#print(cos(1,2))
def sin(node1, node2,type):
    if node1>node2:
        node1,node2 = node2,node1
    return ((type[node2,1]-type[node1,1])/ln(node1,node2,type))
dcord = np.array([i for i in cord])#np.copy(cord)

strain = np.array([((ln(*i,dcord) - ln(*i,cord))/(ln(*i,cord))) for i in m])
q = const*strain

trans_mats = np.array([[cos(*b,dcord),sin(*b,dcord),-cos(*b,dcord),-sin(*b,dcord)] for b in m])
#print(trans_mats)
#print(q)
g_nodal_force = np.array([trans_mats[c]*q[c] for c in range(len(m))])

int_force_mat=force_matrix(0,[0,0])
for f,val in enumerate(g_nodal_force):
    co = m[f] 
    if co[0]>co[1]:
        co[0],co[1] = co[1],co[0]
    force_mat += force_matrix(co[0],val[:2])
    force_mat+=force_matrix(co[1],val[2:])

force_mat -= int_force_mat
print(force_mat)
#print(np.array(stress_matrix())*a)
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
        
        #print(del_lst,'DEL LIST')
        np.set_printoptions(precision=2)
        return matrix

bc_force_matrix = boundary_condtions(force_mat,lst)

def iteration():
    global dcord
    global force_mat
    def local_stiff_mat(element,p):
        l= cos(*element, dcord)
        m= sin(*element, dcord)
        matrix= np.array([[  l**2 ,  l*m  , -l**2 , -l*m  ],
                        [  l*m  ,  m**2 , -l*m  , -m**2 ],
                        [ -l**2 , -l*m  ,  l**2 ,  l*m  ],
                        [ -l*m  , -m**2 ,  l*m  ,  m**2 ]]
                    )+  q[p]*np.array([[-m**2 , l*m , m**2 , -l*m],
                                    [l*m , -l**2 , -l*m , l**2],
                                    [m**2 , -l*m , -m**2 , l*m],
                                    [-l*m , l**2 , l*m , -l**2]])

        #print(matrix)
        return matrix/ln(*element,cord)
    #print(local_stiff_mat(d))
    np.set_printoptions(linewidth=np.inf)
    def global_stiff_mat(coord_list):
        matrix= np.zeros([len(cord)*2,len(cord)*2])
        global e_len
        e_len= []
        #print(matrix)
        for p,i in enumerate(coord_list):
            t1= np.array([[i[0],i[0]], [i[0],i[1]], [i[1],i[0]], [i[1],i[1]]])
            t2= np.array([[0,0], [0,2], [2,0], [2,2]])
            #print(t1)
            k= local_stiff_mat(i,p)
            for j in range(4):
                matrix[2*t1[j,0]:2*t1[j,0]+2 , 2*t1[j,1]:2*t1[j,1]+2] += k[t2[j,0]:t2[j,0]+2 , t2[j,1]:t2[j,1]+2]
            e_len.append(ln(*i,cord))
            #print(k*youngs_modulus*a)
        np.set_printoptions(precision=2)
        return matrix

    #print(global_stiff_mat(m)*const)


    #print(force_mat)


    #AFTER MULTIPYING WITH CONSTANT
    global_stiff_matrix= global_stiff_mat(m)*const


    #print(global_stiff_matrix)
    #APPLYING BOUNDARY CONDITIONS
    bc_stiff_mat= boundary_condtions(global_stiff_matrix,lst)
    #print(bc_stiff_mat,'bc stiff')
    bc_force_matrix= boundary_condtions(force_mat,lst)
    #print(bc_force_matrix, 'bc force')

    #INVERTING STIFFNESS MATRIX
    bc_stiff_mat_inv= np.linalg.inv(bc_stiff_mat)

    #MULTIPLYING INVERSE TO FORCE MATRIX
    nodal_disp= bc_stiff_mat_inv@bc_force_matrix
    #print(nodal_disp)
    for i in del_lst:
        nodal_disp= np.insert(nodal_disp,i,0)
    nodal_disp1 = nodal_disp.copy()
    dcord += nodal_disp1.reshape(len(cord), 2)

    np.set_printoptions(precision=4)
    #print(nodal_disp, 'UNDEFORMED')


    def stress_matrix():
        
        global def_e_len
        
        def_e_len= []
        for i in m:
            def_e_len.append(ln(*i,cord))
        #print(def_e_len)
        global strain
        strain= []
        for j in range(len(e_len)):
            strain.append((def_e_len[j]-e_len[j])/e_len[j])
        stress= [youngs_modulus*k for k in strain]
        #np.set_printoptions(precision= 2)
        return stress


    print(g_nodal_force)

    print(bc_force_matrix)
    force_mat = bc_force_matrix
    #print(global_stiff_mat(m)*const)
    bc_stiff_mat= boundary_condtions(global_stiff_matrix,lst)
    #print(bc_stiff_mat/10**3)
    def color_plotter(coordinate,strs):
        j=0
        for k in m:
            k=[coordinate[k[0]],coordinate[k[1]]]
            c=[k[0][0],k[1][0]]
            d=[k[0][1],k[1][1]]
            
            fos= strs[j]/yield_stress
            #print(fos)
            if fos**2>=1:
                r,g,b,a = 0 ,1 ,0 ,1
                plt.annotate('material will fail at this load',(3000,-5000),color= 'white')
                plt.plot(c,d,color=(r,g,b,a))
            elif fos<=0:
                r,g,b,a = 1+fos, 1+fos, 1,  1
                plt.plot(c,d,color=(r,g,b,a))
            elif fos>0:
                r,g,b,a= 1, 1-fos, 1-fos,  1         
                plt.plot(c,d,color=(r,g,b,a))
            j+=1
        
        #plt.ylim(-6000,6000)
        #plt.xlim(-1000,12000)
        #plt.xlim(cord[0,0]*1.2,cord[-1,0]*1.2)
        
        #plt.show()

    plt.axes().set_facecolor('black')

    color_plotter(dcord, stress_matrix())

for i in range(3):
    iteration()






""" a = 1.0
f = 4.8
l=5
f=0

si = lambda e:20*(e-e**2)
et = lambda e:20*(1-2*e)
kt = a*et(0)*l
r = 2.4
du = r/kt
de = (1/l)*du
e = 
j = np.linspace(0,9,30)
#plt.plot(j,[si(i) for i in j])
#plt.show()
 """