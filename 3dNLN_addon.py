bl_info = {
    "name": "Non_linear_Truss",
    "author": "K@rth!kN",
    "version": (1, 0),
    "blender": (2, 80, 0),
    "location": "3D Viewport > Sidebar > 3DTrussFEA",
    "description": "Easy NON-LINEAR 3D Truss simulation in blender!",
    "category": "Development",
}

import bpy
import math
from bpy.props import FloatProperty, IntProperty, BoolProperty, StringProperty
import subprocess
import os
import csv
import time
script_dir = os.path.dirname(os.path.abspath(__file__))
dir_path = os.path.dirname(script_dir)
output_bytes = subprocess.check_output('where python', shell=True)

output_str = output_bytes.decode('utf-8')
python_paths = output_str.strip().split('\n')
main_path = os.path.join(dir_path, 'Non_linear_3D_truss.py')
if python_paths:
    python_interpreter = python_paths[0].strip()

#verts = bpy.context.active_object.data.vertices
#edges = bpy.context.active_object.data.edges

def define_mesh():
    so = bpy.context.active_object
    if so and so.type == 'MESH' :
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='SELECT')
        bpy.ops.mesh.delete(type = 'ONLY_FACE')
        bpy.ops.object.mode_set(mode='OBJECT')
    ###################################
    bpy.context.view_layer.update()
    global verts
    verts = so.data.vertices
    global edges
    edges = so.data.edges
    global faces
    faces = so.data.polygons
    
    coord_path = os.path.join(dir_path, 'truss_coord.csv')
    #print(script_dir)
    with open(coord_path,'w') as cord:
        try:
            for i in verts:
                global_cord = so.matrix_world @ i.co
                cord.write(f'{global_cord.x},{global_cord.y},{global_cord.z}\n')
        except Exception:
            print('failed')
    
    con_path = os.path.join(dir_path, 'truss_connect.csv')
    with open(con_path,'w') as connect:
        for j in edges:
            connect.write(f'{j.vertices[0]},{j.vertices[1]}\n')

#print([[i.vertices[0] , i.vertices[1]] for i in bpy.context.active_object.data.edges])
        

    
def call_sparser():
    scene = bpy.context.scene
    ip = scene.my_addon_properties.interp_path
    if ip == '':
        subprocess.run([python_interpreter, main_path, 'sparser'  ], check = True)
    else:   
        subprocess.run([ip, main_path, 'sparser'  ], check = True)

def call_solver():
    scene = bpy.context.scene
    ip = scene.my_addon_properties.interp_path
    if ip == '':
        subprocess.run([python_interpreter, main_path, 'main_sequence'  ], check = True)
    else:   
        subprocess.run([ip, main_path, 'main_sequence'  ], check = True)
        
class MyAddonProperties(bpy.types.PropertyGroup):
    opt_x: bpy.props.BoolProperty(name='UX', default=False)
    opt_y: bpy.props.BoolProperty(name='UY', default=False)
    opt_z: bpy.props.BoolProperty(name='UZ', default=False)
    

    
    force_X: FloatProperty(
    name="FX",
    default=0.0,
    description="load vector in x direction"
    )
    force_Y: FloatProperty(
    name="FY",
    default=0.0,
    description="load vector in y direction"
    )
    force_Z: FloatProperty(
    name="FZ",
    default=0.0,
    description="load vector in z direction"
    )

    youngs_modulus : FloatProperty(
    name = 'Youngs Modulus',
    default = 0.0,
    description = 'Enter youngs modulus of material in N/m^2',
    min = 0
    )
    distribute : BoolProperty(
        name = 'Distribute load to selection',
        default = False
    )

    
    area : FloatProperty(
    name = 'Section Area(m^2)',
    default = 0.0,
    description = 'Enter area of the cross section in m^2',
    min = 0
    )
    

    
    interp_path : StringProperty(
    name = 'Interpreter Path',
    description = 'bpy does not support matplotlib. Enter path of different python interpreter with matplotlib installed',
    default = ''
    )
    
    yield_stress : FloatProperty(
    name = 'Yield Stress',
    default = 0.0,
    description = 'Enter yield stress of material in N/m**2',
    min = 0,
    )

    convergence: FloatProperty(
    name = 'Floating Point Convergence',
    default = 0.001,
    description = 'Enter tolerance',
    min = 0,
    max = 0.1
    )

    minimum_opalescence: FloatProperty(
    name = 'Minimum Plot Opalescence',
    default = 0.2,
    description = 'minimum opalescece value of the plots from 0 for lowest to 1 for highest',
    min = 0,
    max = 1
    )
    
    plot_enum : bpy.props.EnumProperty(
        name = 'Plot library',
        description = 'select the python library to be used for plotting',
        items = [
            ('OP1' , 'Matplotlib' , 'Standard. Takes less time'),
            ('OP2' , 'Mayavi' , 'Looks better. Takes longer time'),])
            
class MESH_OT_properties(bpy.types.Operator):
    bl_idname = 'mesh.properties'
    bl_label = 'add properties'
    
    def execute(self,context):
        scene = bpy.context.scene
        y = scene.my_addon_properties.youngs_modulus
        a = scene.my_addon_properties.area
        ys = scene.my_addon_properties.yield_stress
        conv = scene.my_addon_properties.convergence
        op = scene.my_addon_properties.minimum_opalescence
        plop = scene.my_addon_properties.plot_enum
        #slfwt = scene.my_addon_properties.slfwt

        #if slfwt == False:
            #de = 0.0
        pl = 0
        if plop == 'OP1':
            pl = 0
        else: 
            pl = 1
        prop_path = os.path.join(dir_path, 'truss_properties.csv')
        with open(prop_path, 'w') as prop:
            prop.write(f'{y},{a},{ys},{conv},{op},{pl}')
        return {'FINISHED'}
    


    
class MESH_OT_mesh_creator(bpy.types.Operator):
    """adds the vertex and edge data into the list"""
    bl_idname = 'mesh.meshdata_add'
    bl_label = 'Apply Mesh'
    
    @classmethod
    def poll(cls, context):
        return context.mode == 'OBJECT'
    def execute(self,context):
        t1 = time.time()
        bpy.ops.object.mode_set(mode='EDIT')
        #bpy.ops.mesh.sort_elements(type='RANDOMIZE')
        bpy.ops.object.mode_set(mode='OBJECT')
        methods = ['VIEW_ZAXIS', 'VIEW_XAXIS', 'CURSOR_DISTANCE']
        print('Node numbering Methods: ',methods)
        def rnge(method):
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.mesh.sort_elements(type = method)
            edges = bpy.context.active_object.data.edges
            bpy.ops.object.mode_set(mode='OBJECT')
            max_arr = [abs(i.vertices[0]-i.vertices[1]) for i in edges ]
            #print(max(max_arr)*6+6)
            return max(max_arr)*3+3
        
    
        minima = [rnge(pl) for pl in methods]
        print('Possible Bandwidths: ',minima)
        print('Method for Minimum Bandwidth: ', methods[minima.index(min(minima))], f':{min(minima)}' )
        print('')
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.sort_elements(type=methods[minima.index(min(minima))])
        bpy.ops.object.mode_set(mode='OBJECT')
        define_mesh()
        t2 = time.time()
        #print(t2-t1)
        return {'FINISHED'}
    


#CREATING BOUNDARY CONDITION OPERATOR
class MESH_OT_add_bc_fr_vertex_list(bpy.types.Operator): # CATEGORY_TYPE_NAME
    """adds the selected vertices to the boundary conditions list. needs the current active object """
    
    
    bl_idname = 'mesh.bc_fr_vertices_add' 
    bl_label  = 'add boundary conditions'
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        return context.mode == 'OBJECT'
    def execute(self,context):
        scene = bpy.context.scene
            
        #global bound_list
        #bound_list =[]
        so = bpy.context.active_object
        verts = so.data.vertices
        x = scene.my_addon_properties.opt_x
        y = scene.my_addon_properties.opt_y
        z = scene.my_addon_properties.opt_z


        
        #all = scene.my_addon_properties.opt_xyz 
        bc_select_vert = [list(verts).index(v) for v in verts if v.select]
 
        bc_path = os.path.join(dir_path, 'truss_bound_con.csv')
    
        with open(bc_path, 'a+') as bc:
            def letter(letter_to_find,vert):
                if os.path.getsize(bc_path) == 0 :
                    return False
                else:
                    with open(bc_path, 'r') as bc:
                        for row in csv.reader(bc) :
                            if letter_to_find in row and str(vert) in row:
                                return True
                    return False
                
            for i in bc_select_vert:                                

                if  x and letter('x',i)==False:
                        
                        bc.write(f'{i},x\n')      

                if y and letter('y',i)==False:
                        
                        bc.write(f'{i},y\n')
                                          
                if z and letter('z',i)==False:
                        
                        bc.write(f'{i},z\n')
                        

                        

        return {'FINISHED'}

    
  
    
class MESH_OT_reset_fr_bc_list(bpy.types.Operator):
    bl_idname = 'mesh.reset_fr_bound_con'
    bl_label = 'reset to none'
    
    def execute(self,context):
        bc_path = os.path.join(dir_path, 'truss_bound_con.csv')
        with open(bc_path, 'w') as clr:
           clr.write('') 

        return {'FINISHED'}
       


# SELECTING THE VERTICES FOR FORCE APPLICATION
class MESH_OT_fr_forces(bpy.types.Operator):
    """Applies the given forces on the selected vertices"""
    bl_idname = 'mesh.fr_load_vector_add'
    bl_label = 'add forces'
    bl_options = {'REGISTER', 'UNDO'}
    
    @classmethod
    def poll(cls, context):
        return context.mode == 'OBJECT'
    
    def execute(self,context):
        scene = bpy.context.scene
        verts = bpy.context.active_object.data.vertices
        f1 = scene.my_addon_properties.force_X
        f2 = scene.my_addon_properties.force_Y
        f3 = scene.my_addon_properties.force_Z
        dist = scene.my_addon_properties.distribute

        force_select_vert = [list(verts).index(v) for v in verts if v.select]
        
        force_path = os.path.join(dir_path, 'truss_force_list.csv')
        dist = scene.my_addon_properties.distribute
        force_select_vert = [list(verts).index(v) for v in verts if v.select]
        


        with open(force_path, 'a+') as force:     
            no = len(force_select_vert)       
            nums = [f1,f2,f3]
            if dist :
                nums = [f/no for f in nums]
            for i in force_select_vert:
                nums.append(i)
            csv.writer(force).writerow(nums)
        return {'FINISHED'}
    
    
class MESH_OT_reset_fr_force_list(bpy.types.Operator):
    bl_idname = 'mesh.reset_fr_forces'
    bl_label = 'reset to none'
    
    def execute(self,context):
        force_path = os.path.join(dir_path, 'truss_force_list.csv')
        with open(force_path, 'w') as clr:
           clr.write('') 
        return {'FINISHED'}
    
    
    # Add the path to the directory containing your module
    #sys.path.append('C:\\Users\\karth\\Documents\\CODE')
    
class MESH_OT_solver(bpy.types.Operator):
    import sys
    bl_idname = 'mesh.solver'
    bl_label = 'Show'
    bl_options = {'REGISTER', 'UNDO'}
    def execute(self,context):
        call_solver()
        return {'FINISHED'}
    

class MESH_OT_sparser(bpy.types.Operator):
    import sys
    bl_idname = 'mesh.sparser'
    bl_label = 'run_code'
    bl_options = {'REGISTER', 'UNDO'}
    def execute(self,context):
        call_sparser()
        return {'FINISHED'}

 
#MASTER PARENT CLASS
class VIEW3D_PT_FEA_panel_1(bpy.types.Panel): 
    bl_space_type = 'VIEW_3D'  #3d viewport
    bl_region_type = 'UI'   #sidebar region
    
    #ADDING LABELS
    bl_label = 'PRE-PROCESSING'
    
    #ADDING THE CATEGORY
    bl_category = 'TrussFEA'
    
    def draw(self,context):
        pass

class VIEW3D_PT_preprocessor(bpy.types.Panel):
    bl_space_type = 'VIEW_3D'  #3d viewport
    bl_region_type = 'UI'   #sidebar region
    bl_parent_id = "VIEW3D_PT_FEA_panel_1"
    #ADDING LABELS
    bl_label ='MATERIAL PROERTIES'
    
    #ADDING THE CATEGORY
    bl_category = 'TrussFEA'
    
    def draw(self,context):
        scene = context.scene
        addon_props = scene.my_addon_properties
        
        box = self.layout.box()
        
        row = box.row()
        row.label(text = 'Material Properties')
        
        row = box.row()
        box.prop(addon_props, 'youngs_modulus')

        row = box.row()
        row.prop(addon_props, 'area')
        
        row = box.row()
        row.prop(addon_props, 'yield_stress')
        
        row = box.row()
        row.label(text = 'General')
        row = box.row()
        row.prop(addon_props, 'convergence')       
    
        row = box.row()
        row.prop(addon_props, 'minimum_opalescence') 
        
        row= box.row()
        row.prop(addon_props, 'plot_enum')
        
        row = box.row()
        row.operator('mesh.properties', text = 'UPDATE')

    


class VIEW3D_PT_mesher(bpy.types.Panel):
    bl_space_type = 'VIEW_3D'  #3d viewport
    bl_region_type = 'UI'   #sidebar region
    bl_parent_id = "VIEW3D_PT_FEA_panel_1"
    #ADDING LABELS
    bl_label = 'MESHING'
    
    #ADDING THE CATEGORY
    bl_category = 'TrussFEA'
    
    def draw(self,context):
        scene = context.scene
        addon_props = scene.my_addon_properties

        
        row = self.layout.row()
        row.operator('mesh.meshdata_add' , text = 'APPLY')
        
        layout = self.layout

        ob = context.object
        group = ob.vertex_groups.active

        rows = 3
        if group:
            rows = 5

        row = layout.row()
        row.template_list("MESH_UL_vgroups", "", ob, "vertex_groups", ob.vertex_groups, "active_index", rows=rows)

        col = row.column(align=True)

        col.operator("object.vertex_group_add", icon='ADD', text="")
        props = col.operator("object.vertex_group_remove", icon='REMOVE', text="")
        props.all_unlocked = props.all = False

        col.separator()

        col.menu("MESH_MT_vertex_group_context_menu", icon='DOWNARROW_HLT', text="")

        if group:
            col.separator()
            col.operator("object.vertex_group_move", icon='TRIA_UP', text="").direction = 'UP'
            col.operator("object.vertex_group_move", icon='TRIA_DOWN', text="").direction = 'DOWN'

        if (
                ob.vertex_groups and
                (ob.mode == 'EDIT' or
                 (ob.mode == 'WEIGHT_PAINT' and ob.type == 'MESH' and ob.data.use_paint_mask_vertex))
        ):
            row = layout.row()

            sub = row.row(align=True)
            sub.operator("object.vertex_group_assign", text="Assign")
            sub.operator("object.vertex_group_remove_from", text="Remove")

            sub = row.row(align=True)
            sub.operator("object.vertex_group_select", text="Select")
            sub.operator("object.vertex_group_deselect", text="Deselect")



class VIEW3D_PT_conditions(bpy.types.Panel):
    bl_space_type = 'VIEW_3D'  #3d viewport
    bl_region_type = 'UI'   #sidebar region
    #bl_parent_id = "VIEW3D_PT_FEA_panel_1"
    #ADDING LABELS
    bl_label ='LOADS'
    
    #ADDING THE CATEGORY
    bl_category = 'TrussFEA'
    
    def draw(self,context):
        pass
        
        
class VIEW3D_PT_boundcon(bpy.types.Panel):
    bl_category = 'TrussFEA'
    bl_space_type = 'VIEW_3D'  #3d viewport
    bl_region_type = 'UI'
    #bl_context = ".objectmode"  # dot on purpose (access from topbar)
    bl_label = "BOUNDARY CONDITIONS"
    bl_parent_id = "VIEW3D_PT_conditions"

    def draw(self, context):
        scene = context.scene
        addon_props = scene.my_addon_properties
        
        ########################
        layout = self.layout

        col = self.layout.column()
        col.label(text='BOUNDARY CONDITIONS')
        
        row = self.layout.row()
        row.prop(addon_props, "opt_x")

        
        row = self.layout.row()
        row.prop(addon_props, "opt_y")

        
        row = self.layout.row()
        row.prop(addon_props, "opt_z")


        col = self.layout.column()

        
        
        col.operator('mesh.bc_fr_vertices_add', text = 'Add Boundary Conditions')
        col.operator('mesh.reset_fr_bound_con', text = 'Reset Boundary Conditions')
        
        bc_path = os.path.join(dir_path, 'truss_bound_con.csv')
        with open(bc_path, 'r') as fb:
            l1 = []
            for i in csv.reader(fb):
                l1.append(i)
            ln1 = len(l1)
            col.label(text = f'Active Boundary Conditions: {ln1}')
            
            
class VIEW3D_PT_forces(bpy.types.Panel):
    bl_category = 'TrussFEA'
    bl_space_type = 'VIEW_3D'  #3d viewport
    bl_region_type = 'UI'
    #bl_context = ".objectmode"  # dot on purpose (access from topbar)
    bl_label = "FORCES"
    bl_parent_id = "VIEW3D_PT_conditions"
    def draw(self,context):
        scene = context.scene
        addon_props = scene.my_addon_properties
        
        ########################
        col = self.layout.column()
        col.separator()
        col.separator()
        
        col = self.layout.column()
        col.label(text = 'FORCES')
        
        row = self.layout.row()
        row.prop(addon_props, "force_X")

        
        row = self.layout.row()
        row.prop(addon_props, "force_Y")

        
        row = self.layout.row()
        row.prop(addon_props, "force_Z")
        
        row = self.layout.row()
        row.prop(addon_props, 'distribute')
        
        col = self.layout.column()
        col.operator('mesh.fr_load_vector_add', text = 'Add Forces and Moments')
        col.operator('mesh.reset_fr_forces', text = 'Reset Forces')
        
        force_path = os.path.join(dir_path, 'truss_force_list.csv')
        with open(force_path, 'r') as ff:
            l2 = []
            for j in csv.reader(ff):
                l2.append(j)
            ln2 = len(l2)
            col.label(text = f'Active Forces: {int(ln2/2)}')    
            
            
class VIEW3D_PT_solverpanel(bpy.types.Panel):
    bl_space_type = 'VIEW_3D'  #3d viewport
    bl_region_type = 'UI'   #sidebar region
    #bl_parent_id = "VIEW3D_PT_FEA_panel_1"
    #ADDING LABELS
    bl_label ='SOLUTION'
    
    #ADDING THE CATEGORY
    bl_category = 'TrussFEA'
    
    def draw(self,context):  
        scene = context.scene
        addon_props = scene.my_addon_properties


        col = self.layout.column()
        col.separator()
        
        row = self.layout.row()
        row.prop(addon_props, 'interp_path')

        row = self.layout.row()
        row.operator('mesh.solver', text = 'SOLVE and SHOW!')

        
        row = self.layout.row()
        row.operator('mesh.sparser', text = 'Show Sparsity')
  
classes = [VIEW3D_PT_FEA_panel_1, VIEW3D_PT_preprocessor ,VIEW3D_PT_mesher,MESH_OT_add_bc_fr_vertex_list, MESH_OT_mesh_creator, MESH_OT_fr_forces,VIEW3D_PT_conditions,VIEW3D_PT_boundcon ,VIEW3D_PT_forces,VIEW3D_PT_solverpanel,MESH_OT_solver,MESH_OT_reset_fr_bc_list,MESH_OT_reset_fr_force_list, MyAddonProperties, MESH_OT_properties, MESH_OT_sparser ]

#REGISTER THE PANEL CLASS IN BLENDER
def register():
    for cls in classes:
        bpy.utils.register_class(cls)
        
    bpy.types.Scene.my_addon_properties = bpy.props.PointerProperty(type=MyAddonProperties)
  
#UNREGISTER THE PANEL CLASS IN BLENDER
def unregister():

    for clss in classes:            
        bpy.utils.unregister_class(clss)
        
    del bpy.types.Scene.my_addon_properties
        
if __name__ == '__main__':
    register()

    


