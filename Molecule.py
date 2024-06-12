import pygame
import colorsys
import numpy as np
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
from pygame.locals import *
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
import ctypes
import json


np.random.seed(1)
#_______________________________________#
class Molecula:
    conv = 0.002
    #conv = 0.01
    with open('json_file/all_vdw_radius.json') as f:
        #all_r = json.load(f)
        all_r = json.loads(json.load(f))
    with open('json_file/all_color_type.json') as f:
        all_color =json.loads(json.load(f))
    all_type_mass     = {"O":16      ,"H":1       ,"C":12              ,"N":14}
    esp_on = False
    pick_out_flag = False
    

    def __init__(self,pdb_name,log_name="",del_mol_proc=0):
        pdb = PDBFile(pdb_name)
        self.all_atom     = 10*np.array([[i.x,i.y,i.z] for i in pdb.positions       ])
        self.all_type     = [                  i.name[0]for i in pdb.topology.atoms()]
        self.all_mass     = np.array([ self.all_type_mass[i] for i in self.all_type])
        self.all_constr   = np.array([[i.atom1.index,i.atom2.index]  for i in pdb.topology.bonds()])
        self.center_mass  = self.get_center_mass() 


        if log_name !="":
            #print("work")
            self.esp_on = log_name !=""
            esp_coord ,esp_ch = self.esp_coord_ch(log_name)
            self.esp_coord,self.esp_ch = self.del_exes(esp_coord,esp_ch,del_mol_proc)
            self.esp_color_ch = self.ch_to_color(self.esp_ch)
            self.update_esp()

    def pick_out_mol(self, flag):
        self.pick_out_flag = flag


    def get_center_mass(self):
        self.center_mass  = np.einsum('ij,i',self.all_atom,self.all_mass)/np.sum(self.all_mass)
        return self.center_mass

    def get_J(self):
        (x, y, z)    = self.all_atom.T
        x2, y2, z2 = x**2, y**2, z**2
        m          = self.all_mass
        Jx , Jy , Jz  = m*(y2+z2), m*(x2+z2), m*(x2+y2)
        Jxy, Jyz, Jxz = m*(-x*y) , m*(-y*z) , m*(-x*z)
        J = np.array([	[Jx, Jxy, Jxz],
        		[Jxy, Jy, Jyz],
        		[Jxz, Jyz, Jz]
        	    ])
        self.J = np.sum(J,axis=2)
        return self.J

    def rotate_to_abc(self):
        self.get_J()
        Q = np.linalg.eig(self.J)[1]
        self.all_atom = np.squeeze(np.matmul(self.all_atom,Q))
        self.esp_coord = np.squeeze(np.matmul(self.esp_coord,Q))
        self.update_esp()

    def del_exes(self,coord,ch,proc):
        a = np.arange(len(ch))
        np.random.shuffle(a)
        a = a[:int(len(ch)*(1-proc))]
        a.sort()
        return coord[a],ch[a]

    def esp_coord_ch(self,name):
        f = open(name)
        xyz = [];ch  = []
        for i in f:
            if "Fit" in i:
                sp = i.split()
                if sp[0]=="ESP":
                    xyz.append([float(j) for j in sp[6:9]])
                if sp[1]=="Fit" and len(sp)==3:
                    ch.append(float(sp[2]))
        return np.array(xyz),np.array(ch)

    def ch_to_color(self,esp_ch,pick_out_flag=False):
        a = (esp_ch - min(esp_ch))/(max(esp_ch)-min(esp_ch))
        d = np.zeros((len(a),3))+1
        d[:,0] = a
        esp_color = []
        for i in d:
            esp_color.append(colorsys.hsv_to_rgb(*i))
        return np.array(esp_color )

    def show_molecule(self,esp=False):
        dark = 0.3
        for i,count in enumerate(self.all_atom):
            if self.pick_out_flag == True:
                h , s , v  = colorsys.rgb_to_hsv(*self.all_color[self.all_type[i]])
                rgb = colorsys.hsv_to_rgb(h,s*dark,v*dark) 
                self.show_Sphere(count,self.all_r[self.all_type[i]]*self.conv,rgb)
            else:
                self.show_Sphere(count,self.all_r[self.all_type[i]]*self.conv,self.all_color[self.all_type[i]])
        for i in self.all_constr:
            if self.pick_out_flag == True:
                h , s , v  = colorsys.rgb_to_hsv(*(1,1,1))
                rgb = colorsys.hsv_to_rgb(h,s*dark,v*dark) 
                self.show_constr(i,self.all_atom,rgb)
            else:
                self.show_constr(i,self.all_atom,(1,1,1))
        if esp==True:
            self.show_ESP()

    def show_constr(self,a,all_coord,color):
        a1,a2 = a
        glLineWidth(7.)
        glBegin(GL_LINES)
        glColor3f(*color);glVertex3d(*all_coord[a1]);glVertex3d(*all_coord[a2])
        glEnd()
        
    def show_Sphere(self,xyz,r,color):
        glPushMatrix()
        glColor3f(*color)
        quad = gluNewQuadric()
        glTranslatef(*xyz);
        grid = 18
        gluSphere(quad,r,grid,grid)
        glTranslatef(*(-xyz));
        glPopMatrix()

    def show_ESP(self):
        glPointSize(2);
        glVertexPointer(3, GL_FLOAT, 12, self.buffer_offset(self.vertices.ctypes.data))
        glColorPointer( 3, GL_FLOAT, 12, self.buffer_offset(self.v_color.ctypes.data ))
        glEnableClientState(GL_VERTEX_ARRAY)
        glEnableClientState(GL_COLOR_ARRAY )
        glDrawArrays(GL_POINTS, 0, len(self.esp_coord))
        glDisableClientState(GL_VERTEX_ARRAY)
        glDisableClientState(GL_COLOR_ARRAY )

    def update_esp(self):
        coord = self.esp_coord.reshape(len(self.esp_coord)*3)
        color = self.esp_color_ch.reshape(len(self.esp_color_ch)*3)
        self.vertices = np.array(coord, dtype = np.float32)
        self.v_color  = np.array(color, dtype = np.float32)
        self.buffer_offset = ctypes.c_void_p
        self.float_size = ctypes.sizeof(ctypes.c_float)

    def shift_coord(self,xyz):
        xyz = np.array(xyz)
        self.all_atom = self.all_atom + xyz
        if self.esp_on:
            self.esp_coord = self.esp_coord + xyz
            self.update_esp()
    def set_normal_coord(self):
        self.shift_coord(-self.get_center_mass())
        self.rotate_to_abc()

#___________________________________________#




