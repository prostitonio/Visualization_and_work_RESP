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
from Molecule import *
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from numba import jit
import networkx as nx


class System:
    coord_visible = True
    index_selection = 0
    flag_selection = False
    time_step = 0
    
    def __init__(self,conf,trj):
        self.u            = mda.Universe(conf,trj)
        self.res_num      = self.u.atoms.resnums-1
        self.hbonds = HBA(universe=self.u)

        self.uniq         = np.unique(self.res_num)
        self.num_mol      = len(self.uniq)
        self.num_atom     = len(self.res_num)
        self.visible_mol  = np.arange(self.num_mol)
        self.visible_atom = np.arange(self.num_atom)
        self.size_box     = self.u.universe.dimensions[:3]
        self.constrains   = []
        for i in self.u.atoms.bonds:
            #self.constrains.append([i.atoms[0].id-1,i.atoms[1].id-1])
            self.constrains.append([i.atoms[0].id,i.atoms[1].id])
        self.constrains = np.array(self.constrains)

        self.visible_constrain = self.constrains.copy() 
        self.set_frame()

    def get_network_hbond(self):
        G = nx.Graph()
        G.add_nodes_from(np.arange(self.num_mol))
        G.add_edges_from(self.hb_mol)
        t = nx.connected_components(G)
        self.all_claster_hbond = [list(c) for c in t ]

    def get_hbond_mol(self):
        self.hbonds.run(start=self.time_step,stop=self.time_step+1)
        hb_atoms = np.array(self.hbonds.results.hbonds[:,2:4],dtype=int)
        #print(hb_atoms)
        hb_mol   = [] 
        for i,j in hb_atoms:
            hb_mol.append([self.res_num[i],self.res_num[j]])
        self.hb_mol = np.array(hb_mol)
        #print(hb_mol)


    def set_frame(self):
        self.all_molecule=[]
        self.frame = self.u.trajectory[self.time_step]
        for i in range(len(self.uniq)):
            sel = np.where(self.res_num==self.uniq[i])[0]
            pos =  self.frame.positions[sel]
            name = self.u.atoms.names[sel]
            name = np.array([i[0] for i in name])
            m: Molecula = Molecula(all_atom_coord=pos,all_atom_type=name,all_atom_constr=np.array([]))
            self.all_molecule.append(m)

        self.update()
        self.shift_coord(-self.get_center_mass_sistem())
        self.update()

        self.get_hbond_mol()
        self.get_network_hbond()

    def next(self):
        self.time_step = self.time_step + 1
        self.set_frame()

    def get_center_mass_sistem(self):
        return np.einsum('ij,i',self.all_atom_c_mass,self.all_atom_mass)/np.sum(self.all_atom_mass)

    def shift_coord(self,xyz):
        xyz = np.array(xyz)
        for mol in self.all_molecule:
            mol.all_atom_coord = mol.all_atom_coord + xyz
            mol.center_mass    = mol.center_mass    + xyz

    def update(self):
        self.all_atom_coord = self.all_molecule[0].all_atom_coord
        self.all_atom_r     = self.all_molecule[0].all_atom_r
        self.all_atom_color = self.all_molecule[0].all_atom_color
        self.all_atom_c_mass= np.array([self.all_molecule[0].center_mass])
        self.all_atom_mass  = np.array([self.all_molecule[0].mass_mol])
        for mol in self.all_molecule[1:]:
            self.all_atom_coord = np.append(self.all_atom_coord ,mol.all_atom_coord,axis=0)
            self.all_atom_r     = np.append(self.all_atom_r     ,mol.all_atom_r    ,axis=0)
            self.all_atom_color = np.append(self.all_atom_color ,mol.all_atom_color,axis=0)
            self.all_atom_c_mass= np.append(self.all_atom_c_mass,np.array([mol.center_mass]),axis=0)
            self.all_atom_mass  = np.append(self.all_atom_mass  ,np.array([mol.mass_mol   ]),axis=0)

    def append(self,mol):
        if type(mol) == Molecula: 
            self.all_molecule.append(mol)
        else:
            pass

    def show_Sphere(self,xyz,r,color):
        glPushMatrix()
        glColor3f(*color)
        quad = gluNewQuadric()
        glTranslatef(*xyz);
        grid = 7
        gluSphere(quad,r,grid,grid)
        glTranslatef(*(-xyz));
        glPopMatrix()

    def show_system(self):
        if self.coord_visible == True:
            self.show_sys_coord()
        for i in self.visible_atom:
            self.show_Sphere(self.all_atom_coord[i],self.all_atom_r[i],self.all_atom_color[i])
        for i in self.visible_constrain:
            self.show_constr(i,self.all_atom_coord,(1,1,1))

    def show_network_hbonds(self):
        np.random.seed(1)
        if self.coord_visible == True:
            self.show_sys_coord()
        for clas in self.all_claster_hbond:
            color = np.random.rand(3)
            if len(clas)<10:
                for mol in clas:
                    self.show_Sphere(self.all_atom_c_mass[mol],0.4,color)
        for i in self.hb_mol:
            #self.show_constr_hbons(i,self.all_atom_c_mass,(0.9,0.1,0))
            pass

    def show_constr(self,a,all_coord,color):
        a1,a2 = a
        #print(np.sum((all_coord[a1]-all_coord[a2])**2))
        glLineWidth(4)
        glBegin(GL_LINES)
        glColor3f(*color);glVertex3d(*all_coord[a1]);glVertex3d(*all_coord[a2])
        glEnd()
    def show_constr_hbons(self,a,all_coord,color):
        a1,a2 = a
        #print(np.sum((all_coord[a1]-all_coord[a2])**2))
        if np.sum((all_coord[a1]-all_coord[a2])**2) < 30:
            glLineWidth(1)
            glBegin(GL_LINES)
            glColor3f(*color);glVertex3d(*all_coord[a1]);glVertex3d(*all_coord[a2])
            glEnd()
        else:
            pass 

    def hide_random_molecules(self,num):
        temp       = self.visible_mol.copy()
        np.random.shuffle(temp)
        select_mol = np.sort(temp[:num])
        self.set_visible_mol(select_mol)

    def set_visible_all(self):
        self.visible_mol       = np.arange(self.num_mol)
        self.visible_atom      = np.arange(self.num_atom)
        self.visible_constrain = self.constrains

    def set_visible_mol(self,select_mol):
        self.visible_mol = select_mol
        a = np.array([])
        for i in select_mol:
            a = np.append(np.where(self.res_num==i)[0],a)
        self.visible_atom = np.array(a,dtype=int)
        test_vis = list(self.visible_atom)
        vis_const = []
        for i in self.constrains:
            if i[0] in test_vis:
                vis_const.append(i)
        self.visible_constrain = np.array(vis_const)


    def show_sys_coord(self):
        _x = np.array([1,0,0])
        _y = np.array([0,1,0])
        _z = np.array([0,0,1])
        self.show_Sphere(_x,0.05,_x)
        self.show_Sphere(_y,0.05,_y)
        self.show_Sphere(_z,0.05,_z)
        glPushMatrix()
        glLineWidth(1.)
        glBegin(GL_LINES)
        glColor3f(*_x);glVertex3d(0,0,0);glVertex3d(*_x)
        glColor3f(*_y);glVertex3d(0,0,0);glVertex3d(*_y)
        glColor3f(*_z);glVertex3d(0,0,0);glVertex3d(*_z)
        glEnd()
        glPopMatrix()

'''
    def show_system(self):
        #print("show sys")
        if self.coord_visible == True:
            self.show_sys_coord() 
        #print(self.all_molecule)
        for mol in self.all_molecule:
            mol.show_molecule()
    def start_selection(self):
        flag = True
        for mol in self.all_molecule:
            mol.pick_out_mol(flag)
        self.all_molecule[self.index_selection].pick_out_mol(False)
    def next_selection(self):
        self.all_molecule[self.index_selection].pick_out_mol(True)
        self.index_selection = self.index_selection + 1
        if self.index_selection > len(self.all_molecule)-1:
            self.index_selection = 0
        self.all_molecule[self.index_selection].pick_out_mol(False)
    def previous_selection(self):
        self.all_molecule[self.index_selection].pick_out_mol(True)
        self.index_selection = self.index_selection - 1
        if self.index_selection < 0:
            self.index_selection = len(self.all_molecule)-1
        self.all_molecule[self.index_selection].pick_out_mol(False)

    def end_selection(self):
        flag = False
        for mol in self.all_molecule:
            mol.pick_out_mol(flag)


'''
