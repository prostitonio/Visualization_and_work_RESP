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



class System:
    all_molecule=[]
    coord_visible = True
    index_selection = 0
    flag_selection = False
    def __init__(self,*all_mol):
        for mol in all_mol:
            if type(mol) == Molecula:
                self.all_molecule.append(mol)
            else:
                pass

    def show_Sphere(self,xyz,r,color):
        glPushMatrix()
        glColor3f(*color)
        quad = gluNewQuadric()
        glTranslatef(*xyz);
        grid = 18
        gluSphere(quad,r,grid,grid)
        glTranslatef(*(-xyz));
        glPopMatrix()

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

    def append(self,mol):
        if type(mol) == Molecula or type(mol) == Molecula_test:
            self.all_molecule.append(mol)
        else:
            #print("add")
            pass

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

    def show_system(self):
        #print("show sys")
        if self.coord_visible == True:
            self.show_sys_coord() 
        #print(self.all_molecule)
        for mol in self.all_molecule:
            mol.show_molecule()
