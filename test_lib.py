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
from Console import *
from System import *
from Keyboard_Joystick import *
from tess import Container
import MDAnalysis as mda


def set_camera(key):
    phi, the, r = key.get_coord()
    glTranslatef(0,0 ,r)
    glRotatef( the,1, 0, 0)
    glRotatef( phi,0, 0, 1)

        
s     = System("topol.tpr","traj_comp.xtc")
#s.get_hbond_mol()
s.get_network_hbond()
key_B = Keyboard_Joystick(s,-113, 310, -32.5)
cons  = Console((1200, 800),50.0)
cons.init()
'''
sel = []

for i, mol in enumerate(s.all_molecule):
    pos = mol.center_mass
    if np.sum(pos**2)<50:
        sel.append(i)

sel = np.array(sel)
s.set_visible_mol(sel)
'''
count = 0
while True:
    cons.clear()

    glPushMatrix()

    set_camera(key_B)
    if count % 170 ==0:
        #s.next()
        pass
    count+=1
    s.show_network_hbonds()
    #s.show_system()
    glPopMatrix()

    cons.update()
