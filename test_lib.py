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
from System import *
from Keyboard_Joystick import *
np.random.seed(1)


def set_camera(key):
    phi, the, r = key.get_coord()
    glTranslatef(0,0 ,r)
    glRotatef( the,1, 0, 0)
    glRotatef( phi,0, 0, 1)

        
mol = Molecula('UNK_D65243.pdb')#,"UNK_C7DFF1.log",0.0)
mol3 = Molecula('UNK_D65243.pdb')#,"UNK_C7DFF1.log",0.0)
mol2 = Molecula('UNK_D65243.pdb')#,"UNK_C7DFF1.log",0.0)
mol.shift_coord([0,0,3])
mol2.shift_coord([0,3,0])
mol3.shift_coord([3,0,0])
sys = System()

sys.append(mol)
sys.append(mol2)
sys.append(mol3)
#sys.start_selection()

key_B = Keyboard_Joystick(sys)

pygame.init()
display = (800, 600)
pygame.display.set_mode(display, DOUBLEBUF|OPENGL)

glEnable(GL_DEPTH_TEST)
glLoadIdentity()
gluPerspective(90, (display[0]/display[1]), 0.1, 50.0)
glClearColor(0.01,0.01,0.1,0) 

while True:
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
    glPushMatrix()

    set_camera(key_B)
    sys.show_system()

    glPopMatrix()

    #get_pressed()
    pygame.display.flip()
    pygame.time.wait(10)

