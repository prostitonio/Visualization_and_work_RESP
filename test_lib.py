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
from tess import Container

def save_com_file(name,sys):
    in_file = ""
    in_file += "%mem=6GB\n"
    in_file += "%NProcShared=32\n"
    in_file += "# HF/6-31G* SCF=Tight Pop=MK IOp(6/33=2,6/41=10,6/42=17) Nosymm\n"
    in_file += "\ncomment\n\n"
    in_file += "0 1\n"
    for mol in sys.all_molecule:
        for i in range(len(mol.all_atom)):
            in_file+=mol.all_type[i]+" "
            for j in mol.all_atom[i]:
                in_file+=str(j)+" "
            in_file+="\n"
    in_file+="\n"
    with open(name, 'w') as f:
        f.write(in_file)
    #print(in_file) 


def set_camera(key):
    phi, the, r = key.get_coord()
    glTranslatef(0,0 ,r)
    glRotatef( the,1, 0, 0)
    glRotatef( phi,0, 0, 1)

        
mol = Molecula('UNK_D65243.pdb')
#mol4 = Molecula('UNK_D65243.pdb',"test.log",0.0)
#mol3 = Molecula('UNK_D65243.pdb')
#mol2 = Molecula('UNK_D65243.pdb')
#mol4.visible_mol = False
mol.shift_coord( [3,0.55,2.2])
#mol2.shift_coord([0,3,0])
#mol3.shift_coord([3,0,0])
sys = System()

sys.append(mol)
#sys.append(mol2)
#sys.append(mol3)
#sys.append(mol4)

cntr = Container(mol.all_atom, limits=(6,3,4.1), periodic=True)
ver  = np.array([v.face_vertices() for v in cntr])
ver_coord = np.array([v.vertices() for v in cntr])
def gen_color():
    return np.random.rand(3)
def show_constr(a1,a2,color):
    glLineWidth(2.)
    glBegin(GL_LINES)
    glColor3f(*color);glVertex3d(*a1);glVertex3d(*a2)
    glEnd()
def show_couter(all_atom,all_num,color):
    for place in all_num:
        for i in range(len(place)-1):
            show_constr(all_atom[place[i]],all_atom[place[i+1]],color)


#sys.start_selection()
#save_com_file("test.com",sys)
key_B = Keyboard_Joystick(sys)

pygame.init()
display = (800, 600)
pygame.display.set_mode(display, DOUBLEBUF|OPENGL)

glEnable(GL_DEPTH_TEST)
glLoadIdentity()
gluPerspective(90, (display[0]/display[1]), 0.1, 50.0)
glClearColor(0.01,0.01,0.1,0) 

while True:
    np.random.seed(1)
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
    glPushMatrix()

    set_camera(key_B)
    sys.show_system()
    for i in range(len(ver)):
        show_couter(ver_coord[i],ver[i],mol.all_color[i])

    
    glPopMatrix()

    #get_pressed()
    pygame.display.flip()
    pygame.time.wait(10)

