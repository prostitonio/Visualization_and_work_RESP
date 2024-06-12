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
np.random.seed(1)


phi = -137
the = 310
dphi = 3
r = -5
dr = 0.5
joysticks = {}
joystick_on_off = False

def show_Sphere(xyz,r,color):
    glPushMatrix()
    glColor3f(*color)
    quad = gluNewQuadric()
    glTranslatef(*xyz);
    grid = 18
    gluSphere(quad,r,grid,grid)
    glTranslatef(*(-xyz));
    glPopMatrix()
def show_sys_coord():
    _x = np.array([1,0,0])
    _y = np.array([0,1,0])
    _z = np.array([0,0,1])
    show_Sphere(_x,0.05,_x)
    show_Sphere(_y,0.05,_y)
    show_Sphere(_z,0.05,_z)
    glPushMatrix()
    glLineWidth(1.)
    glBegin(GL_LINES)
    glColor3f(*_x);glVertex3d(0,0,0);glVertex3d(*_x)
    glColor3f(*_y);glVertex3d(0,0,0);glVertex3d(*_y)
    glColor3f(*_z);glVertex3d(0,0,0);glVertex3d(*_z)
    glEnd()
    glPopMatrix()
def set_camera():
    global phi,the,r
    glTranslatef(0,0 ,r)
    glRotatef( the,1, 0, 0)
    glRotatef( phi,0, 0, 1)
def get_pressed():
    global phi,the,dthe,r,dr,joystick_on_off,joysticks

    keys = pygame.key.get_pressed()
    if joystick_on_off == False:
        if keys[pygame.K_LEFT]:
            phi = phi + dphi
        if keys[pygame.K_RIGHT]:
            phi = phi - dphi
        if keys[pygame.K_UP]:
            the = the + dphi
        if keys[pygame.K_DOWN]:
            the = the - dphi

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            pygame.quit()
            quit()

        if event.type == pygame.JOYDEVICEADDED:
            joystick_on_off = True
            joy = pygame.joystick.Joystick(event.device_index)
            joysticks[joy.get_instance_id()] = joy
            print(f"Joystick {joy.get_instance_id()} connencted")
        if event.type == pygame.JOYDEVICEREMOVED:
            joystick_on_off = False
            del joysticks[event.instance_id]
            print(f"Joystick {event.instance_id} disconnected")

        if joystick_on_off == False:
            if event.type == pygame.MOUSEBUTTONDOWN:
                if event.button == 4:
                    r+=dr
                if event.button == 5:
                    r-=dr
        else:
            if event.type == pygame.JOYBUTTONDOWN:
                #print("Joystick button pressed.")
                if event.button == 0:
                    joystick = joysticks[event.instance_id]
                    if joystick.rumble(0, 0.7, 500):
                        #print(f"Rumble effect played on joystick {event.instance_id}")
                        pass
            if event.type == pygame.JOYBUTTONUP:
                print("Joystick button released.")
            #text_print.reset()
    joystick_count = pygame.joystick.get_count()
    #text_print.tprint(screen, f"Number of joysticks: {joystick_count}")
    #text_print.indent()
    for joystick in joysticks.values():
        jid = joystick.get_instance_id()
        #text_print.tprint(screen, f"Joystick {jid}")
        #text_print.indent()
        name = joystick.get_name()
        #text_print.tprint(screen, f"Joystick name: {name}")
        guid = joystick.get_guid()
        #text_print.tprint(screen, f"GUID: {guid}")
        power_level = joystick.get_power_level()
        #text_print.tprint(screen, f"Joystick's power level: {power_level}")
        axes = joystick.get_numaxes()
        #text_print.tprint(screen, f"Number of axes: {axes}")
        #text_print.indent()
        d_r   = round(joystick.get_axis(1),2)
        d_phi = round(joystick.get_axis(2),2)
        d_the = round(joystick.get_axis(3),2)
        k_al = 0.7
        k_r  = 0.1
        the = the + k_al*d_the
        phi = phi + k_al*d_phi
        r   = r   - k_r*d_r
        #for i in range(axes):
        #    axis = joystick.get_axis(i)
            #text_print.tprint(screen, f"Axis {i} value: {axis:>6.3f}")
        #text_print.unindent()
        buttons = joystick.get_numbuttons()
        #text_print.tprint(screen, f"Number of buttons: {buttons}")
        #text_print.indent()
        for i in range(buttons):
            button = joystick.get_button(i)
            #text_print.tprint(screen, f"Button {i:>2} value: {button}")
        #text_print.unindent()
        hats = joystick.get_numhats()
        #text_print.tprint(screen, f"Number of hats: {hats}")
        #text_print.indent()
        for i in range(hats):
            hat = joystick.get_hat(i)
            #text_print.tprint(screen, f"Hat {i} value: {str(hat)}")
#_______________________________________#
mol = Molecula('UNK_D65243.pdb',"UNK_C7DFF1.log",0.0)
#mol.set_normal_coord()


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

    set_camera()
    show_sys_coord()
    mol.show_molecule()
#    mol1.show_molecule()
    mol.show_ESP()

    glPopMatrix()

    get_pressed()
    pygame.display.flip()
    pygame.time.wait(10)

