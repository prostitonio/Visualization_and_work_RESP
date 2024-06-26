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
import MDAnalysis as mda



class Console:
    def __init__(self,size_window,distanse):
        self.size_window = size_window
        self.distanse    = distanse

    def init(self):
        pygame.init()
        display = self.size_window
        pygame.display.set_mode(display, DOUBLEBUF|OPENGL)

        glEnable(GL_DEPTH_TEST)
        glLoadIdentity()
        gluPerspective(90, (display[0]/display[1]), 0.1, self.distanse)
        glClearColor(0.01,0.01,0.1,0) 
    def clear(self):
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
    
    def update(self):
        pygame.display.flip()
        pygame.time.wait(10)

