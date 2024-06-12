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
import time 

class Keyboard_Joystick:
    phi = -137
    the = 310
    r = -5

    dr_mouse = 0.5
    dphi_mouse = 3
    
    dr_joysticks = 0.1
    dphi_joysticks = 1.2

    mode = "rotate_main_system"
    joysticks = {}
    joystick_on_off = False
    


    def __init__(self,sys,phi=-137,the=310,r=-5):
        self.phi = phi
        self.the = the
        self.r   = r
        self.sys = sys

    def set_mouse_speed(self,dphi,dr):
        self.dr_mouse = dr
        self.dphi_mouse = dphi 

    def set_joystick_speed(self,dphi,dr):
        self.dr_joysticks = dr
        self.dphi_joysticks = dphi 

    def get_coord(self):
        self.set_coord()
        return self.phi, self.the, self.r

    def set_coord(self):
        all_event = pygame.event.get()
        for event in all_event:
            if event.type == pygame.QUIT:
                pygame.quit()
                quit()

            if event.type == pygame.JOYDEVICEADDED:
                self.joystick_on_off = True
                joy = pygame.joystick.Joystick(event.device_index)
                self.joysticks[joy.get_instance_id()] = joy
                print(f"Joystick {joy.get_instance_id()} connencted")

            if event.type == pygame.JOYDEVICEREMOVED:
                self.joystick_on_off = False
                del joysticks[event.instance_id]
                print(f"Joystick {event.instance_id} disconnected")

            if event.type == pygame.KEYDOWN:
                #print(event.key)
                if event.key == 115: # S
                    if self.mode != "selection_molecules":
                        self.sys.start_selection()
                        self.mode = "selection_molecules"
                    else:
                        self.sys.end_selection()
                        self.mode = "rotate_main_system"
        if self.joystick_on_off == True:
            for event in all_event:
                if event.type == pygame.JOYBUTTONDOWN:
                    joystick_count = pygame.joystick.get_count()
                    #text_print.tprint(screen, f"Number of joysticks: {joystick_count}")
                    for joystick in self.joysticks.values():
                        button = joystick.get_button(0)
    #                    time.sleep(0.1)
                        if button == 1:
                            if self.mode != "selection_molecules":
                                self.sys.start_selection()
                                self.mode = "selection_molecules"
                            else:
                                self.sys.end_selection()
                                self.mode = "rotate_main_system"

        match self.mode:
            case "rotate_main_system":
                self.rotate_main_system(all_event)
            case "selection_molecules":
                self.selection_molecules(all_event)
    def selection_molecules(self,all_event):
        for event in all_event:
            if event.type == pygame.KEYDOWN:
                if event.key == 110: # N
                    self.sys.next_selection()
                if event.key == 98: # B
                    self.sys.previous_selection()
                

    def rotate_main_system(self,all_event):
        if self.joystick_on_off == False:
            
            keys = pygame.key.get_pressed()
            if keys[pygame.K_LEFT]:
                self.phi = self.phi + self.dphi_mouse
            if keys[pygame.K_RIGHT]:
                self.phi = self.phi - self.dphi_mouse
            if keys[pygame.K_UP]:
                self.the = self.the + self.dphi_mouse
            if keys[pygame.K_DOWN]:
                self.the = self.the - self.dphi_mouse

            for event in all_event:
                if event.type == pygame.MOUSEBUTTONDOWN:
                    if event.button == 4:
                        self.r=self.r + self.dr_mouse
                    if event.button == 5:
                        self.r=self.r -self.dr_mouse

        else:
            joystick_count = pygame.joystick.get_count()
            #text_print.tprint(screen, f"Number of joysticks: {joystick_count}")
            for joystick in self.joysticks.values():
                jid = joystick.get_instance_id()
                #text_print.tprint(screen, f"Joystick {jid}")
                name = joystick.get_name()
                #text_print.tprint(screen, f"Joystick name: {name}")
                guid = joystick.get_guid()
                #text_print.tprint(screen, f"GUID: {guid}")
                power_level = joystick.get_power_level()
                #text_print.tprint(screen, f"Joystick's power level: {power_level}")
                axes = joystick.get_numaxes()
                #text_print.tprint(screen, f"Number of axes: {axes}")
                self.r   =  self.r - self.dr_joysticks    *round(joystick.get_axis(1),2)
                self.phi =  self.phi+ self.dphi_joysticks *round(joystick.get_axis(2),2)
                self.the = self.the + self.dphi_joysticks *round(joystick.get_axis(3),2)
                #for i in range(axes):
                #    axis = joystick.get_axis(i)
                    #text_print.tprint(screen, f"Axis {i} value: {axis:>6.3f}")
                buttons = joystick.get_numbuttons()
                #text_print.tprint(screen, f"Number of buttons: {buttons}")
                for i in range(buttons):
                    button = joystick.get_button(i)
                    #text_print.tprint(screen, f"Button {i:>2} value: {button}")
                hats = joystick.get_numhats()
                #text_print.tprint(screen, f"Number of hats: {hats}")
                for i in range(hats):
                    hat = joystick.get_hat(i)
                    #text_print.tprint(screen, f"Hat {i} value: {str(hat)}")
