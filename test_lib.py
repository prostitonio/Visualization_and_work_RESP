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
    global phi
    global the
    global r
    glTranslatef(0,0 ,r)
    glRotatef( the,1, 0, 0)
    glRotatef( phi,0, 0, 1)
def get_pressed():
    global phi
    global the
    global dthe
    global r
    global dr
    global joysticks
    global joystick_on_off

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
class Molecula:
    conv = 0.002
    #conv = 0.01
    all_r             = {"O":152*conv,"H":110*conv,"C":170*conv        ,"N":155*conv}
    all_color         = {"O":(1,0,0) ,"H":(1,1,1) ,"C":(0.22,0.22,0.22),"N":(0.22,0.16,0.88)}
    all_type_mass     = {"O":16      ,"H":1       ,"C":12              ,"N":14}
    def __init__(self,pdb_name,log_name="",del_mol_proc=0):
        pdb = PDBFile(pdb_name)
        self.all_atom     = 10*np.array([[i.x,i.y,i.z] for i in pdb.positions       ])
        self.all_type     = [                  i.name[0]for i in pdb.topology.atoms()]
        self.all_mass     = np.array([ self.all_type_mass[i] for i in self.all_type])
        self.all_constr   = np.array([[i.atom1.index,i.atom2.index]  for i in pdb.topology.bonds()])
        self.center_mass  = self.get_center_mass() 


        if log_name !="":
            self.esp_on = log_name !=""
            esp_coord ,esp_ch = self.esp_coord_ch(log_name)
            self.esp_coord,self.esp_ch = self.del_exes(esp_coord,esp_ch,del_mol_proc)
            self.esp_color_ch = self.ch_to_color(self.esp_ch)
            self.update_esp()

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

    def ch_to_color(self,esp_ch):
        a = (esp_ch - min(esp_ch))/(max(esp_ch)-min(esp_ch))
        d = np.zeros((len(a),3))+1
        d[:,0] = a
        esp_color = []
        for i in d:
            esp_color.append(colorsys.hsv_to_rgb(*i))
        return np.array(esp_color )

    def show_molecule(self):
        for i,count in enumerate(self.all_atom):
            self.show_Sphere(count,self.all_r[self.all_type[i]],self.all_color[self.all_type[i]])
        for i in self.all_constr:
            self.show_constr(i,self.all_atom,(1,1,1))

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




mol = Molecula('UNK_D65243.pdb',"UNK_C7DFF1.log",0.0)
mol.set_normal_coord()
#mol1 = Molecula('UNK_F4033B.pdb')
#mol1.shift_coord([0,2,2])


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
 
