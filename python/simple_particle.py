from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
from OpenGL import GL as gl
import random, time, numpy, math

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pymongo import Connection

class ParticleRenderer():

    def __init__(self):
        self.vert = '''
        void main(){
            gl_FrontColor = gl_Color;
            gl_Position = ftransform();
        }
        '''

        self.frag = """
        void main(){
            gl_FragColor = vec4( 1, 1, 1, 1 );
        }
        """

        self.theta = 0.0
        self.delta = 0.0
        self.TIMEOUTFACTOR = 60.0

        # view params
        self.ox = 0
        self.oy = 0
        self.buttonState      = 0
        self.camera_trans     = [0, -50, -25]
        self.camera_rot       = [5, 125, 0]
        self.camera_trans_lag = [0, -2, -25]
        self.camera_rot_lag   = [0, 0, 0]
        self.inertia          = 0.1

        self.max_distance = 5000

        self.points = None

        self.shader_program = None

    def idle(self):
        t = time.clock()
        passed = t - self.delta
        self.theta += passed
        self.delta = t  
        glutPostRedisplay()

    def define_shaders(self):
        self.shader_program = gl.glCreateProgram()   

        vobj = glCreateShader( GL_VERTEX_SHADER )
        glShaderSource( vobj, self.vert )
        glCompileShader( vobj )
        print glGetShaderInfoLog(vobj)
        glAttachShader( self.shader_program, vobj )

        fobj = glCreateShader( GL_FRAGMENT_SHADER )
        glShaderSource( fobj, self.frag)    
        glCompileShader( fobj )
        print glGetShaderInfoLog(fobj)
        glAttachShader( self.shader_program, fobj )
        
        glLinkProgram( self.shader_program )
        print glGetProgramInfoLog(self.shader_program)

    def reshape(self, width, height ):
        glViewport(0, 0, width, height);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluPerspective(5.0, width / float(height), 1, self.max_distance );
        glMatrixMode(GL_MODELVIEW);

    def mouse(self, button, state, x, y):
        if (state == GLUT_DOWN):
            self.buttonState |= 1 << button
        elif (state == GLUT_UP):
            self.buttonState = 0;

        mods = glutGetModifiers()
        if (mods & GLUT_ACTIVE_SHIFT):
            self.buttonState = 2
        elif (mods & GLUT_ACTIVE_CTRL):
            self.buttonState = 3;

        self.ox = x
        self.oy = y;

        glutPostRedisplay()

    def motion(self, x, y):
        dx = x - self.ox
        dy = y - self.oy

        if self.buttonState == 3:
            # left+middle = zoom
            self.camera_trans[2] += (dy / 100.0) * 0.5 * abs(self.camera_trans[2])
        elif self.buttonState & 2:
            # middle = translate
            self.camera_trans[0] += dx / 100.0
            self.camera_trans[1] -= dy / 100.0
        elif self.buttonState & 1:
            # left = rotate
            self.camera_rot[0] += dy / 5.0
            self.camera_rot[1] += dx / 5.0

        #print self.camera_rot[0]
        self.ox = x
        self.oy = y
        glutPostRedisplay()

    def display(self):   
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT )  
        glLoadIdentity()   

        c = 0
        for c in xrange(3):
          self.camera_trans_lag[c] += (self.camera_trans[c] - self.camera_trans_lag[c]) * self.inertia;
          self.camera_rot_lag[c] += (self.camera_rot[c] - self.camera_rot_lag[c]) * self.inertia;
        glTranslatef(self.camera_trans_lag[0], self.camera_trans_lag[1], self.camera_trans_lag[2])
        glRotatef(self.camera_rot_lag[0], 1.0, 0.0, 0.0)
        glRotatef(self.camera_rot_lag[1], 0.0, 1.0, 0.0)
        glBlendFunc( GL_SRC_ALPHA, GL_ONE )
        glUseProgram( self.shader_program )
        glEnableClientState( GL_VERTEX_ARRAY )
        glVertexPointer(3, GL_FLOAT, 0, self.points )
        glDrawArrays( 0, 0, len( self.points ) )
     
        glutSwapBuffers()

    def init(self):
        glutInit(sys.argv)
        glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE )
        glutInitWindowSize(1280, 720)
        # Open a window
        glutCreateWindow("Augury")
        
        glutReshapeFunc( self.reshape )
        glutDisplayFunc( self.display )
        glutMouseFunc( self.mouse )
        glutMotionFunc( self.motion )
        glutIdleFunc( self.idle )
        glEnable( GL_DEPTH_TEST )
        glClearColor(0, 0, 0, 0)

        self.define_shaders()

        connection = Connection('localhost', 27017)
        db = connection['test']
        collection = db['bigFOF']

        result = collection.find({ "snapnum" : 85 }, { 'x' : 1, 'y' : 1, 'z' : 1 })
        rows = result.count()
        print "rows: ", rows

        self.points = numpy.zeros((rows, 3))

        for i, row in enumerate(result):
            self.points[i][0] = float(row['x']) * 5.0
            self.points[i][1] = float(row['y']) * 5.0
            self.points[i][2] = float(row['z']) * 5.0

        glutMainLoop()

if __name__ == '__main__':
    pr = ParticleRenderer()
    pr.init()
