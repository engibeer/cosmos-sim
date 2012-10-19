from OpenGL.GL import *
from OpenGL.arrays import vbo
from OpenGL.GL import shaders
from OpenGL.GLUT import *
from OpenGL.GLU import *

import numpy
from numpy import array

class PointRenderer():

    def __init__(self):
        self.m_vbo = None

        self.vert = """
        void main(){
            gl_FrontColor = gl_Color;
            gl_Position = ftransform();
        }
        """

        # green
        self.frag = """
        void main(){
            gl_FragColor = vec4( 0, 1, 0, 1 );
        }
        """

        self.shader_program = None

        self.bl = 0

        # view params
        self.ox = 0
        self.oy = 0
        self.buttonState      = 0
        self.camera_trans     = [0, -2, -25]
        self.camera_rot       = [0, 0, 0]
        self.camera_trans_lag = [0, -2, -25]
        self.camera_rot_lag   = [0, 0, 0]
        self.inertia          = 0.1

    def define_shaders(self):
        self.shader_program = glCreateProgram()

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

    def reshape(self, width, height):
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(60.0, width / height, 0.1, 1000.0)
        glMatrixMode(GL_MODELVIEW)
        glViewport(0, 0, width, height)

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

    def motion(self,x, y):
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
        
        self.ox = x
        self.oy = y
        glutPostRedisplay()

    def display(self):      
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT )
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        camera_trans_lag = self.camera_trans_lag
        c = 0
        for c in xrange(3):
            self.camera_trans_lag[c] += (self.camera_trans[c] - camera_trans_lag[c]) * self.inertia;
            self.camera_rot_lag[c] += (self.camera_rot[c] - self.camera_rot_lag[c]) * self.inertia;
        glTranslatef(camera_trans_lag[0], camera_trans_lag[1], camera_trans_lag[2])
        glRotatef(self.camera_rot_lag[0], 1.0, 0.0, 0.0)
        glRotatef(self.camera_rot_lag[1], 0.0, 1.0, 0.0)
        glUseProgram(self.shader_program)

        self.m_vbo.bind()

        glEnableClientState( GL_VERTEX_ARRAY )
        glVertexPointerf(self.m_vbo)
        glDrawArrays(GL_TRIANGLES, 0, len(self.m_vbo))

        glutSwapBuffers()
        glutReportErrors()

    def idle(self):
        self.m_vbo = vbo.VBO(
            array([
                [  0, 1, 0 ], # triangle
                [ -1,-1, 0 ],
                [  1,-1, 0 ],
                [  2,-1, 0 ], # square
                [  4,-1, 0 ],
                [  4, 1, 0 ],
                [  2,-1, 0 ],
                [  4, 1, 0 ],
                [  2, 1, 0 ],
            ],'f')
        )
        glutPostRedisplay()

    def begin(self):
        glutInit()
        glutInitWindowSize(400, 400)
        glutCreateWindow("Simple OpenGL")
        glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB)
        glutDisplayFunc(self.display)
        glutReshapeFunc(self.reshape)
        glutMouseFunc(self.mouse)
        glutMotionFunc(self.motion)
        glutIdleFunc(self.idle)

        self.define_shaders()

        glutMainLoop()

if __name__ == '__main__':
    pr = PointRenderer()
    pr.begin()