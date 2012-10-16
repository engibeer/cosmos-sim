from OpenGL.GL import *
from OpenGL.arrays import vbo
from OpenGL.GL import shaders
from OpenGL.GLUT import *
from OpenGL.GLU import *

from numpy import array

vert = """
void main(){
    gl_FrontColor = gl_Color;
    gl_Position = ftransform();
}
"""

frag = """
void main(){
    gl_FragColor = vec4( 0, 1, 0, 1 );
}
"""

shader_program = None

def define_shader():
    global shader_program
    shader_program = glCreateProgram()

    vobj = glCreateShader( GL_VERTEX_SHADER )
    glShaderSource( vobj, vert )
    glCompileShader( vobj )
    print glGetShaderInfoLog(vobj)
    glAttachShader( shader_program, vobj )

    fobj = glCreateShader( GL_FRAGMENT_SHADER )
    glShaderSource( fobj, frag)    
    glCompileShader( fobj )
    print glGetShaderInfoLog(fobj)
    glAttachShader( shader_program, fobj )

    glLinkProgram( shader_program )

def reshape(width, height):
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(65.0, width / float(height), 1, 1000 );
    glMatrixMode(GL_MODELVIEW);

def display():
    m_vbo = vbo.VBO(
            array( [
                [  0, 1, 0 ],
                [ -1,-1, 0 ],
                [  1,-1, 0 ],
                [  2,-1, 0 ],
                [  4,-1, 0 ],
                [  4, 1, 0 ],
                [  2,-1, 0 ],
                [  4, 1, 0 ],
                [  2, 1, 0 ],
            ],'f')
        )
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT )  
    glLoadIdentity()   
    gluLookAt( 60.0,60.0,25.0,
              0.0,0.0,0.0,
              0.0,0.0,1.0 )

    glTranslatef(25.0, 25.0, 0)

    glUseProgram(shader_program)

    m_vbo.bind()

    glEnableClientState( GL_VERTEX_ARRAY )
    glVertexPointerf(m_vbo)
    glDrawArrays(GL_TRIANGLES, 0, 9)

    glutSwapBuffers()

if __name__ == '__main__':
    glutInit()
    glutInitWindowSize(400,400)
    glutCreateWindow("Simple OpenGL")
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB)
    glutDisplayFunc(display)
    glutReshapeFunc(reshape)
    
    define_shader()

    glutMainLoop()