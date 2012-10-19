from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
from OpenGL import GL as gl
import random
import time
import numpy
import math

#These three defines exist in OpenGL.GL, but does not correspond to those used here
GL_GEOMETRY_SHADER_EXT       = 0x8DD9
GL_GEOMETRY_INPUT_TYPE_EXT   = 0x8DDB
GL_GEOMETRY_OUTPUT_TYPE_EXT  = 0x8DDC
GL_GEOMETRY_VERTICES_OUT_EXT = 0x8DDA

vert = '''
void main(){
    gl_FrontColor = gl_Color;
    gl_Position = ftransform();
}
'''

geom = '''#version 120 
#extension GL_EXT_geometry_shader4 : enable 
const float radius = 0.5;
varying out vec2 coord;

void main() 
{
  gl_FrontColor = gl_FrontColorIn[0];

  coord = vec2( -1,-1 );
  gl_Position = (gl_PositionIn[0] + gl_ProjectionMatrix * vec4(-radius,-radius,0,0) );
  EmitVertex();
  coord = vec2( -1,1 );
  gl_Position = (gl_PositionIn[0] + gl_ProjectionMatrix * vec4(-radius,radius, 0,0) );
  EmitVertex();
  coord = vec2( 1,-1 );
  gl_Position = (gl_PositionIn[0] + gl_ProjectionMatrix * vec4( radius,-radius, 0,0) );
  EmitVertex();
  coord = vec2( 1,1 );
  gl_Position = (gl_PositionIn[0] + gl_ProjectionMatrix * vec4( radius,radius, 0,0) );
  EmitVertex();  
  EndPrimitive();
}
'''
frag = '''
varying in vec2 coord;
void main(){
    vec4 color = gl_Color;
    color.a = 1.0-length( coord );
    if (color.a<0.0) discard;

    // A VERY FAKE "lighting" model
    float d = dot( normalize(vec3(coord,1.0)), vec3( 0.19, 0.19, 0.96225 ) );
    color.rgb *= d*d;
    // end "lighting"
    
    gl_FragColor = color;
}
'''

theta = 0.0
delta = 0.0
TIMEOUTFACTOR = 5.0

# view params
ox = 0
oy = 0
buttonState      = 0
camera_trans     = [0, -2, -25]
camera_rot       = [0, 0, 0]
camera_trans_lag = [0, -2, -25]
camera_rot_lag   = [0, 0, 0]
inertia          = 0.1

def my_idle( ):
    global theta
    global delta 
    t = time.clock()
    passed = t-delta
    theta += passed
    delta = t

    dirvector = numpy.array( [ -1-(math.sin( theta/2.0 )), -1-(math.cos(theta/2.0)), 1.0 ] )
    dirvector *=  (1.0 / math.sqrt( (dirvector**2).sum() ) )
    

    global POINTS
    global VEL
    global TIMES

    #Add "gravity" to all the velocities
    VEL[:,2] -= passed * 9.81 * 5
    #Increment positions according to velocities
    POINTS += passed*VEL
    
    pbelow_zero = POINTS[:,2] < 0.0
    # "drag" all points below zero (in z) up again
    POINTS[pbelow_zero] *= [1,1,0]
    # Change the z velocity for all "collisions" with the zero z plane
    VEL[ pbelow_zero ] *= [1,1,-1.0]

    # add passed seconds to all points timeout
    TIMES += passed

    timeouts = TIMES > TIMEOUTFACTOR
    TIMES[timeouts] = numpy.random.random( timeouts.sum() ) * 0.1 * TIMEOUTFACTOR
    VEL[timeouts] = dirvector + 0.5*numpy.random.random( timeouts.sum() * 3 ).reshape( (-1,3) )
    VEL[timeouts] *= 50.0
    POINTS[timeouts] = [0.0,0.0,0.0]
    
    glutPostRedisplay()

shader_program = None
def define_shader():
    global shader_program
    shader_program = gl.glCreateProgram()   
    glProgramParameteri(shader_program, GL_GEOMETRY_INPUT_TYPE_EXT, gl.GL_POINTS )
    glProgramParameteri(shader_program, GL_GEOMETRY_OUTPUT_TYPE_EXT, gl.GL_TRIANGLE_STRIP )
    glProgramParameteri(shader_program, GL_GEOMETRY_VERTICES_OUT_EXT, 4 )

    vobj = glCreateShader( GL_VERTEX_SHADER )
    glShaderSource( vobj, vert )
    glCompileShader( vobj )
    print glGetShaderInfoLog(vobj)
    glAttachShader( shader_program, vobj )
    
    gobj = glCreateShader( GL_GEOMETRY_SHADER_EXT )
    glShaderSource( gobj, geom)    
    glCompileShader( gobj )
    print glGetShaderInfoLog(gobj)
    glAttachShader( shader_program, gobj )

    fobj = glCreateShader( GL_FRAGMENT_SHADER )
    glShaderSource( fobj, frag)    
    glCompileShader( fobj )
    print glGetShaderInfoLog(fobj)
    glAttachShader( shader_program, fobj )
    
    glLinkProgram( shader_program )
    print glGetProgramInfoLog(shader_program)

def reshape( width, height ):
   glViewport(0, 0, width, height);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   gluPerspective(65.0, width / float(height), 1, 1000 );
   glMatrixMode(GL_MODELVIEW);

def mouse(button, state, x, y):
    global buttonState
    global ox, oy
    if (state == GLUT_DOWN):
        buttonState |= 1 << button
    elif (state == GLUT_UP):
        buttonState = 0;

    mods = glutGetModifiers()
    if (mods & GLUT_ACTIVE_SHIFT):
        buttonState = 2
    elif (mods & GLUT_ACTIVE_CTRL):
        buttonState = 3;

    ox = x
    oy = y;

    glutPostRedisplay()

def motion(x, y):
    global buttonState
    global ox, oy
    global camera_trans
    global camera_rot
    dx = x - ox
    dy = y - oy

    if buttonState == 3:
        # left+middle = zoom
        camera_trans[2] += (dy / 100.0) * 0.5 * abs(camera_trans[2])
    elif buttonState & 2:
        # middle = translate
        camera_trans[0] += dx / 100.0
        camera_trans[1] -= dy / 100.0
    elif buttonState & 1:
        # left = rotate
        camera_rot[0] += dy / 5.0
        camera_rot[1] += dx / 5.0
    
    ox = x
    oy = y
    glutPostRedisplay()

POINTS = None
COLORS = None
VEL = None
TIMES = None
def display( ):   
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT )  
   glLoadIdentity()   
   gluLookAt( 60.0,60.0, 25.0,
              0.0,0.0,0.0,
              0.0,0.0,1.0 )

   c = 0
   for c in xrange(3):
     camera_trans_lag[c] += (camera_trans[c] - camera_trans_lag[c]) * inertia;
     camera_rot_lag[c] += (camera_rot[c] - camera_rot_lag[c]) * inertia;
   glTranslatef(camera_trans_lag[0], camera_trans_lag[1], camera_trans_lag[2])
   glRotatef(camera_rot_lag[0], 1.0, 0.0, 0.0)
   glRotatef(camera_rot_lag[1], 0.0, 1.0, 0.0)
   glBlendFunc( GL_SRC_ALPHA, GL_ONE )
   glUseProgram( shader_program )
   glEnableClientState( GL_COLOR_ARRAY )
   glEnableClientState( GL_VERTEX_ARRAY )
   glColorPointer( 3,GL_FLOAT,0,COLORS )
   glVertexPointer(3,GL_FLOAT,0,POINTS )
   glDrawArrays( 0,0, len( POINTS ) )
  
   glutSwapBuffers()


def init():
    glutInit(sys.argv)
    glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE )
    glutInitWindowSize(1280, 720)
    # Open a window
    glutCreateWindow("Glut test window")
    
    glutReshapeFunc( reshape )
    glutDisplayFunc( display )
    glutMouseFunc( mouse )
    glutMotionFunc( motion )
    glutIdleFunc( my_idle )
    glEnable( GL_DEPTH_TEST )
    glClearColor( 1,1,1,0 )

    define_shader()

    # 100k Balls
    count = 100000

    global POINTS
    global COLORS
    global VEL
    global TIMES    
    COLORS =numpy.random.random( count*3 ).reshape( (-1,3) )
    POINTS = 100000*numpy.ones( count*3 ).reshape( (-1,3) )
    VEL = numpy.zeros( count*3 ).reshape( (-1,3) )
    TIMES = numpy.random.random( count ) * TIMEOUTFACTOR
    glutMainLoop();

init()



