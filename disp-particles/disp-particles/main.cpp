
#include <GL/glew.h>
#include <GL/freeglut.h>

#include <paramgl.h>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <assert.h>
#include <math.h>
#include <cutil_inline.h>
#include <cutil_gl_inline.h>
#include <shrUtils.h>

#include "bodysystemcuda.h"
#include "render_particles.h"
#include "cuda_runtime_api.h"

ParticleRenderer::DisplayMode displayMode = ParticleRenderer::PARTICLE_SPRITES_COLOR;

ParticleRenderer *m_renderer;

BodySystem<float>     *m_nbody;
BodySystemCUDA<float> *m_nbodyCuda;

float* m_hPos;
float* m_hVel;
float* m_hColor;

int buttonState = 0;
int ox = 0, oy = 0;

// view params
float camera_trans[]     = {0, -2, -100};
float camera_rot[]       = {0, 0, 0};
float camera_trans_lag[] = {0, -2, -100};
float camera_rot_lag[]   = {0, 0, 0};
const float inertia      = 0.1f;

int p = 256;
int q = 1;
bool useHostMem = false;
int  numDevsRequested = 1;
bool displayEnabled = true;
bool bPause = false;
bool bFullscreen = false;
bool bDispInteractions = false;
bool bSupportDouble = false;
int flopsPerInteraction = 20;

char deviceName[100];

enum { M_VIEW = 0, M_MOVE };

int numBodies = 16384;

int numIterations = 0; // run until exit

// fps
static int fpsCount = 0;
static int fpsLimit = 5;
cudaEvent_t startEvent, stopEvent;
cudaEvent_t hostMemSyncEvent;

struct NBodyParams
{       
    float m_timestep;
    float m_clusterScale;
    float m_velocityScale;
    float m_softening;
    float m_damping;
    float m_pointSize;
    float m_x, m_y, m_z;
};

NBodyParams activeParams = { 0.016f, 1.54f, 8.0f, 0.1f, 1.0f, 1.0f, 0, -2, -100}

void computePerfStats(double &interactionsPerSecond, double &gflops, 
                      float milliseconds, int iterations)
{
    // double precision uses intrinsic operation followed by refinement,
    // resulting in higher operation count per interaction.
    // (Note Astrophysicists use 38 flops per interaction no matter what,
    // based on "historical precident", but they are using FLOP/s as a 
    // measure of "science throughput". We are using it as a measure of 
    // hardware throughput.  They should really use interactions/s...
    // const int flopsPerInteraction = fp64 ? 30 : 20; 
    interactionsPerSecond = (float)numBodies * (float)numBodies;
    interactionsPerSecond *= 1e-9 * iterations * 1000 / milliseconds;
    gflops = interactionsPerSecond * (float)flopsPerInteraction;
}

// check for OpenGL errors
void checkGLErrors(const char *s)
{
    GLenum error;
    while ((error = glGetError()) != GL_NO_ERROR) {
        fprintf(stderr, "%s: error - %s\n", s, (char *) gluErrorString(error));
    }
}

void initGL(int *argc, char** argv)
{  
    // First initialize OpenGL context, so we can properly set the GL for CUDA.
    // This is necessary in order to achieve optimal performance with OpenGL/CUDA interop.
    glutInit(argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize(720, 480);
    glutCreateWindow("CUDA n-body system");
    if (bFullscreen)
        glutFullScreen();
    GLenum err = glewInit();

    if (GLEW_OK != err)
    {
        shrLog("GLEW Error: %s\n", glewGetErrorString(err));
        cutilDeviceReset();
        exit(-1);
    }
    else if (!glewIsSupported("GL_VERSION_2_0 "
                         "GL_VERSION_1_5 "
			             "GL_ARB_multitexture "
                         "GL_ARB_vertex_buffer_object")) 
    {
        fprintf(stderr, "Required OpenGL extensions missing.");
        exit(-1);
    }
    else
    {
#if   defined(WIN32)
        wglSwapIntervalEXT(0);
#elif defined(LINUX)
        glxSwapIntervalSGI(0);
#endif      
    }

    glEnable(GL_DEPTH_TEST);
    glClearColor(0.0, 0.0, 0.0, 1.0);
   
    checkGLErrors("initGL");
}

void displayNBodySystem()
{
    m_renderer->setSpriteSize(activeParams.m_pointSize);
	m_renderer->setPBO(m_nbody->getCurrentReadBuffer(), m_nbody->getNumBodies(), (sizeof(float) > 4));
	m_renderer->display(displayMode);
}

void _resetRenderer()
{
    float color[4] = { 1.0f, 0.6f, 0.3f, 1.0f};
    m_renderer->setBaseColor(color);

    m_renderer->setColors(m_hColor, m_nbody->getNumBodies());
    m_renderer->setSpriteSize(activeParams.m_pointSize);  
}

void display()
{
    static double gflops = 0;
    static double ifps = 0;
    static double interactionsPerSecond = 0;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);  
    if (displayEnabled)
    {
        // view transform for camera
        {
            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();

            for (int c = 0; c < 3; ++c)
            {
                camera_trans_lag[c] += (camera_trans[c] - camera_trans_lag[c]) * inertia;
                camera_rot_lag[c] += (camera_rot[c] - camera_rot_lag[c]) * inertia;
            }

            glTranslatef(camera_trans_lag[0], camera_trans_lag[1], camera_trans_lag[2]);
            glRotatef(camera_rot_lag[0], 1.0, 0.0, 0.0);
            glRotatef(camera_rot_lag[1], 0.0, 1.0, 0.0);
        }

        displayNBodySystem();

        if (bFullscreen)
        {
            beginWinCoords();
            char msg0[256], msg1[256], msg2[256];

            if (bDispInteractions) {
               sprintf(msg1, "%0.2f billion interactions per second", interactionsPerSecond);
            } else {
               sprintf(msg1, "%0.2f GFLOP/s", gflops);
            }
            sprintf(msg0, "%s", deviceName);
            sprintf(msg2, "%0.2f FPS [%s | %d bodies]", 
                    ifps, "single precision", numBodies);  

            glBlendFunc(GL_ONE_MINUS_DST_COLOR, GL_ZERO); // invert color
            glEnable(GL_BLEND);
            glColor3f(0.46f, 0.73f, 0.0f);
            glPrint(80, glutGet(GLUT_WINDOW_HEIGHT) - 122, msg0, GLUT_BITMAP_TIMES_ROMAN_24);
            glColor3f(1.0f, 1.0f, 1.0f);
            glPrint(80, glutGet(GLUT_WINDOW_HEIGHT) - 96, msg2, GLUT_BITMAP_TIMES_ROMAN_24);
            glColor3f(1.0f, 1.0f, 1.0f);
            glPrint(80, glutGet(GLUT_WINDOW_HEIGHT) - 70, msg1, GLUT_BITMAP_TIMES_ROMAN_24);
            glDisable(GL_BLEND);

            endWinCoords();          
        }

        glutSwapBuffers();
    }

    fpsCount++;
    // this displays the frame rate updated every second (independent of frame rate)
    if (fpsCount >= fpsLimit) 
    {
        char fps[256];
    
        float milliseconds = 1;
        // stop timer
        cutilSafeCall(cudaEventRecord(stopEvent, 0));  
        cutilSafeCall(cudaEventSynchronize(stopEvent));
        cutilSafeCall( cudaEventElapsedTime(&milliseconds, startEvent, stopEvent));
        
        milliseconds /= (float)fpsCount;
        computePerfStats(interactionsPerSecond, gflops, milliseconds, 1);

        ifps = 1.f / (milliseconds / 1000.f);
        sprintf(fps, 
		        "CUDA N-Body (%d bodies): "
		        "%0.1f fps | %0.1f BIPS | %0.1f GFLOP/s | %s", 
		        numBodies, ifps, interactionsPerSecond, gflops,
                "single precision");  

        glutSetWindowTitle(fps);
        fpsCount = 0; 
        fpsLimit = (ifps > 1.f) ? (int)ifps : 1;
        if (bPause) fpsLimit = 0;
        
        // restart timer
        cutilSafeCall(cudaEventRecord(startEvent, 0));
    }

    glutReportErrors();
}

void reshape(int w, int h)
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0, (float) w / (float) h, 0.1, 1000.0);

    glMatrixMode(GL_MODELVIEW);
    glViewport(0, 0, w, h);
}

void mouse(int button, int state, int x, int y)
{
    int mods;

    if (state == GLUT_DOWN)
        buttonState |= 1<<button;
    else if (state == GLUT_UP)
        buttonState = 0;

    mods = glutGetModifiers();
    if (mods & GLUT_ACTIVE_SHIFT) 
    {
        buttonState = 2;
    } 
    else if (mods & GLUT_ACTIVE_CTRL) 
    {
        buttonState = 3;
    }

    ox = x; oy = y;

    glutPostRedisplay();
}

void motion(int x, int y)
{
    float dx = (float)(x - ox);
    float dy = (float)(y - oy);

    if (buttonState == 3) 
    {
        // left+middle = zoom
        camera_trans[2] += (dy / 100.0f) * 0.5f * fabs(camera_trans[2]);
    } 
    else if (buttonState & 2) 
    {
        // middle = translate
        camera_trans[0] += dx / 100.0f;
        camera_trans[1] -= dy / 100.0f;
    }
    else if (buttonState & 1) 
    {
        // left = rotate
        camera_rot[0] += dy / 5.0f;
        camera_rot[1] += dx / 5.0f;
    }
    
    ox = x; oy = y;
    glutPostRedisplay();
}

void idle(void)
{
    glutPostRedisplay();
}

//////////////////////////////////////////////////////////////////////////////
// Program main
//////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv) 
{
	initGL(&argc, argv);

	m_nbodyCuda = new BodySystemCUDA<float>(numBodies, 1, p, q, 1, 0);
    m_nbody = m_nbodyCuda;

    // allocate host memory
    m_hPos = new float[numBodies*4];
    m_hVel = new float[numBodies*4];
    m_hColor = new float[numBodies*4];

    m_nbody->setSoftening(activeParams.m_softening);
    m_nbody->setDamping(activeParams.m_damping);
		
    cutilSafeCall( cudaEventCreate(&startEvent) );
    cutilSafeCall( cudaEventCreate(&stopEvent) );
    cutilSafeCall( cudaEventCreate(&hostMemSyncEvent) );

    m_renderer = new ParticleRenderer;
    _resetRenderer();
	
	glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutIdleFunc(idle);

	glutMainLoop();
}