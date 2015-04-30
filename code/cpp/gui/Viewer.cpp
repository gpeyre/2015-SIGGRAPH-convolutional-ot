#include "Viewer.h"
#include "Image.h"

#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glext.h>
#else
#include <GL/gl.h>
#include <GL/glext.h>
#include <GL/glut.h>
#endif

int  Viewer::windowSize[2] = { 800, 800 };
bool Viewer::renderWireframe = false;
bool Viewer::renderSelected = false;

std::set<unsigned> Viewer::selectedVert;
ot::ConvSolver* Viewer::solverPtr = 0;
TriMesh*        Viewer::meshPtr = 0;
VectorXd        Viewer::values;
double Viewer::timestep = 0.0;
int    Viewer::verbose = 0;

GLuint Viewer::surfaceDL = 0;
Shader Viewer::shader;
Camera Viewer::camera;

void 
Viewer::keyboard(unsigned char c, int /*x*/, int /*y*/)
{
    switch(c)
    {
        case 27:
            exit(0);
            break;
        case 'w':
            renderWireframe = !renderWireframe;
            updateDisplayList();
            break;
        case 's':
            renderSelected = !renderSelected;
            updateDisplayList();
            break;
        case 'r':
            clearData();
            updateDisplayList();
            break;
        case 'g':
            runGeodesics();
            updateDisplayList();
            break;
        case 'i':
            runInterpolation();
            updateDisplayList();
            break;
        case '/':
            takeScreenshot();
            break;
        case '=':
            timestep += 0.1;
            break;
        case '-':
            timestep -= 0.1;
            break;
        default:
            break;
    }
}
   
void 
Viewer::init(int argc, char** argv)
{
    initGLUT(argc, argv);
    initGLSL();
    setGL();
    updateDisplayList();
    glutMainLoop();
}
   
void 
Viewer::initGLUT(int argc, char** argv)
{
    glutInitWindowSize(windowSize[0], windowSize[1]);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInit(&argc, argv);
    glutCreateWindow("heatOT");

    glutMouseFunc   (Viewer::mouse   );
    glutKeyboardFunc(Viewer::keyboard);
    glutSpecialFunc (Viewer::special );
    glutMotionFunc  (Viewer::motion  );
    glutDisplayFunc (Viewer::display );
    glutIdleFunc    (Viewer::idle    );
}
   
void 
Viewer::initGLSL()
{
    shader.loadVertex  ("../gui/vertex.glsl"  );
    shader.loadFragment("../gui/fragment.glsl");
}
   
void 
Viewer::special(int i, int /*x*/, int /*y*/)
{
    switch (i)
    {
        case GLUT_KEY_UP:
            camera.zoomIn();
            break;
        case GLUT_KEY_DOWN:
            camera.zoomOut();
            break;
        case 27:
            exit(0);
        break;
        default:
            break;
    }
}

void 
Viewer::mouse(int button, int state, int x, int y)
{
    if ((glutGetModifiers() & GLUT_ACTIVE_SHIFT) && (state == GLUT_UP))
    {
        pickVertex(x, y);
        return;
    }
    camera.mouse(button, state, x, y);
}

void 
Viewer::motion(int x, int y)
{
    camera.motion(x, y);
}

void 
Viewer::idle()
{
    glutPostRedisplay();
}

void 
Viewer::display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    shader.enable();

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    double aspect = double(viewport[2]) / double(viewport[3]);
    const double fovy = 50.;
    const double clipNear = .01;
    const double clipFar = 1000.;
    gluPerspective(fovy, aspect, clipNear, clipFar);
      
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    Camera::Quaternion    eye(0., 0., 0.,-2.5*camera.zoom());
    Camera::Quaternion center(0., 0., 0., 0. );
    Camera::Quaternion     up(0., 0., 1., 0. );
    gluLookAt(   eye.x(),    eye.y(),    eye.z(),
              center.x(), center.y(), center.z(),
                  up.x(),     up.y(),     up.z() );
            
    GLint uniformEye = glGetUniformLocation(shader, "eye");
    Camera::Quaternion r = camera.computeCurrentRotation();
    eye = r.conjugate() * eye * r;
    glUniform3f(uniformEye, eye.x(), eye.y(), eye.z());
      
    GLint uniformLight = glGetUniformLocation( shader, "light" );
    Camera::Quaternion light(0., -1., 1., -10.);
    light = r.conjugate() * light * r;
    glUniform3f(uniformLight, light.x(), light.y(), light.z());

    camera.setGLModelView();
    callDisplayList();

    shader.disable();
    glutSwapBuffers();
}
   
void 
Viewer::updateDisplayList()
{
    if (surfaceDL)
    {
        glDeleteLists(surfaceDL, 1);
        surfaceDL = 0;
    }

    surfaceDL = glGenLists(1);
    glNewList(surfaceDL, GL_COMPILE);
    drawScene();
    glEndList();
}
   
void Viewer::setGL()
{
    glClearColor( 1., 1., 1., 1. );
    setLighting();
}
   
void Viewer::setLighting()
{
    GLfloat position[4] = { 20., 30., 40., 0. };
    glLightfv(GL_LIGHT0, GL_POSITION, position);
    glEnable(GL_NORMALIZE);
    glEnable(GL_LIGHT0);

    GLfloat  diffuse[4] = { .8, .5, .3, 1. };
    GLfloat specular[4] = { .3, .3, .3, 1. };
    GLfloat  ambient[4] = { .2, .2, .5, 1. };
      
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,   diffuse );
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,   ambient );
    glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, 16.     );
}
   
void 
Viewer::callDisplayList()
{
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glCallList(surfaceDL);
    glPopAttrib();
}
   
void 
Viewer::drawScene()
{
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1., 1.);
    drawPolygons();
    glDisable(GL_POLYGON_OFFSET_FILL);
    if (renderWireframe) drawWireframe();
    if (renderSelected ) drawSelectedVerts();
    drawSelectedVerts();
    glPopAttrib();
}

void 
Viewer::drawPolygons()
{
    double minVal = values.minCoeff();
    double maxVal = values.maxCoeff();
    double range  = maxVal - minVal;
    if (range == 0.0) range = 1.0;

    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);

    std::vector<Vec3> vertNormal;
    meshPtr->computeVertNormals(vertNormal);

    glBegin(GL_TRIANGLES);
    for (unsigned face = 0; face < meshPtr->numFaces(); ++face) 
    {
        if (renderWireframe) 
        {
            Vec3 n = meshPtr->computeFaceNormal(face);
            glNormal3dv(n.data());
        }

        std::vector<unsigned> fVerts;
        meshPtr->getFaceVerts(face, fVerts);
        assert(fVerts.size() == 3);

        for (unsigned i = 0; i < 3; ++i) 
        {
            unsigned vert = fVerts[i];
            Vec3 p = meshPtr->getVertPos(vert);
            if (!renderWireframe) glNormal3dv(vertNormal[vert].data());

            double alpha = (values(vert) - minVal) / range;         
            glColor3f(alpha, alpha, alpha);
            glVertex3dv(p.data());
        }
    }
    glEnd();

    glDisable(GL_COLOR_MATERIAL);
}
   
void 
Viewer::drawWireframe()
{
    shader.disable();   
    glDisable(GL_LIGHTING);
    glColor4f(0., 0., 0., 1.);

    glBegin(GL_LINES);
    for (unsigned face = 0; face < meshPtr->numFaces(); ++face) 
    {
        std::vector<unsigned> fVerts;
        meshPtr->getFaceVerts(face, fVerts);
        assert(fVerts.size() == 3);

        for (unsigned i = 0; i < 3; ++i) 
        {
            Vec3 p0 = meshPtr->getVertPos(fVerts[(i+1)%3]);
            Vec3 p1 = meshPtr->getVertPos(fVerts[(i+2)%3]);
            glVertex3dv(p0.data());
            glVertex3dv(p1.data());
        }
    }
    glEnd();
}
   
void 
Viewer::drawSelectedVerts( )
{
    shader.disable();
    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    glPushAttrib(GL_ALL_ATTRIB_BITS);

    glEnable(GL_COLOR_MATERIAL);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glColor3f(0.0, 0.5, 0.5);
    double h = 0.5 * meshPtr->computeMaxEdgeLength();
    for (std::set<unsigned>::const_iterator 
        it = selectedVert.begin(); it != selectedVert.end(); ++it) 
    {
        Vec3 p = meshPtr->getVertPos(*it);
        glPushMatrix();
        glTranslated(p.x(), p.y(), p.z());
        glutSolidSphere(h, 10, 10);
        glPopMatrix();
    }

    glPopAttrib();
}

void 
Viewer::drawVerts()
{
    for (unsigned i = 0; i < meshPtr->numVerts(); ++i)
    {
        Vec3 p = meshPtr->getVertPos(i);
        glLoadName(i);
        glBegin(GL_POINTS);
        glVertex3dv(p.data());
        glEnd();
    }
}

void 
Viewer::pickVertex(int x, int y)
{
    int width  = glutGet(GLUT_WINDOW_WIDTH );
    int height = glutGet(GLUT_WINDOW_HEIGHT);
    if (x < 0 || x >= width || y < 0 || y >= height) return;

    int bufSize = meshPtr->numVerts();
    GLuint* buf = new GLuint[bufSize];
    glSelectBuffer(bufSize, buf);

    GLint viewport[4];
    GLdouble projection[16];
    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetDoublev(GL_PROJECTION_MATRIX, projection);

    glRenderMode(GL_SELECT);
    glInitNames();
    glPushName(0);

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluPickMatrix(x, viewport[3]-y, 10, 10, viewport);
    glMultMatrixd(projection);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    drawVerts();
    glPopMatrix();

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();

    glMatrixMode(GL_MODELVIEW);
    long hits = glRenderMode(GL_RENDER);

    int index = -1;
    double min_z = 1.0e100;
    for( long i = 0; i < hits; ++i ) 
    {
        double distance = buf[4*i + 1];
        if (distance < min_z) 
        {
            index = buf[4*i + 3];  
            min_z = distance;
        }
    }
    delete[] buf;
    if (index < 0) return;

    std::set<unsigned>::iterator it = selectedVert.find(index);
    if (it == selectedVert.end())
        selectedVert.insert(index);
    else
        selectedVert.erase(it);
    updateDisplayList();
}

void
Viewer::clearData()
{
    values = VectorXd::Zero(meshPtr->numVerts());
    selectedVert.clear();
}

void
Viewer::pickCenter()
{
    unsigned index = meshPtr->numVerts();
    double minVal = std::numeric_limits<int>::max();
    for (unsigned vert = 0; vert < meshPtr->numVerts(); ++vert) 
    {
        Vec3 p = meshPtr->getVertPos(vert);
        double val = p.norm();
        if (val < minVal) 
        {
            minVal = val;
            index = vert;
        }
    }
    selectedVert.insert(index);    
}

void 
Viewer::runGeodesics()
{
    // pick at least one source vert
    if (selectedVert.empty()) pickCenter();

    // precompute vert area
    VectorXd vertArea;
    meshPtr->computeVertArea(vertArea);

    // set source
    VectorXd src = VectorXd::Zero(meshPtr->numVerts());
    for(std::set<unsigned>::const_iterator
        it = selectedVert.begin(); it != selectedVert.end(); ++it) 
    {
        unsigned vert = *it;
        src(vert) = 1.0 / vertArea(vert);
    }
    src = solverPtr->applyKernel(src);

    // init values
    values = VectorXd::Zero(meshPtr->numVerts());

    #pragma omp parallel for
    for (int vert = 0; vert < (int) meshPtr->numVerts(); ++vert) 
    {
        // per vert
        VectorXd dst = VectorXd::Zero(meshPtr->numVerts());
        dst(vert) = 1.0 / vertArea(vert);
        dst = solverPtr->applyKernel(dst);

        // call OT
        values(vert) = solverPtr->computeDistance(src, dst, verbose);
    }
    
    // debug
    std::cout << "MinDist = " << values.minCoeff() << std::endl; 
    std::cout << "MaxDist = " << values.maxCoeff() << std::endl; 
    //
}

void 
Viewer::runInterpolation()
{
    if (selectedVert.size() != 2)
    {
        std::cout << "Interpolation needs just two verts" << std::endl;
        return;
    }

    std::set<unsigned>::const_iterator it = selectedVert.begin();
    unsigned vert0 = *it++;
    unsigned vert1 = *it;
    std::cout << "Verts: " << vert0 << " ; " << vert1 << std::endl;

    VectorXd vertArea;
    meshPtr->computeVertArea(vertArea);

    std::vector<VectorXd> p(2);
    p[0] = VectorXd::Zero(meshPtr->numVerts());
    p[0](vert0) = 1.0 / vertArea(vert0);
    p[0] = solverPtr->applyKernel(p[0]);
    p[1] = VectorXd::Zero(meshPtr->numVerts());
    p[1](vert1) = 1.0 / vertArea(vert1);
    p[1] = solverPtr->applyKernel(p[1]);

    double entropy0 = solverPtr->computeMarginalEntropy(p[0]);
    double entropy1 = solverPtr->computeMarginalEntropy(p[1]);
    ot::Options& opt = solverPtr->options();
    opt.upperEntropy = std::max(entropy0, entropy1);
    std::cout << "MaxEntropy: " << opt.upperEntropy << std::endl;
    
    VectorXd alpha(2);
    alpha(0) = 1. - timestep;
    alpha(1) = timestep;

    solverPtr->computeBarycenter(p, alpha, values, false, verbose);
    std::cout << "interpolation " << timestep << " done" << std::endl;
    
    double entropy = solverPtr->computeMarginalEntropy(values);
    std::cout << "BarycenterEntropy = " << entropy << std::endl;

    // debug
    std::cout << "MinDist = " << values.minCoeff() << std::endl; 
    std::cout << "MaxDist = " << values.maxCoeff() << std::endl; 
    //
}

void Viewer::takeScreenshot()
{
    static int index = 0;
  
    GLint view[4];
    glGetIntegerv(GL_VIEWPORT, view);
    int w = view[2];
    int h = view[3];

    // get pixels
    Image image(w, h);
    glReadPixels(0, 0, w, h, GL_BGR, GL_FLOAT, &image(0,0));

    std::stringstream filename;
    filename << "snapshot" << index << ".tga";
    image.write(filename.str().c_str());
    std::cout << "snapshot " << index << " taken" << std::endl;

    index++;
}
