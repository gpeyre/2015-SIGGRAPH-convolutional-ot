#ifndef VIEWER_H
#define VIEWER_H

#include "Camera.h"
#include "Shader.h"
#include "ot/types.h"
#include "ot/ConvSolver.h"
#include "ot/TriMesh.h"
#include <set>

class Viewer
{
public:
   // Parameters
   static int  windowSize[2];
   static bool renderWireframe;
   static bool renderSelected;

   // Data
   static std::set<unsigned> selectedVert;
   static TriMesh*        meshPtr;
   static ot::ConvSolver* solverPtr;
   static VectorXd        values;
   static double timestep;
   static int    verbose;
   
   // GL Data
   static Camera camera;
   static Shader shader;
   static GLuint surfaceDL;  

   // init
   static void init(int argc, char** argv);
   static void initGLUT(int argc, char** argv);
   static void initGLSL();

   // GLUT callbacks
   static void mouse(int button, int state, int x, int y);
   static void keyboard(unsigned char c, int x, int y);
   static void special(int i, int x, int y);
   static void motion(int x, int y);
   static void display();
   static void idle();

   // GL setup   
   static void setGL();
   static void setLighting();
   static void callDisplayList();
   static void updateDisplayList();

   // draw routines
   static void drawScene();
   static void drawVerts();
   static void drawPolygons();
   static void drawWireframe();
   static void drawSelectedVerts();
   static void pickVertex(int x, int y);
   static void takeScreenshot();

   // Application
   static void clearData();
   static void pickCenter();
   static void runGeodesics();
   static void runInterpolation();
};

#endif // VIEWER_H
