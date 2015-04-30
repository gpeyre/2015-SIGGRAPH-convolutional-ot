#include "Camera.h"

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glut.h>
#endif

Camera::Quaternion 
Camera::projectClickToSphere(int x, int y)
{
   GLint viewport[4];
   glGetIntegerv(GL_VIEWPORT, viewport);   
   int w = viewport[2];
   int h = viewport[3];

   Quaternion p(0.,
                2. * double(x) / double(w) - 1.,
                2. * double(y) / double(h) - 1.,
                0.);
   if(p.squaredNorm() > 1.) 
   {
      p.normalize();
      p.vec().z() = 0.;
   }
   else
   {
      p.vec().z() = std::sqrt(1. - p.squaredNorm());
   }
   return p;
}

Camera::Quaternion 
Camera::computeCurrentRotation() const
{
   return (mDrag * mClick.conjugate()) * mRot;
}

void 
Camera::setGLModelView() const
{
   Quaternion r = computeCurrentRotation();
   double w = r.w();
   double x = r.x();
   double y = r.y();
   double z = r.z();

   GLdouble M[16] = 
   {
      1.-2.*y*y-2.*z*z, 2.*x*y+2.*w*z, 2.*x*z-2.*w*y, 0.,
      2.*x*y-2.*w*z, 1.-2.*x*x-2.*z*z, 2.*y*z+2.*w*x, 0.,
      2.*x*z+2.*w*y, 2.*y*z-2.*w*x, 1.-2.*x*x-2.*y*y, 0.,
      0., 0., 0., 1.
   };

   glMatrixMode(GL_MODELVIEW);
   glMultMatrixd(M);
}

void 
Camera::mouse(int button, int state, int x, int y)
{
   // left button
   if (button == 0) {
      if (state == GLUT_DOWN)
      {
         mClick = mDrag = mLast = projectClickToSphere(x, y);
      }
      if (state == GLUT_UP)
      {

         mRot = computeCurrentRotation();
         mClick = mDrag = Quaternion(1.,0.,0.,0.);
      }
      return;
   }

   // mouse wheel
   if (button == 3) {
      zoomIn();
      return;
   } 
   if (button == 4) {
      zoomOut();
      return;
   }
}

void 
Camera::motion(int x, int y)
{
   mLast = mDrag;
   mDrag = projectClickToSphere(x, y);
}

void 
Camera::zoomIn()
{
   mZoom -= 0.1;
}

void 
Camera::zoomOut()
{
   mZoom += 0.1;
}
