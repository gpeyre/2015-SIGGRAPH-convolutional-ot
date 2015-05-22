#ifndef CAMERA_H
#define CAMERA_H

#include <Eigen/Geometry>

class Camera
{
public:
   typedef Eigen::Quaternion<double> Quaternion;

   Camera() : 
   mClick( 1.,0.,0.,0.),
   mDrag ( 1.,0.,0.,0.),
   mLast ( 1.,0.,0.,0.),
   mRot  ( 0.,0.,1.,0.),
   mZoom (0.5)
   { }

   double zoom() const { return mZoom; }
   
   void setGLModelView() const;
   
   void mouse(int button, int state, int x, int y);
   
   void motion( int x, int y );
   
   void zoomIn();
   
   void zoomOut();
   
   Quaternion computeCurrentRotation() const;

private:
   Quaternion projectClickToSphere(int x, int y);

   Quaternion mClick;
   // mouse coordinates of current click
   
   Quaternion mDrag;
   // mouse coordinates of current drag
   
   Quaternion mLast;
   // mouse coordinates of previous drag
   
   Quaternion mRot;
   // camera rotation

   double mZoom;
   // zoom
};

#endif // CAMERA_H
