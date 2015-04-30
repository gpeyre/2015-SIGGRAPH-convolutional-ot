#include <iostream>
#include <math.h>
#include <map>
#include <cstring>
#include <vector>
#include "blVector3D.hpp"
#include "mex.h"
using namespace std;

extern void _main();

// vertices, faces, vectors, bases, radius, nRadialSamples
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray *prhs[]) {
	if (nrhs != 6)
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin", "ADDCYLINDERSHELPER requires six inputs.");
    else if (nlhs != 3) // vertices, faces, indicator
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout", "ADDCYLINDERSHELPER produces three outputs.");

    // Read list of vertices
    double *verticesMtx = mxGetPr(prhs[0]);
    int nv = (int)mxGetM(prhs[0]);
    int three = (int)mxGetN(prhs[0]);
    
    if (three != 3) mexErrMsgTxt("Vertices must have three columns.");
    
    double *vertexCoords[3];
    vertexCoords[0] = verticesMtx;
    vertexCoords[1] = verticesMtx+nv;
    vertexCoords[2] = verticesMtx+nv*2;
    
    vector< blVector3d<double> > vertices;
    blVector3d<double> v;
    for (int i = 0; i < nv; i++) {
        for (int j = 0; j < 3; j++) v[j] = vertexCoords[j][i];
        vertices.push_back(v);
    }
    
    // Read list of faces
    double *facesMtx = mxGetPr(prhs[1]);
    int nf = (int)mxGetM(prhs[1]);
    three = (int)mxGetN(prhs[1]);
    
    if (three != 3) mexErrMsgTxt("Faces must have three columns.");
    
    double *faceCoords[3];
    faceCoords[0] = facesMtx;
    faceCoords[1] = facesMtx+nf;
    faceCoords[2] = facesMtx+nf*2;
    
    vector< blVector3d<int> > faces;
    blVector3d<int> f;
    for (int i = 0; i < nf; i++) {
        for (int j = 0; j < 3; j++) f[j] = (int)faceCoords[j][i];
        faces.push_back(f);
    }
    
    // Read vector field
    double *vfMtx = mxGetPr(prhs[2]);
    int nvf = (int)mxGetM(prhs[2]);
    three = (int)mxGetN(prhs[2]);
    
    if (three != 3) mexErrMsgTxt("Vector field must have three columns.");
    
    double *vfCoords[3];
    vfCoords[0] = vfMtx;
    vfCoords[1] = vfMtx+nvf;
    vfCoords[2] = vfMtx+nvf*2;
    
    vector< blVector3d<double> > vectorField;
    for (int i = 0; i < nvf; i++) {
        for (int j = 0; j < 3; j++) v[j] = vfCoords[j][i];
        vectorField.push_back(v);
    }
    
    // Read vector field base points
    double *baseMtx = mxGetPr(prhs[3]);
    int nvf2 = (int)mxGetM(prhs[3]);
    three = (int)mxGetN(prhs[3]);
    
    if (three != 3) mexErrMsgTxt("Bases must have three columns.");
    if (nvf != nvf2) mexErrMsgTxt("Bases must have same size as vector field.");
    
    double *vfBaseCoords[3];
    vfBaseCoords[0] = baseMtx;
    vfBaseCoords[1] = baseMtx+nvf;
    vfBaseCoords[2] = baseMtx+nvf*2;
    
    vector< blVector3d<double> > vectorFieldBase;
    for (int i = 0; i < nvf; i++) {
        for (int j = 0; j < 3; j++) v[j] = vfBaseCoords[j][i];
        vectorFieldBase.push_back(v);
    }
    
    double radius = *mxGetPr(prhs[4]);
    int nRadialSamples = (int)*mxGetPr(prhs[5]);
    
    // DO WORK HERE ////////////////////////////////////////////////////////////////////
    
    for (int i = 0; i < nvf; i++) {
        blVector3d<double> vec = vectorField[i];
        blVector3d<double> base = vectorFieldBase[i];
        
        if (vec.Norm2() < 1e-9) continue; // skip zero vectors
        
        blVector3d<double> v1, v2;
        if (fabs(vec.y) > fabs(vec.x)) {
            v1.x = 0;
            v1.y = vec.z;
            v1.z = -vec.y;
        } else {
            v1.x = -vec.z;
            v1.y = 0;
            v1.z = vec.x;
        }
        v1.NormalizeVector();
        
        v2 = CrossProduct(vec, v1);
        v2.NormalizeVector();
        
        int startVtx = (int)vertices.size()+1;
        for (int j = 0; j < nRadialSamples; j++) {
            double theta = (double)j/(double)nRadialSamples * 2*3.14159265359;
            blVector3d<double> displacement = (v1*cos(theta)+v2*sin(theta))*radius;
            vertices.push_back(base+displacement);
            vertices.push_back(base+displacement+vec);
            
            blVector3d<int> newFace;
            newFace[0] = startVtx + (j*2)%(2*nRadialSamples);
            newFace[1] = startVtx + (j*2 + 1)%(2*nRadialSamples);
            newFace[2] = startVtx + (j*2 + 2)%(2*nRadialSamples);
            faces.push_back(newFace);
            
            newFace[2] = startVtx + (j*2 + 3)%(2*nRadialSamples);
            newFace[1] = startVtx + (j*2 + 1)%(2*nRadialSamples);
            newFace[0] = startVtx + (j*2 + 2)%(2*nRadialSamples);
            faces.push_back(newFace);
            
            // add cap
            if (j>0 && j < nRadialSamples - 1) {
                newFace[0] = startVtx;
                newFace[1] = startVtx + (j*2)%(2*nRadialSamples);
                newFace[2] = startVtx + (j*2 + 2)%(2*nRadialSamples);
                faces.push_back(newFace);
                
                newFace[0] = 1+startVtx;
                newFace[1] = startVtx + (j*2 + 1)%(2*nRadialSamples);
                newFace[2] = startVtx + (j*2 + 3)%(2*nRadialSamples);
                faces.push_back(newFace);
            }
        }
    }
    
    ////////////////////////////////////////////////////////////////////////////////////
    
    
    double *verticesResult = (double*)mxCalloc(3*vertices.size(),sizeof(double));
    int pos = 0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < vertices.size(); j++)
            verticesResult[pos++] = vertices[j][i];
    plhs[0] = mxCreateNumericMatrix(0,0,mxDOUBLE_CLASS,mxREAL);
    mxSetData(plhs[0],verticesResult);
    mxSetM(plhs[0],vertices.size());
    mxSetN(plhs[0],3);
    
    double *facesResult = (double*)mxCalloc(3*faces.size(),sizeof(double));
    pos = 0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < faces.size(); j++)
            facesResult[pos++] = (double)faces[j][i];
    plhs[1] = mxCreateNumericMatrix(0,0,mxDOUBLE_CLASS,mxREAL);
    mxSetData(plhs[1],facesResult);
    mxSetM(plhs[1],faces.size());
    mxSetN(plhs[1],3);
    
    double *indicatorResult = (double*)mxCalloc(vertices.size(),sizeof(double));
    pos = 0;
    for (int i = 0; i < nv; i++)
        indicatorResult[pos++] = 0;
    for (int i = nv; i < vertices.size(); i++)
        indicatorResult[pos++] = 1;
    
    plhs[2] = mxCreateNumericMatrix(0,0,mxDOUBLE_CLASS,mxREAL);
    mxSetData(plhs[2],indicatorResult);
    mxSetM(plhs[2],vertices.size());
    mxSetN(plhs[2],1);
    
}