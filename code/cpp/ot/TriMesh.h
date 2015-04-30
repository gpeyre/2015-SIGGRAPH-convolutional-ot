#ifndef TRIMESH_H
#define TRIMESH_H

#include "types.h"
#include <vector>

class TriMesh
{
public:
    TriMesh() { }

    TriMesh(const TriMesh& mesh)
    : mPoints(mesh.mPoints)
    , mFaceToVert(mesh.mFaceToVert)
    { }

    TriMesh& operator = (const TriMesh& mesh)
    {
        mPoints = mesh.mPoints;
        mFaceToVert = mesh.mFaceToVert;
        return *this;
    }

    void clear();

    void normalize();

    void read(const char* filename);

    // count

    inline unsigned numVerts() const { return mPoints.size(); }

    inline unsigned numFaces() const { return mFaceToVert.size(); }
    
    inline const Vec3& getVertPos(unsigned vert) const { return mPoints[vert]; }

    // Face

    void getFaceVerts(unsigned face, std::vector<unsigned>& verts) const;

    unsigned getVertIndexInFace(unsigned face, unsigned vert) const;

    unsigned getFaceVert(unsigned face, unsigned index) const;

    double computeFaceAngle(unsigned face, unsigned vertInFace) const;

    double computeFaceCotan(unsigned face, unsigned vertInFace) const;

    double computeFaceArea(unsigned face) const;

    Vec3 computeFaceNormal(unsigned face) const;

    // Methods querying geometry

    double computeTotalArea() const;
    
    double computeMaxEdgeLength() const;

    void computeVertArea(VectorXd& vertArea) const;

    void computeVertNormals(std::vector<Vec3>& vertNormal) const;

    void buildHeatKernel(SparseMatrix& A, double timestep) const;

    double computeEdgeLength(unsigned vert0, unsigned vert1) const;

    void computeFaceGrad(const VectorXd& vertForm, std::vector<Vec3>& faceVector) const;

    void computeDivergence(const std::vector<Vec3>& faceVector, VectorXd& vertForm) const;

private:
    std::vector<Vec3> mPoints;
    std::vector< std::vector<unsigned> > mFaceToVert;
};

#endif // TRIMESH_H

