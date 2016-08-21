#include "TriMesh.h"
#include "Triangle.h"

#include <map>
#include <fstream>
#include <sstream>

void TriMesh::getFaceVerts(unsigned face, std::vector<unsigned>& verts) const
{
    verts = mFaceToVert[face];
}

unsigned TriMesh::getVertIndexInFace(unsigned face, unsigned vert) const
{
    std::vector<unsigned> fVert;
    getFaceVerts(face, fVert);
    assert(fVert.size() == 3);

    for (unsigned i = 0; i < 3; ++i) 
        if (fVert[i] == vert) return i;

    assert(false);
    return 0;    
}

unsigned TriMesh::getFaceVert(unsigned face, unsigned index) const
{
    return mFaceToVert[face][index];
}

double TriMesh::computeFaceAngle(unsigned face, unsigned vertInFace) const 
{
    std::vector<unsigned> fVert;
    getFaceVerts(face, fVert);
    assert(fVert.size() == 3);

    Vec3 p[3];
    for (unsigned i = 0; i < 3; ++i)
        p[i] = getVertPos(fVert[(vertInFace+i)%3]);
    return Triangle::angle(p[0], p[1], p[2]);
}

double TriMesh::computeFaceCotan(unsigned face, unsigned vertInFace) const
{
    std::vector<unsigned> fVert;
    getFaceVerts(face, fVert);
    assert(fVert.size() == 3);

    Vec3 p[3];
    for (unsigned i = 0; i < 3; ++i)
        p[i] = getVertPos(fVert[(vertInFace+i)%3]);
    return Triangle::cotan(p[0], p[1], p[2]);
}

double TriMesh::computeFaceArea(unsigned face) const
{
    std::vector<unsigned> fVert;
    getFaceVerts(face, fVert);
    assert(fVert.size() == 3);

    Vec3 p0 = getVertPos(fVert[0]);
    Vec3 p1 = getVertPos(fVert[1]);
    Vec3 p2 = getVertPos(fVert[2]);
    return Triangle::area(p0, p1, p2);
}

Vec3 TriMesh::computeFaceNormal(unsigned face) const
{
    std::vector<unsigned> fVert;
    getFaceVerts(face, fVert);
    assert(fVert.size() == 3);

    Vec3 p0 = getVertPos(fVert[0]);
    Vec3 p1 = getVertPos(fVert[1]);
    Vec3 p2 = getVertPos(fVert[2]);
    return (p1-p0).cross(p2-p1).normalized();
}

double TriMesh::computeTotalArea() const
{
    double sum = 0.0;
    for (unsigned face = 0; face < numFaces(); ++face) 
        sum += computeFaceArea(face);
    return sum;
}

double TriMesh::computeMaxEdgeLength() const
{
    double maxLen = 0.0;
    for (unsigned face = 0; face < numFaces(); ++face) 
    {
        std::vector<unsigned> fVerts;
        getFaceVerts(face, fVerts);
        assert(fVerts.size() == 3);

        for (unsigned i = 0; i < 3; ++i) 
        {
            unsigned vert0 = fVerts[(i+1)%3];
            unsigned vert1 = fVerts[(i+2)%3];
            double len = computeEdgeLength(vert0, vert1);
            maxLen = std::max(maxLen, len);
        }
    }
    return maxLen;
}

void TriMesh::computeVertArea(VectorXd& vertArea) const
{
    vertArea = VectorXd::Zero(numVerts());
    for (unsigned face = 0; face < numFaces(); ++face) 
    {
        std::vector<unsigned> fVert;
        getFaceVerts(face, fVert);
        assert(fVert.size() == 3);

        double faceArea = computeFaceArea(face);
        for (unsigned i = 0; i < 3; ++i) vertArea(fVert[i]) += faceArea;
    }
    vertArea /= 3.0;
}

void TriMesh::computeVertNormals(std::vector<Vec3>& vertNormal) const
{
    vertNormal = std::vector<Vec3>(numVerts(), Vec3::Zero());
    
    for (unsigned face = 0; face < numFaces(); ++face) 
    {
        std::vector<unsigned> fVert;
        getFaceVerts(face, fVert);
        assert(fVert.size() == 3);

        double faceArea = computeFaceArea(face);
        Vec3 faceNormal = computeFaceNormal(face);
        for (unsigned i = 0; i < 3; ++i) 
            vertNormal[fVert[i]] += faceArea * faceNormal;
    }

    for (unsigned vert = 0; vert < numVerts(); ++vert) 
    {
        vertNormal[vert].normalize();
    }
}

void TriMesh::buildHeatKernel(SparseMatrix& A, double timestep) const
{
    VectorXd diag;
    computeVertArea(diag);
    std::map< std::pair<unsigned,unsigned>, double > offdiag; 
    for (unsigned face = 0; face < numFaces(); ++face) 
    {
        std::vector<unsigned> fVert;
        getFaceVerts(face, fVert);
        assert(fVert.size() == 3);

        for (unsigned i = 0; i < 3; ++i)
        {
            unsigned vert0 = fVert[(i+1)%3];
            unsigned vert1 = fVert[(i+2)%3];
            if (vert0 > vert1) std::swap(vert0, vert1);
            std::pair<unsigned,unsigned> edge(vert0,vert1);

            double value = timestep * 0.5 * computeFaceCotan(face, i);
            diag(vert0) += value;
            diag(vert1) += value;
            if (offdiag.find(edge) == offdiag.end())
                offdiag[edge] = -value;
            else
                offdiag[edge] -= value;
        }
    }

    std::vector<Triplet> triplet;
    for(unsigned i = 0; i < diag.size(); ++i)
    {
        triplet.push_back(Triplet(i,i,diag(i)));
    }
    for(std::map< std::pair<unsigned,unsigned>, double >::const_iterator 
        it = offdiag.begin(); it != offdiag.end(); ++it)
    {
        triplet.push_back(Triplet(it->first.first,it->first.second,it->second));
        triplet.push_back(Triplet(it->first.second,it->first.first,it->second));        
    } 
    
    A.resize(numVerts(), numVerts());
    A.setFromTriplets(triplet.begin(),triplet.end());
}

double TriMesh::computeEdgeLength(unsigned vert0, unsigned vert1) const
{
    Vec3 p0 = getVertPos(vert0);
    Vec3 p1 = getVertPos(vert1);
    return (p1-p0).norm();
}

void TriMesh::computeFaceGrad(const VectorXd& vertForm, std::vector<Vec3>& faceVector) const
{
    faceVector.resize(numFaces());
    for (unsigned face = 0; face < numFaces(); ++face) 
    {
        std::vector<unsigned> fVert;
        getFaceVerts(face, fVert);
        assert(fVert.size() == 3);

        Vec3 pos[3];
        double val[3];
        for (unsigned i = 0; i < 3; ++i) 
        {
            val[i] = vertForm(fVert[i]);
            pos[i] = getVertPos(fVert[i]);
        }

        Vec3 X = Vec3::Zero();
        for (unsigned i = 0; i < 3; ++i) 
        {
            unsigned j = (i + 1) % 3;
            unsigned k = (i + 2) % 3;
            X += val[i] * (pos[k] - pos[j]);
        }

        Vec3 n = (pos[1]-pos[0]).cross(pos[2]-pos[1]);
        double twice_area = n.norm();
        n /= twice_area;

        faceVector[face] = n.cross(X) / twice_area;
    }
}

void TriMesh::computeDivergence(const std::vector<Vec3>& faceVector, VectorXd& vertForm) const
{
    vertForm = VectorXd::Zero(numVerts());
    for (unsigned face = 0; face < numFaces(); ++face) 
    {
        std::vector<unsigned> fVert;
        getFaceVerts(face, fVert);
        assert(fVert.size() == 3);
        
        Vec3 n = computeFaceNormal(face);
        for (unsigned i = 0; i < 3; ++i) 
        {
            Vec3 p1 = getVertPos(fVert[(i+1)%3]);
            Vec3 p2 = getVertPos(fVert[(i+2)%3]);
            vertForm[fVert[i]] += 0.5*faceVector[face].dot(n.cross(p2-p1));
        }
    }
}

void TriMesh::clear()
{
    mPoints.clear();
    mFaceToVert.clear();
}

void TriMesh::read(const char* filename)
{
    clear();
    std::string line;
    std::ifstream in(filename);
    while(getline(in, line))
    {
        std::stringstream ss(line);
        std::string token;
        ss >> token;

        if (token == "v") 
        {
            double x, y, z;
            ss >> x >> y >> z;
            mPoints.push_back(Vec3(x,y,z));
            continue;
        }

        if (token == "f") 
        {
            std::vector<unsigned> face;
            while (ss >> token) 
            {
                unsigned index;
                std::string indexstring;
                std::stringstream tokenstream(token);
                getline(tokenstream, indexstring, '/');
                std::stringstream indexstream(indexstring);
                indexstream >> index;
                face.push_back(index-1);
            }
            mFaceToVert.push_back(face);
        }
    }
    in.close();
}

void TriMesh::normalize()
{
    Vec3 c = Vec3::Zero();
    for (unsigned vert = 0; vert < numVerts(); ++vert) {
        c += mPoints[vert];
    }
    c /= numVerts();
    for (unsigned vert = 0; vert < numVerts(); ++vert) {
        mPoints[vert] -= c;
    }

    double scale = std::sqrt(computeTotalArea());
    for (unsigned vert = 0; vert < numVerts(); ++vert) {
        mPoints[vert] /= scale;
    }
}
