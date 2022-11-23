#pragma once

#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "LinAlgTools.h"
#include "tools/InputParser.h"
#include "tools/CommonOperations.h"
#include "tools/Algorithm.h"
#include "happly.h"

#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <memory>
#include <numeric>
#include <unordered_map>
#include <iomanip>

struct vertex{
    using Real = CommonTypes::Real;
    using Real3= CommonTypes::Real3;

    Real3 position_ = {{0,0,0}};
    Real3 normals_ = {{0,0,0}};

    Real& operator[](int i){return position_[i];}
};

struct triangle{
    using  INT3 = CommonTypes::index3;

    INT3 triangleindices_;

    // each triangle has 3 vertices
    std::array<vertex,3> vertices_;

    int& operator[](int i){return triangleindices_[i];}
    int operator[](int i)const {return triangleindices_[i];}
};

class Mesh{
    public:
        using Real3 = CommonTypes::Real3;
        using Real  = CommonTypes::Real;
        using INT2  = CommonTypes::index2;
        using INT3  = CommonTypes::index3;

        // default constructor and destructor
        Mesh() {};
        ~Mesh(){};

                                        /*******************************************
                                         * ******** setting Function **************
                                         * ****************************************/
        void setBoxLength(Real3 box) { boxLength_ = box;setPBC(true);}
        void setPBC(bool isPbc) {isPeriodic_=isPbc;}

                                        /********************************************
                                         * ******** accessing Function **************
                                         * *****************************************/
        std::vector<vertex>& accessvertices() {return vertices_;}
        std::vector<triangle>& accesstriangles() {return triangles_;}
        std::vector<Real>& accessTriangleArea() {return triangleArea_;}
        std::vector<Real3>& accessPerVertexDir1() {return PerVertexdir1_;}
        std::vector<Real3>& accessPerVertexDir2() {return PerVertexdir2_;}

                                        /********************************************
                                         * **********  getter Function **************
                                         * *****************************************/
        const std::vector<vertex>& getvertices() const {return vertices_;}
        std::vector<Real3> getVertexPositions();
        const std::vector<triangle>& gettriangles() const {return triangles_;}
        const std::vector<Real>& getTriangleArea() const {return triangleArea_;}
        const std::vector<Real3>& getPerVertexDir1() const {return PerVertexdir1_;}
        const std::vector<Real3>& getPerVertexDir2() const {return PerVertexdir2_;}
        const std::vector<Real3>& getFaceNormals() const {return facetNormals_;}
        const std::vector<std::vector<int>>& getMapVertexToFace() const {return MapVertexIndicesToFaceIndices_;}
        const std::vector<Real3>& getCornerAreas() const {return cornerArea_;}
        int getNumVertices() const {return vertices_.size();}
        int getNumTriangles() const {return triangles_.size();}
        Real3 getBoxLength() const {return boxLength_;}
        const std::map<INT2, std::vector<int>> getMapEdgeToFace() const {return MapEdgeToFace_;}

        // Scale all the vertices by a single number 
        void scaleVertices(Real num);

        // Calculate the corner area according to trimesh code C++
        void CalculateCornerArea();

        // calculate volume
        Real calculateVolume();

        // function that calculates the normals of each of the vertex --> weighted by area 
        void CalcVertexNormalsAreaWeighted();

        // function that calculates the normals of each of the vertex --> not weighted by anything and updates the normals in each vertices
        void CalcVertexNormals();

        // Find the vertex direction 
        void CalcPerVertexDir();

        // update the triangles/vertex/edges if they were to be changed
        void update();
                                        ///////////////////////////////
                                        /////       PBC stuff   ///////
                                        ///////////////////////////////

        // Check whether or not the mesh is periodic 
        bool isPeriodic() {return isPeriodic_;}
        // find PBC distance between 2 vertices 
        void getVertexDistance(const vertex& v1, const vertex& v2, Real3& distVec, Real& dist);
        void getVertexDistance(const Real3& v1, const Real3& v2, Real3& distVec, Real& dist);
        void CalculateShift(const Real3& v1, const Real3& v2, Real3& shiftVec);
        // move a vertex into pbc box
        void MoveVertexIntoBox(const Real3& OldVertPos, Real3& NewVertexPos);
        // find shifted vertex position
        Real3 getShiftedVertexPosition(const vertex& v1, const vertex& v2);
        // find the shit that takes a position into the box 
        Real3 getShiftIntoBox(const Real3& v1);
    
    private:
        // vertices and triangles in the mesh 
        std::vector<vertex> vertices_;
        std::vector<triangle> triangles_;

        // the Areas and normals
        std::vector<Real> triangleArea_;
        std::vector<Real3> cornerArea_;
        std::vector<Real3> facetNormals_;
        std::vector<Real3> vertexNormals_;
        std::vector<Real3> PerVertexdir1_;
        std::vector<Real3> PerVertexdir2_;

        // map the edge to face 
        std::map<INT2, std::vector<int>> MapEdgeToFace_;

        // a map from vertex indices to the face indices  
        std::vector<std::vector<int>> MapVertexIndicesToFaceIndices_;

        // map vertex to edges 
        std::vector<std::vector<INT2>> MapVertexToEdges_;

        // the factor to which all the vertices on the mesh will be multiplied by
        Real factor_=1.0;

        // whether or not the mesh is periodic
        bool isPeriodic_=false;
        Real3 boxLength_;
};


namespace MeshTools
{
    using Real3 = CommonTypes::Real3;
    using Real  = CommonTypes::Real;
    using INT3  = CommonTypes::index3;
    using INT2  = CommonTypes::index2;

    bool readPLY(std::string& filename, Mesh& mesh);
    bool readPLYlibr(std::string& filename, Mesh& mesh);

    // write PLY file
    void writePLY(std::string filename, const std::vector<Real3>& Vertices, const std::vector<INT3>& faces, Real factor=1.0);
    void writePLY(std::string filename, const std::vector<Real3>& Vertices, const std::vector<INT3>& faces, const std::vector<Real3>& normals);
    void writePLYRGB(std::string filename, const std::vector<Real3>& Vertices, const std::vector<INT3>& faces, const std::vector<Real3>& RGB);
    void writePLY(std::string filename, Mesh& mesh);

    // write non pbc mesh 
    void writeNonPBCMesh(std::string filename, Mesh& mesh);
    void writeNonPeriodicTriangleIndices(std::string name, Mesh& mesh);

    // auxiliary functions for writing mesh areas
    void writeMeshArea(std::string filename, Mesh& mesh);
    void writeCuttedMesh(std::string filename, Mesh& mesh, Real3& volume);

    // write STL file
    void writeSTL(std::string filename, Mesh& mesh);
    void writeSTL(std::string filename, const std::vector<Real3>& vertices, const std::vector<INT3>& faces);

    // calculate PBC distance
    Real3 calculateShift(const Real3& vec1, const Real3& vec2, const Real3& boxLength);

    // calculate pbc distance 
    void calculateDistance(const Real3& vec1, const Real3& vec2, const Real3& boxlength, Real3& distance, Real& distsq);

    // find if an edge is periodic 
    bool isPeriodicEdge(const Real3& vec1, const Real3& vec2, Real3& newarr, const Real3& boxLength);

    // map vertices to faces
    void MapVerticesToFaces(Mesh& mesh, std::vector<std::vector<int>>& map);

    // calculate triangle areas and Face normals 
    void CalculateTriangleAreasAndFaceNormals(Mesh& mesh, std::vector<Real>& Areas, std::vector<Real3>& Normals);

    // find vertex neighbors 
    void CalculateVertexNeighbors(Mesh& mesh, std::vector<std::vector<int>>& neighborIndices);

    // map Edge to faces
    void MapEdgeToFace(Mesh& mesh, std::map<INT2,std::vector<int>>& magEdgeToFace, std::vector<std::vector<INT2>>& MapVertexToEdge);

    // map edge to opposing vertices --> needs boundary indicators  
    void MapEdgeToOpposingVertices(Mesh& mesh, std::map<INT2, std::vector<int>>& mapEdgeToFace,std::map<INT2, std::vector<int>>& MapEdgeToOppoVertices);

    // find boundary vertices 
    void CalculateBoundaryVertices(Mesh& mesh, const std::map<INT2, std::vector<int>>& mapEdgeToFace, std::vector<bool>& boundaryIndicator);

    // check if point is on boundary
    bool IsBoundary(int Index, const std::vector<bool>& boundaryIndicator);

    // convert a PBC mesh to non PBC mesh --> usually for visualization purposes
    void ConvertToNonPBCMesh(Mesh& mesh, std::vector<Real3>& vertices, std::vector<INT3>& faces);

    // check if a particular triangle is periodic
    bool IsPeriodicTriangle(std::vector<vertex>& Vertices,INT3& face, Real3 BoxLength);
    bool IsPeriodicTriangle(const Mesh& mesh, int faceindex);

    // shift a triangle 
    void ShiftPeriodicTriangle(const std::vector<vertex>& Vertices, const INT3& faces, Real3 BoxLength, Real3& A, Real3& B, Real3& C);

    // make an Edge, edge is simply the 2 indices of {{minIndex, maxIndex}}
    INT2 makeEdge(int i, int j);

    // calculate the corner area 
    void CalculateCornerArea(Mesh& mesh, std::vector<Real3>& CornerArea, std::vector<Real>& VertexArea);

    // cut mesh and give a new mesh 
    void CutMesh(Mesh& mesh, std::vector<INT3>& face, std::vector<Real3>& vertices, Real3 volume);

    // cut mesh and give a new mesh
    void CutMesh(Mesh& mesh, Real3 volume);

    // Moller Trumbore Ray-Triangle intersection method 
    bool MTRayTriangleIntersection(Real3& A, Real3& B, Real3& C, Real3& O, Real3& D, Real& t, Real& u, Real& v);

    // check degenerate triangles 
    void CheckDegenerateTriangle(Mesh& mesh, \
                                             std::vector<int>& MergeFaces,\
                                             std::vector<INT2>& MergeVertices);

    // decimate degenerate triangles
    bool decimateDegenerateTriangle(Mesh& mesh);

    // calculate barycenter of all the boundary vertices 
    void CalculateBoundaryBarycenter(Mesh& mesh, std::vector<bool>& boundaryIndicator, Real3& barycenter);

    // correct  Mesh --> Many times, we need to cut certain triangles out for certain reasons
    // this function takes care of rearranging the mesh vertices and triangles
    void CorrectMesh(Mesh& mesh, std::vector<int>& FaceIndices);

    // check if a triangle is isolated --> all of its edges are boundaries
    bool IsIsolatedFace(Mesh& mesh, int faceIndex, const std::map<INT2, std::vector<int>>& mapEdgeToFace);

    // check if a triangle is teeth like --> 2 of its edges are boundaries
    bool IsTeethlikeFace(Mesh& mesh, int faceIndex, const std::map<INT2, std::vector<int>>& mapEdgeToFace);

    // regenerate mesh after some faces are cut
    void ReconstructMeshAfterFaceCut(Mesh& mesh);
};