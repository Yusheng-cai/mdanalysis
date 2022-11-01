#include "Mesh.h"

void Mesh::MoveVertexIntoBox(const Real3& OldVerPos, Real3& NewVerPos)
{
    if (isPeriodic())
    {
        Real3 boxCenter;
        boxCenter = boxLength_ * 0.5;

        Real3 diff;
        for (int i=0;i<3;i++)
        {
            diff[i] = OldVerPos[i] - boxCenter[i];

            if (diff[i] >= boxLength_[i] * 0.5) diff[i] = diff[i] - boxLength_[i];
            if (diff[i] <= -boxLength_[i] * 0.5) diff[i] = diff[i] + boxLength_[i]; 

            NewVerPos[i] = boxCenter[i] + diff[i];
        }
    }
}

void Mesh::CalcPerVertexDir()
{
    PerVertexdir1_.resize(vertices_.size());
    PerVertexdir2_.resize(vertices_.size());

    for (int i=0;i<triangles_.size();i++)
    {
        auto& t = triangles_[i];
 
        int index0 = t.triangleindices_[0];
        int index1 = t.triangleindices_[1];
        int index2 = t.triangleindices_[2];

        Real3 diff0, diff1, diff2;
        Real diffsq0, diffsq1, diffsq2;

        getVertexDistance(vertices_[index1], vertices_[index0], diff0, diffsq0);
        getVertexDistance(vertices_[index2], vertices_[index1], diff1, diffsq1);
        getVertexDistance(vertices_[index0], vertices_[index2], diff2, diffsq2);

        PerVertexdir1_[index0] = diff0;
        PerVertexdir1_[index1] = diff1;
        PerVertexdir1_[index2] = diff2;
    }

    for (int i=0;i<PerVertexdir1_.size();i++)
    {
        PerVertexdir1_[i] = LinAlg3x3::CrossProduct(PerVertexdir1_[i], vertices_[i].normals_);
        LinAlg3x3::normalize(PerVertexdir1_[i]);

        Real3 B = LinAlg3x3::CrossProduct(vertices_[i].normals_, PerVertexdir1_[i]) ;
        LinAlg3x3::normalize(B);

        PerVertexdir2_[i] = B;
    } 
}

Mesh::Real Mesh::calculateVolume()
{
    Real volume_ =0.0;
    for (int i=0;i<triangles_.size();i++)
    {
        auto& t = triangles_[i];

        Real3 g;
        g.fill(0);
        for (int j=0;j<3;j++)
        {
            int index= t.triangleindices_[j];
            for (int k=0;k<3;k++)
            {
                g[k] += vertices_[index].position_[k] /3.0;
            }
        }

        Real3 vec1, vec2;
        int index1 = t.triangleindices_[0];
        int index2 = t.triangleindices_[1];
        int index3 = t.triangleindices_[2];
        for (int j=0;j<3;j++)
        {
            vec1[j] = vertices_[index2].position_[j] - vertices_[index1].position_[j];
            vec2[j] = vertices_[index3].position_[j] - vertices_[index1].position_[j];
        }

        Real3 N = LinAlg3x3::CrossProduct(vec1, vec2);

        volume_ += 1.0/6.0 * LinAlg3x3::DotProduct(g,N);
    }

    return volume_;
}

void Mesh::scaleVertices(Real num)
{
    #pragma omp parallel for
    for (int i=0;i<vertices_.size();i++)
    {
        for (int j=0;j<3;j++)
        {
            vertices_[i].position_[j] = num * vertices_[i].position_[j];
        }
    }

    update();
}

void Mesh::getVertexDistance(const Real3& v1, const Real3& v2, Real3& distVec, Real& dist)
{
    distVec.fill(0);
    distVec = v1 - v2;

    if (isPeriodic()){
        for (int i=0;i<3;i++){
            if (distVec[i] > boxLength_[i] * 0.5) distVec[i] -= boxLength_[i];
            if (distVec[i] < -boxLength_[i] * 0.5) distVec[i] += boxLength_[i];
        }
    }
    dist = LinAlg3x3::norm(distVec);
}

void Mesh::getVertexDistance(const vertex& v1, const vertex& v2, Real3& distVec, Real& dist)
{
    getVertexDistance(v1.position_, v2.position_, distVec, dist);
}

void Mesh::CalculateShift(const Real3& v1, const Real3& v2, Real3& shiftVec)
{
    for (int i=0;i<3;i++)
    {
        Real d = v1[i] - v2[i];

        if (d > (0.5 * boxLength_[i])) shiftVec[i] = -boxLength_[i];
        else if (d < (- 0.5 * boxLength_[i])) shiftVec[i] = boxLength_[i];
        else shiftVec[i] = 0.0;
    }
}

Mesh::Real3 Mesh::getShiftedVertexPosition(const vertex& v1, const vertex& v2)
{
    Real3 dist;
    Real distsq;

    // v1 - v2
    getVertexDistance(v1, v2, dist, distsq);
    Real3 ret = v2.position_ + dist;

    return ret;
}

Mesh::Real3 Mesh::getShiftIntoBox(const Real3& v1)
{
    Real3 shift;
    Real3 center;
    center = boxLength_ * 0.5;

    for (int i=0;i<3;i++)
    {
        Real diff = v1[i] - center[i];

        if (diff >= boxLength_[i] * 0.5) shift[i] = -boxLength_[i];
        else if (diff <= - boxLength_[i] * 0.5) shift[i] = boxLength_[i];
        else shift[i] = 0.0;
    }

    return shift;
}

// Compute per-vertex point areas
void Mesh::CalculateCornerArea()
{
    triangleArea_.clear();
    triangleArea_.resize(vertices_.size(),0.0);

    int nf = triangles_.size();

    cornerArea_.clear();
    cornerArea_.resize(nf);

	for (int i = 0; i < nf; i++) {
		// Edges
        auto& t = triangles_[i];
        int index1 = t.triangleindices_[0];
        int index2 = t.triangleindices_[1];
        int index3 = t.triangleindices_[2];

        Real3 edge1, edge2, edge3;
        Real edge1sq, edge2sq, edge3sq;

        getVertexDistance(vertices_[index3], vertices_[index2], edge1, edge1sq);
        getVertexDistance(vertices_[index1], vertices_[index3], edge2, edge2sq);
        getVertexDistance(vertices_[index2], vertices_[index1], edge3, edge3sq);

		// Compute corner weights
        Real3 cross = LinAlg3x3::CrossProduct(edge1, edge2);
        Real  area  = 0.5 * LinAlg3x3::norm(cross);

        Real e1 = LinAlg3x3::norm(edge1);
        Real e2 = LinAlg3x3::norm(edge2);
        Real e3 = LinAlg3x3::norm(edge3);

		Real3 l2 = { e1*e1, e2*e2, e3*e3};

		// Barycentric weights of circumcenter
		Real3 bcw = { l2[0] * (l2[1] + l2[2] - l2[0]),
		                 l2[1] * (l2[2] + l2[0] - l2[1]),
		                 l2[2] * (l2[0] + l2[1] - l2[2]) };

		if (bcw[0] <= 0.0f) {
			cornerArea_[i][1] = -0.25f * l2[2] * area/LinAlg3x3::DotProduct(edge1, edge3);
			cornerArea_[i][2] = -0.25f * l2[1] * area/LinAlg3x3::DotProduct(edge1, edge2);
			cornerArea_[i][0] = area - cornerArea_[i][1] - cornerArea_[i][2];
		} else if (bcw[1] <= 0.0f) {
			cornerArea_[i][2] = -0.25f * l2[0] * area/LinAlg3x3::DotProduct(edge2, edge1);
			cornerArea_[i][0] = -0.25f * l2[2] * area/LinAlg3x3::DotProduct(edge2, edge3);
			cornerArea_[i][1] = area - cornerArea_[i][2] - cornerArea_[i][0];
		} else if (bcw[2] <= 0.0f) {
			cornerArea_[i][0] = -0.25f * l2[1] * area/LinAlg3x3::DotProduct(edge3, edge2);
			cornerArea_[i][1] = -0.25f * l2[0] * area/LinAlg3x3::DotProduct(edge3, edge1);
			cornerArea_[i][2] = area - cornerArea_[i][0] - cornerArea_[i][1];
		} else {
			float scale = 0.5f * area / (bcw[0] + bcw[1] + bcw[2]);
			for (int j = 0; j < 3; j++)
            {
                int next = j - 1;
                int nextnext = j -2;

                if (next < 0)
                {
                    next += 3;
                }

                if (nextnext < 0)
                {
                    nextnext += 3;
                }

				cornerArea_[i][j] = scale * (bcw[next] +
				                             bcw[nextnext]);
            }
		}

		triangleArea_[t.triangleindices_[0]] += cornerArea_[i][0];
		triangleArea_[t.triangleindices_[1]] += cornerArea_[i][1];
		triangleArea_[t.triangleindices_[2]] += cornerArea_[i][2];
	}
}

void Mesh::update()
{
    for (int i=0;i<triangles_.size();i++)
    {
        auto& t = triangles_[i];

        for (int j=0;j<3;j++)
        {
            int index = t.triangleindices_[j];
            t.vertices_[j] = vertices_[index];
        }
    }
}

void Mesh::CalcVertexNormals()
{
    MeshTools::CalculateTriangleAreasAndFaceNormals(*this, triangleArea_, facetNormals_);

    // calculate vertex normals
    std::vector<Real3> VertNorms(getNumVertices(), {{0,0,0}});

    // calculate the normals of each vertices
    for (int i=0;i<triangles_.size();i++){
        auto& t = triangles_[i];
        Real3 normals = facetNormals_[i];

        for (int j=0;j<3;j++){
            int index = t.triangleindices_[j];

            for (int k=0;k<3;k++){
                VertNorms[index][k] += normals[k];
            }
        }
    }

    #pragma omp parallel for
    for (int i=0;i<getNumVertices();i++)
    {
        Real norm = LinAlg3x3::norm(VertNorms[i]);
        vertices_[i].normals_ = VertNorms[i]/norm;
    }
}

void Mesh::CalcVertexNormalsAreaWeighted()
{
    // calculate vertex normals
    vertexNormals_.resize(vertices_.size());
    Real3 zeroArr = {{0,0,0}};
    std::fill(vertexNormals_.begin(), vertexNormals_.end(), zeroArr);

    // calculate the normals of each vertices
    for (int i=0;i<triangles_.size();i++)
    {
        auto& t = triangles_[i];
        Real area = triangleArea_[i];
        Real3 normals = facetNormals_[i];

        for (int j=0;j<3;j++)
        {
            int index = t.triangleindices_[j];
            auto& vNorm = vertexNormals_[index];

            for (int k=0;k<3;k++)
            {
                vNorm[k] += area*normals[k];
            }
        }
    }


    for (int i=0;i<vertexNormals_.size();i++)
    {
        Real norm = LinAlg3x3::norm(vertexNormals_[i]);

        for (int j=0;j<3;j++)
        {
            vertexNormals_[i][j] = vertexNormals_[i][j]/norm;
        }
    }

    for (int i=0;i<vertices_.size();i++)
    {
        vertices_[i].normals_ = vertexNormals_[i];
    }
}

std::vector<Mesh::Real3> Mesh::getVertexPositions()
{
    std::vector<Real3> v;
    for (int i=0;i<vertices_.size();i++)
    {
        v.push_back(vertices_[i].position_);
    }
    return v;
}

                                    /**************************************
                                     ************   MeshTools *************
                                     *************************************/       
void MeshTools::writePLY(std::string filename, Mesh& mesh)
{
    std::vector<Real3> verts;
    std::vector<INT3> faces;

    const auto& v = mesh.getvertices();
    const auto& f = mesh.gettriangles();

    verts.resize(v.size());
    faces.resize(f.size());

    for (int i=0;i<v.size();i++)
    {
        verts[i] = v[i].position_;
    }

    for (int i=0;i<f.size();i++)
    {
        faces[i] = f[i].triangleindices_;
    }

    writePLY(filename, verts, faces);
}

void MeshTools::writePLY(std::string filename, const std::vector<Real3>& vertices, const std::vector<INT3>& faces, const std::vector<Real3>& normals)
{
    std::ofstream ofs;
    ofs.open(filename);

    ofs << "ply" << "\n";
    ofs << "format ascii 1.0\n";
    ofs << "comment Created by Yusheng Cai\n";

    int sizeVertex = vertices.size();
    int sizetriangle = faces.size();

    ofs << "element vertex " << sizeVertex << std::endl; 
    ofs << "property float x\n";
    ofs << "property float y\n";
    ofs << "property float z\n";
    ofs << "property float nx\n";
    ofs << "property float ny\n";
    ofs << "property float nz\n";
    ofs << "element face " << sizetriangle << "\n";
    ofs << "property list uchar uint vertex_indices\n";
    ofs << "end_header\n";

    ofs << std::fixed << std::setprecision(6);
    for (int i=0;i<sizeVertex;i++)
    {
        for (int j=0;j<3;j++)
        {
            ofs << vertices[i][j] << " ";
        }

        for (int j=0;j<3;j++)
        {
            ofs << normals[i][j] << " ";
        }

        ofs << "\n";
    }

    for (int i=0;i<sizetriangle;i++)
    {
        ofs << 3 << " ";
        auto& t = faces[i];

        for (int j=0;j<3;j++)
        {
            ofs << t[j] << " ";
        }
        ofs << "\n";
    }

    ofs.close();
}


void MeshTools::writePLY(std::string filename, const std::vector<Real3>& vertices, const std::vector<INT3>& faces, \
Real factor)
{
    std::ofstream ofs;
    ofs.open(filename);

    ofs << "ply" << "\n";
    ofs << "format ascii 1.0\n";
    ofs << "comment Created by Yusheng Cai\n";

    int sizeVertex = vertices.size();
    int sizetriangle = faces.size();

    ofs << "element vertex " << sizeVertex << std::endl; 
    ofs << "property float x\n";
    ofs << "property float y\n";
    ofs << "property float z\n";
    ofs << "element face " << sizetriangle << "\n";
    ofs << "property list uchar uint vertex_indices\n";
    ofs << "end_header\n";

    ofs << std::fixed << std::setprecision(6);
    for (int i=0;i<sizeVertex;i++)
    {
        for (int j=0;j<3;j++)
        {
            ofs << vertices[i][j] * factor << " ";
        }

        ofs << "\n";
    }

    for (int i=0;i<sizetriangle;i++)
    {
        ofs << 3 << " ";
        auto& t = faces[i];

        for (int j=0;j<3;j++)
        {
            ofs << t[j] << " ";
        }
        ofs << "\n";
    }

    ofs.close();
}

void MeshTools::writePLYRGB(std::string filename, const std::vector<Real3>& vertices, const std::vector<INT3>& faces, const std::vector<Real3>& RGB)
{
    std::ofstream ofs;
    ofs.open(filename);

    ofs << "ply" << "\n";
    ofs << "format ascii 1.0\n";
    ofs << "comment Created by Yusheng Cai\n";

    int sizeVertex = vertices.size();
    int sizetriangle = faces.size();
    int sizeRGB    = RGB.size();

    ofs << "element vertex " << sizeVertex << std::endl; 
    ofs << "property float x\n";
    ofs << "property float y\n";
    ofs << "property float z\n";
    ofs << "property uchar red\n";
    ofs << "property uchar green\n";
    ofs << "property uchar blue\n";
    ofs << "element face " << sizetriangle << "\n";
    ofs << "property list uchar uint vertex_indices\n";
    ofs << "end_header\n";

    ASSERT((sizeRGB == sizeVertex), "The size of RGB value provided must agree with the size of vertex.");

    ofs << std::fixed << std::setprecision(6);
    for (int i=0;i<sizeVertex;i++)
    {
        for (int j=0;j<3;j++)
        {
            ofs << vertices[i][j] << " ";
        }

        for (int j=0;j<3;j++)
        {
            ofs << (int)RGB[i][j] << " ";
        }

        ofs << "\n";
    }

    for (int i=0;i<sizetriangle;i++)
    {
        ofs << 3 << " ";
        auto& t = faces[i];

        for (int j=0;j<3;j++)
        {
            ofs << t[j] << " ";
        }
        ofs << "\n";
    }

    ofs.close();
}

void MeshTools::writeSTL(std::string name, const std::vector<Real3>& vertices, const std::vector<INT3>& faces){

}

void MeshTools::writeSTL(std::string name, Mesh& mesh){
    std::ofstream ofs;
    ofs.open(name);

    const auto& f = mesh.gettriangles();
    const auto& n = mesh.getFaceNormals();
    
    ofs << "solid " << name << "\n";

    for (int i=0;i<f.size();i++)
    {
        ofs << "facet normal " << n[i][0] << " " << n[i][1] << " " << n[i][2] << "\n";

        ofs << "\touter loop\n";
        auto& t = f[i];

        for (int j=0;j<3;j++)
        {
            auto& pos = t.vertices_[j].position_;
            ofs << "vertex " << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
        }

        ofs << "\tendloop\n";
        ofs << "endfacet\n";
    }
    ofs << "endsolid " << name;

}

bool MeshTools::readPLYlibr(std::string& filename, Mesh& mesh)
{
    happly::PLYData plydata(filename);
    std::vector<std::array<double,3>> vPos = plydata.getVertexPositions();
    std::vector<std::vector<size_t>> fInd = plydata.getFaceIndices<size_t>();
    std::vector<std::string> normalNames = {"nx", "ny", "nz"};
    std::vector<bool> hasProperty_(3, false);
    std::vector<Real3> normals_;
    normals_.resize(vPos.size());

    auto names = plydata.getElement("vertex").getPropertyNames();
    bool hasnormals=false;
    for (int i=0;i<3;i++){
        hasProperty_[i] = std::find(names.begin(), names.end(), normalNames[i]) != names.end();
    }

    if (hasProperty_[0] && hasProperty_[1] && hasProperty_[2]){
        hasnormals=true;
    }

    // do something when it does have normals
    if (hasnormals)
    {
        for (int i=0;i<3;i++){
            const auto& n = plydata.getElement("vertex").getProperty<double>(normalNames[i]);
            ASSERT((n.size() == vPos.size()), "The size of nx must equal to the number of vertices.");

            for (int j=0;j<vPos.size();j++){normals_[j][i] = n[j];}
        }
    }

    // first we update the mesh
    auto& vertices = mesh.accessvertices();
    auto& triangles= mesh.accesstriangles();
    vertices.clear();
    triangles.clear();

    vertices.resize(vPos.size());
    triangles.resize(fInd.size());

    for (int i=0;i<vertices.size();i++)
    {
        auto& v = vertices[i];
        for (int j=0;j<3;j++){
            v.position_[j] = vPos[i][j];
            v.normals_[j]  = normals_[i][j];
        }
    }

    for (int i=0;i<triangles.size();i++){
        auto& t = triangles[i];

        ASSERT((fInd[i].size() == 3), "We are reading triangles here while the provided face is " << fInd[i].size());

        for (int j=0;j<3;j++){
            t.triangleindices_[j] = fInd[i][j];
            t.vertices_[j] = vertices[fInd[i][j]];
        }
    }

    if ( ! hasnormals){
        std::cout << "Calculating normals by myself." << std::endl;
        mesh.CalcVertexNormals();

        mesh.update();
    }

    return true;
}

bool MeshTools::readPLY(std::string& filename, Mesh& mesh)
{
    auto& vertices = mesh.accessvertices();
    auto& triangles= mesh.accesstriangles();
    vertices.clear();
    triangles.clear();

    // open the file
    std::ifstream ifs;
    std::stringstream ss;
    ifs.open(filename);

    if (! ifs.is_open()){
        return false;
    }

    std::string sentence;
    int numfaces;
    int numvertex;
    while (std::getline(ifs, sentence)){
        ss.str(sentence);
        std::string token;

        std::vector<std::string> vectorstring;
        while (ss >> token){
            vectorstring.push_back(token);
        }

        if (vectorstring[0] == "end_header"){
            break;
        }

        if (vectorstring[0] == "format"){
            ASSERT((vectorstring[1] == "ascii"), "Currently do not support non-ascii format ply files.");
        }

        if (vectorstring[0] == "element"){
            ASSERT((vectorstring.size() == 3), "The element can only have 3.");
            if (vectorstring[1] == "vertex"){
                numvertex = StringTools::StringToType<int>(vectorstring[2]);
            }

            if (vectorstring[1] == "face"){
                numfaces = StringTools::StringToType<int>(vectorstring[2]);
            }
        }

        ss.clear();
    }

    ss.clear();
    // read in the vertices as well as their normals     
    std::string datasentence;
    for (int i=0;i<numvertex;i++)
    {
        std::getline(ifs, datasentence);
        ss.str(datasentence);
        std::string data;
        std::vector<std::string> vectordata;

        vertex v;

        while (ss >> data)
        {
            vectordata.push_back(data);
        }

        ASSERT((vectordata.size() == 6), "The ply file must contain x, y, z, nx, ny ,nz");

        Real3 position;
        Real3 normals;

        for (int j=0;j<3;j++)
        {
            position[j] = StringTools::StringToType<Real>(vectordata[j]);
            normals[j] = StringTools::StringToType<Real>(vectordata[j+3]);
        }

        v.position_ = position;
        v.normals_ = normals;

        vertices.push_back(v);

        ss.clear();
    }

    ss.clear();
    // read in the triangles
    std::string trianglesentence_;
    for (int i=0;i<numfaces;i++)
    {
        std::getline(ifs,trianglesentence_);
        ss.str(trianglesentence_);
        std::string data;
        std::vector<std::string> vectordata;

        while (ss >> data){
            vectordata.push_back(data);
        }

        ASSERT((vectordata.size() == 4), "The triangle file must contain 3 f1 f2 f3");

        triangle t;

        INT3 faceid;

        for (int j=0;j<3;j++){faceid[j] = StringTools::StringToType<int>(vectordata[j+1]);}

        t.triangleindices_ = faceid;
        for (int j=0;j<3;j++){
            t.vertices_[j] = vertices[faceid[j]];
        }

        triangles.push_back(t);
        ss.clear();
    }

    ifs.close();

    return true;
}

MeshTools::Real3 MeshTools::calculateShift(const Real3& vec1, const Real3& vec2, const Real3& boxLength)
{
    Real3 diff;
    for (int i=0;i<3;i++)
    {
        Real d = vec1[i] - vec2[i];

        if (d > (0.5 * boxLength[i])) diff[i] = -boxLength[i];
        else if (d < (- 0.5 * boxLength[i])) diff[i] = boxLength[i];
        else diff[i] = 0.0;
    }
    
    return diff;
}

void MeshTools::calculateDistance(const Real3& vec1, const Real3& vec2, const Real3& boxLength, Real3& distance, Real& distsq)
{
    distsq = 0.0;
    for (int i=0;i<3;i++)
    {
        Real d = vec1[i] - vec2[i];

        if (d > (0.5 * boxLength[i])) distance[i] = d-boxLength[i];
        else if (d < (- 0.5 * boxLength[i])) distance[i] = d+boxLength[i];
        else distance[i] = d;

        distsq += distance[i] * distance[i];
    }
}

bool MeshTools::isPeriodicEdge(const Real3& vec1, const Real3& vec2, Real3& newarr, const Real3& boxLength)
{
    // calculate pbc shift b/t vec1 and vec2 in some box with boxlength
    Real3 diff = calculateShift(vec1, vec2, boxLength);

    // shifted vector 
    newarr = vec1 + diff;

    // check if this is periodic edge
    for (int i=0;i<3;i++){
        if (std::abs(diff[i]) > (boxLength[i] * 0.5)){
            return true;
        }
    }

    return false;
}

void MeshTools::MapVerticesToFaces(Mesh& mesh, std::vector<std::vector<int>>& map)
{
    auto& vertices = mesh.getvertices();
    auto& triangles= mesh.gettriangles();

    map.clear();
    map.resize(vertices.size());

    for (int i=0;i<triangles.size();i++)
    {
        auto t = triangles[i].triangleindices_;

        for (int j=0;j<3;j++)
        {
            map[t[j]].push_back(i);
        }
    }
}

void MeshTools::CalculateTriangleAreasAndFaceNormals(Mesh& mesh, std::vector<Real>& Areas, std::vector<Real3>& Normals)
{
    auto& triangles = mesh.accesstriangles();
    auto& vertices  = mesh.accessvertices();

    Areas.clear();
    Areas.resize(triangles.size());

    Normals.clear();
    Normals.resize(triangles.size());

    // calculate the area of the triangles as well as the normals of the faces of triangles
    #pragma omp parallel for
    for (int i=0;i<triangles.size();i++){
        auto& t = triangles[i];

        int index1 = t.triangleindices_[0];
        int index2 = t.triangleindices_[1];
        int index3 = t.triangleindices_[2];

        Real3 diff1 = {};
        Real3 diff2 = {};
        Real3 diff3 = {};
        Real norm1, norm2, norm3;

        mesh.getVertexDistance(vertices[index1], vertices[index2], diff1, norm1);
        mesh.getVertexDistance(vertices[index3], vertices[index2], diff2, norm2);
        mesh.getVertexDistance(vertices[index3], vertices[index1], diff3, norm3);

        Real3 crossProduct = LinAlg3x3::CrossProduct(diff1, diff2);
        Real norm = LinAlg3x3::norm(crossProduct);
        Real a    = norm*0.5;
        Real3 n   = crossProduct/norm;

        Normals[i] = n;
        Areas[i]   = a;
    }
}

void MeshTools::CalculateVertexNeighbors(Mesh& mesh, std::vector<std::vector<int>>& neighborIndices)
{
    auto& vertices = mesh.getvertices();
    auto& triangles= mesh.gettriangles();

    neighborIndices.clear();
    neighborIndices.resize(vertices.size());

    for (int i=0;i<triangles.size();i++){
        auto& t = triangles[i];

        for (int j=0;j<3;j++){
            int index1 = t.triangleindices_[j];
            for (int k=0;k<3;k++){
                if (j != k){
                    int index2 = t.triangleindices_[k];
                    auto it = std::find(neighborIndices[index1].begin(), neighborIndices[index1].end(), index2);
                    if (it == neighborIndices[index1].end()){
                        neighborIndices[index1].push_back(index2);
                    }
                }
            }
        }
    }
}

void MeshTools::MapEdgeToFace(Mesh& mesh, std::map<INT2, std::vector<int>>& mapEdgeToFace, std::vector<std::vector<INT2>>& mapVertexToEdge)
{
    // get the triangles and vertices from mesh object
    const auto& triangles = mesh.gettriangles();
    const auto& vertices  = mesh.getvertices();

    // clear the inputs 
    mapVertexToEdge.clear();
    mapVertexToEdge.resize(vertices.size());
    mapEdgeToFace.clear();

    // iterate over all the triangles
    for (int i=0;i<triangles.size();i++){
        // find the triangles indices
        auto& t = triangles[i];

        for (int j=0;j<3;j++){
            // get the 2 adjacent indices
            int idx1 = t.triangleindices_[j];
            int idx2 = t.triangleindices_[(j+1)%3];

            // find the min/max of the 2
            int minIndex = std::min(idx1, idx2);
            int maxIndex = std::max(idx1, idx2);
            INT2 arr = {{minIndex, maxIndex}};

            // insert into map
            auto it = mapEdgeToFace.find(arr);
            if (it == mapEdgeToFace.end()){
                std::vector<int> faceIndices;
                faceIndices.push_back(i);

                mapEdgeToFace.insert(std::make_pair(arr,faceIndices));
            }
            else{
                it -> second.push_back(i);
            }
        }

        // map vertex to edges 
        for (int j=0;j<3;j++){
            int idx1 = t.triangleindices_[j];
            int idx2 = t.triangleindices_[(j+1)%3];
            int minIdx = std::min(idx1, idx2);
            int maxIdx = std::max(idx1, idx2);

            INT2 indices = {{minIdx, maxIdx}};

            for (int k=0;k<2;k++)
            {
                // find if the edge is already in the vector
                auto f = std::find(mapVertexToEdge[indices[k]].begin(), mapVertexToEdge[indices[k]].end(), indices);

                if (f == mapVertexToEdge[indices[k]].end())
                {
                    mapVertexToEdge[indices[k]].push_back(indices);
                }
            }
        }
    }

    // do a check to make sure that each edge is shared by at least 1 and at most 2 faces
    int id = 0;
    for (auto it=mapEdgeToFace.begin(); it != mapEdgeToFace.end();it++){
        if (it -> second.size() > 2){
            std::cout << "This is for edge " << id << " consisting of vertex " << it ->first[0] <<" " << it->first[1] <<"\n";
            for (auto s : it -> second)
            {
                std::cout << "Face it corresponds to : " << s <<"\n";
            }
        }
        ASSERT((it->second.size() == 1 || it -> second.size() ==2), "The number of faces shared by an edge needs to be either 1 or 2 \
        , however the calculated result shows " << it -> second.size());
        id ++;
    }

    return;
}

void MeshTools::CalculateBoundaryVertices(Mesh& mesh, std::map<INT2, std::vector<int>>& mapEdgeToFace, std::vector<bool>& boundaryIndicator)
{
    const auto& vertices = mesh.getvertices();

    boundaryIndicator.clear();
    boundaryIndicator.resize(vertices.size(), false);

    for (auto it = mapEdgeToFace.begin(); it != mapEdgeToFace.end(); it ++){
        int size = it -> second.size();
        ASSERT(( size == 1 || size == 2), "An edge can only be shared by 1 or 2 faces.");
        INT2 arr = it -> first;

        // if the edge size == 1, then it is a boundary edge
        if (size == 1){
            boundaryIndicator[arr[0]] = true;
            boundaryIndicator[arr[1]] = true;
        }
    }
}

bool MeshTools::IsBoundary(int index, const std::vector<bool>& boundaryIndicator)
{
    return boundaryIndicator[index];
}

void MeshTools::ConvertToNonPBCMesh(Mesh& mesh, std::vector<Real3>& vertices, std::vector<INT3>& triangles)
{
    // if periodic, then do something, else do nothing 
    if (mesh.isPeriodic()){
        auto MeshVertices = mesh.getvertices();
        auto MeshTriangles = mesh.gettriangles();

        // first let's copy all the vertices 
        vertices.clear();
        for (auto& v : MeshVertices){
            vertices.push_back(v.position_);
        }

        // check if triangle is periodic 
        for (auto& t : MeshTriangles){
            // find all the edge lengths
            bool periodicTriangle = MeshTools::IsPeriodicTriangle(MeshVertices, t.triangleindices_, mesh.getBoxLength());

            // if this particular triangle is not periodic 
            if (! periodicTriangle){
                triangles.push_back(t.triangleindices_);
            }
            // if it's periodic triangle, then we push back 3 new vertices 
            else{
                Real3 verticesNew1;
                Real3 verticesNew2, verticesDiff2;
                Real distsq2;
                Real3 verticesNew3, verticesDiff3;
                Real distsq3;
                int idx1 = t.triangleindices_[0];
                int idx2 = t.triangleindices_[1];
                int idx3 = t.triangleindices_[2];

                // get the pbc corrected distance 
                mesh.getVertexDistance(MeshVertices[idx2].position_, MeshVertices[idx1].position_,verticesDiff2, distsq2);
                mesh.getVertexDistance(MeshVertices[idx3].position_, MeshVertices[idx1].position_,verticesDiff3, distsq3); 

                // get the new vertices --> with respect to position 1
                verticesNew2 = MeshVertices[idx1].position_ + verticesDiff2;
                verticesNew3 = MeshVertices[idx1].position_ + verticesDiff3;

                // Find approximately the center of the triangle
                Real3 center_of_triangle = (MeshVertices[idx1].position_ + verticesNew2 + verticesNew3) * (1.0/3.0);
                Real3 shift = mesh.getShiftIntoBox(center_of_triangle);

                // update the vertices
                verticesNew1 = MeshVertices[idx1].position_ + shift;
                verticesNew2 = verticesNew2 + shift;
                verticesNew3 = verticesNew3 + shift;

                int NewIndex1 = vertices.size();
                vertices.push_back(verticesNew1);
                int NewIndex2 = vertices.size();
                vertices.push_back(verticesNew2);
                int NewIndex3 = vertices.size();
                vertices.push_back(verticesNew3);

                INT3 NewT = {{NewIndex1, NewIndex2, NewIndex3}};
                triangles.push_back(NewT);
            }
        }
    }
}

bool MeshTools::IsPeriodicTriangle(std::vector<vertex>& Vertices,INT3& face, Real3 BoxLength)
{
    bool IsPeriodic=false;
    for (int i=0;i<3;i++){
        int index1 = face[i];
        int index2 = face[(i+1) % 3]; 
        Real3 diff;
        for (int j=0;j<3;j++){
            diff[j] = Vertices[index1].position_[j] - Vertices[index2].position_[j];

            if (std::abs(diff[j]) >= 0.5 * BoxLength[j]){
                return true;
            }
        }
    }

    return false;
}

bool MeshTools::IsPeriodicTriangle(const Mesh& mesh, int faceIndex){
    const auto& f = mesh.gettriangles();
    const auto& v = mesh.getvertices();
    Real3 boxLength = mesh.getBoxLength();
    INT3 t = f[faceIndex].triangleindices_;
    for (int i=0;i<3;i++){
        int index1 = t[i];
        int index2 = t[(i+1) % 3];
        Real3 diff;
        for (int j=0;j<3;j++){
            diff[j] = v[index1].position_[j] - v[index2].position_[j];

            if (std::abs(diff[j] >= 0.5 * boxLength[j])){return true;}
        }
    }

    return false;
}

void MeshTools::ShiftPeriodicTriangle(const std::vector<vertex>& Vertices, const INT3& face,\
                                        Real3 BoxLength, Real3& A, Real3& B, Real3& C){
    // half box
    Real3 half_box = BoxLength * 0.5;

    // first shift B C wrt to A
    A = Vertices[face[0]].position_;
    B = Vertices[face[1]].position_;
    C = Vertices[face[2]].position_;
    Real3 shiftB = MeshTools::calculateShift(B, A, BoxLength);
    Real3 shiftC = MeshTools::calculateShift(C, A, BoxLength);
    B = B + shiftB;
    C = C + shiftC;

    // calculate the triangle center
    Real3 center = (A + B + C) * 1.0/3.0;

    // calculate shift
    Real3 shiftTriangle = MeshTools::calculateShift(center, half_box, BoxLength);

    A = A + shiftTriangle;
    B = B + shiftTriangle;
    C = C + shiftTriangle;
}

MeshTools::INT2 MeshTools::makeEdge(int i, int j)
{
    int minIndex = std::min(i,j);
    int maxIndex = std::max(i,j);

    INT2 ret = {{minIndex, maxIndex}};

    return ret;
}

void MeshTools::MapEdgeToOpposingVertices(Mesh& mesh, std::map<INT2, std::vector<int>>& mapEdgeToFace, std::map<INT2, std::vector<int>>& MapEdgeToOppoVertices)
{
    const auto& tri = mesh.gettriangles();

    // clear the map at which we want to write to  
    MapEdgeToOppoVertices.clear();

    // iterate over the edges 
    for (auto it = mapEdgeToFace.begin(); it != mapEdgeToFace.end(); it ++){
        std::vector<int> faces = it -> second;
        INT2 edge = it -> first;
        std::vector<int> opposingPoints;

        int numf = faces.size();

        // only perform calculation if we are working with non boundary edge --> (edges shared by 2 faces)
        if (numf == 2){
            // iterate over the 2 faces 
            for (int f : faces){
                auto& TriIndices = tri[f].triangleindices_;

                // iterate over the triangular indices of each of the face 
                for (int id : TriIndices){
                    bool IsInEdge = std::find(edge.begin(), edge.end(), id) != edge.end();

                    if ( ! IsInEdge){
                        opposingPoints.push_back(id);
                    }
                }
            }

            ASSERT((opposingPoints.size() == 2), "The opposing points of an edge needs to be 2 while it is " << opposingPoints.size());
            MapEdgeToOppoVertices.insert(std::make_pair(edge, opposingPoints));
        }
    }
}

void MeshTools::CutMesh(Mesh& mesh, Real3 volume){
    const auto& verts = mesh.getvertices();
    const auto& tri   = mesh.gettriangles();

    std::vector<vertex> vertices;
    std::vector<triangle> faces;

    int index=0;
    std::map<int,int> MapOldIndexToNew;
    for (int i=0;i<verts.size();i++){
        auto& v = verts[i];
        if (v.position_[0] >= volume[0] && v.position_[1] >= volume[1] && v.position_[2] >= volume[2]){
            int newindex = index;
            vertices.push_back(v);
            MapOldIndexToNew.insert(std::make_pair(i, newindex));
            index ++;
        }
    }

    for (auto& t : tri){
        bool it1 = MapOldIndexToNew.find(t.triangleindices_[0]) != MapOldIndexToNew.end();
        bool it2 = MapOldIndexToNew.find(t.triangleindices_[1]) != MapOldIndexToNew.end();
        bool it3 = MapOldIndexToNew.find(t.triangleindices_[2]) != MapOldIndexToNew.end();

        if (it1 && it2 && it3){
            triangle newt;
            for (int i=0;i<3;i++){
                newt[i] = MapOldIndexToNew[t[i]];
            }
            faces.push_back(newt);
        }
    }

    auto& v = mesh.accessvertices();    
    auto& f = mesh.accesstriangles();  
    v.clear();
    v.insert(v.end(),vertices.begin(), vertices.end());
    f.clear();
    f.insert(f.end(), faces.begin(), faces.end());
    
}

void MeshTools::CutMesh(Mesh& mesh, std::vector<INT3>& faces, std::vector<Real3>& vertices, Real3 volume)
{
    vertices.clear();
    faces.clear();

    const auto& verts = mesh.getvertices();
    const auto& tri   = mesh.gettriangles();

    int index=0;
    std::map<int,int> MapOldIndexToNew;
    for (int i=0;i<verts.size();i++){
        auto& v = verts[i];
        if (v.position_[0] >= volume[0] && v.position_[1] >= volume[1] && v.position_[2] >= volume[2]){
            int newindex = index;
            vertices.push_back(v.position_);
            MapOldIndexToNew.insert(std::make_pair(i, newindex));
            index ++;
        }
    }

    for (auto& t : tri){
        bool it1 = MapOldIndexToNew.find(t.triangleindices_[0]) != MapOldIndexToNew.end();
        bool it2 = MapOldIndexToNew.find(t.triangleindices_[1]) != MapOldIndexToNew.end();
        bool it3 = MapOldIndexToNew.find(t.triangleindices_[2]) != MapOldIndexToNew.end();

        if (it1 && it2 && it3){
            INT3 newt;
            for (int i=0;i<3;i++){
                auto it = MapOldIndexToNew.find(t.triangleindices_[i]);
                int newIndex = it -> second;
                newt[i] = newIndex;
            }
            faces.push_back(newt);
        }
    }
}

void MeshTools::CalculateCornerArea(Mesh& mesh, std::vector<Real3>& CornerArea, std::vector<Real>& VertexArea)
{
    const auto& vertices = mesh.getvertices();
    const auto& triangles= mesh.gettriangles();
    int nf = triangles.size();
    int nv = vertices.size();

    // resize corner area
    CornerArea.clear();
    CornerArea.resize(nf);

    // resize vertex area
    VertexArea.clear();
    VertexArea.resize(nv,0.0);

	for (int i = 0; i < nf; i++) {
		// Edges
        auto& t = triangles[i];
        int index1 = t.triangleindices_[0];
        int index2 = t.triangleindices_[1];
        int index3 = t.triangleindices_[2];

        Real3 edge1, edge2, edge3;
        Real edge1sq, edge2sq, edge3sq;

        mesh.getVertexDistance(vertices[index3], vertices[index2], edge1, edge1sq);
        mesh.getVertexDistance(vertices[index1], vertices[index3], edge2, edge2sq);
        mesh.getVertexDistance(vertices[index2], vertices[index1], edge3, edge3sq);

		// Compute corner weights
        Real3 cross = LinAlg3x3::CrossProduct(edge1, edge2);
        Real  area  = 0.5 * LinAlg3x3::norm(cross);

        Real e1 = LinAlg3x3::norm(edge1);
        Real e2 = LinAlg3x3::norm(edge2);
        Real e3 = LinAlg3x3::norm(edge3);

		Real3 l2 = { e1*e1, e2*e2, e3*e3};

		// Barycentric weights of circumcenter
		Real3 bcw = { l2[0] * (l2[1] + l2[2] - l2[0]),
		                 l2[1] * (l2[2] + l2[0] - l2[1]),
		                 l2[2] * (l2[0] + l2[1] - l2[2]) };

		if (bcw[0] <= 0.0f) {
			CornerArea[i][1] = -0.25f * l2[2] * area/LinAlg3x3::DotProduct(edge1, edge3);
			CornerArea[i][2] = -0.25f * l2[1] * area/LinAlg3x3::DotProduct(edge1, edge2);
			CornerArea[i][0] = area - CornerArea[i][1] - CornerArea[i][2];
		} else if (bcw[1] <= 0.0f) {
			CornerArea[i][2] = -0.25f * l2[0] * area/LinAlg3x3::DotProduct(edge2, edge1);
			CornerArea[i][0] = -0.25f * l2[2] * area/LinAlg3x3::DotProduct(edge2, edge3);
			CornerArea[i][1] = area - CornerArea[i][2] - CornerArea[i][0];
		} else if (bcw[2] <= 0.0f) {
			CornerArea[i][0] = -0.25f * l2[1] * area/LinAlg3x3::DotProduct(edge3, edge2);
			CornerArea[i][1] = -0.25f * l2[0] * area/LinAlg3x3::DotProduct(edge3, edge1);
			CornerArea[i][2] = area - CornerArea[i][0] - CornerArea[i][1];
		} else {
			float scale = 0.5f * area / (bcw[0] + bcw[1] + bcw[2]);
			for (int j = 0; j < 3; j++){
                int next = j - 1;
                int nextnext = j -2;

                if (next < 0){
                    next += 3;
                }

                if (nextnext < 0){
                    nextnext += 3;
                }

				CornerArea[i][j] = scale * (bcw[next] +
				                             bcw[nextnext]);
            }
		}

		VertexArea[t.triangleindices_[0]] += CornerArea[i][0];
		VertexArea[t.triangleindices_[1]] += CornerArea[i][1];
		VertexArea[t.triangleindices_[2]] += CornerArea[i][2];
	}
}

bool MeshTools::MTRayTriangleIntersection(Real3& A, Real3& B, Real3& C, Real3& O, Real3& D, Real& t, Real& u, Real& v)
{
    Real epsilon = 1e-8;
    
    // declare the variables 
    Real3 T, E1, E2, P, Q;
    Real det, invdet;

    // find E1=B-A and E2=C-A
    E1 = B - A;
    E2 = C - A;

    P = LinAlg3x3::CrossProduct(D, E2);
    det= LinAlg3x3::DotProduct(P, E1);

    // if det is close to 0, then the ray and triangle are parallel
    if (std::abs(det) < epsilon){
        return false;
    }

    // calculate inverse of det  
    invdet = 1.0 / det;
    T = O - A;

    // find u
    u = LinAlg3x3::DotProduct(T, P) * invdet;
    if (u < 0 || u > 1) return false;

    // find v
    Q = LinAlg3x3::CrossProduct(T,E1);
    v = LinAlg3x3::DotProduct(D, Q) * invdet;
    if (v < 0 || u+v > 1) return false;

    // find t
    t = LinAlg3x3::DotProduct(E2, Q) * invdet;

    return true;
}

void MeshTools::writeNonPBCMesh(std::string name, Mesh& mesh)
{
    if (mesh.isPeriodic())
    {
        std::vector<Real3> tempVertices;
        std::vector<INT3> tempFaces;
        MeshTools::ConvertToNonPBCMesh(mesh, tempVertices, tempFaces);

        MeshTools::writePLY(name, tempVertices, tempFaces);
    }
    else
    {
        std::cout << "WARNING: You asked to print NON PBC Mesh while the Mesh is not PBC." << "\n";
    }
}

void MeshTools::writeNonPeriodicTriangleIndices(std::string name, Mesh& mesh)
{
    std::ofstream ofs;
    ofs.open(name);
    std::vector<int> NonPeriodicTriangleIndices;

    const auto& f = mesh.gettriangles();
    const auto& v = mesh.getvertices();

    if (mesh.isPeriodic()){
        for (int i =0;i<f.size();i++){
            auto t = f[i].triangleindices_;
            if (! MeshTools::IsPeriodicTriangle(const_cast<std::vector<vertex>&>(v), t, mesh.getBoxLength())){
                NonPeriodicTriangleIndices.push_back(i);
            }
        }
    }
    else{
        NonPeriodicTriangleIndices.resize(f.size());
        std::iota(NonPeriodicTriangleIndices.begin(), NonPeriodicTriangleIndices.end(),0);
    }

    for (int index : NonPeriodicTriangleIndices){
        ofs << index << " ";
    }

    ofs.close();
}

void MeshTools::writeMeshArea(std::string filename, Mesh& mesh)
{
    // calculate the area and facet normals 
    std::vector<Real> triangleA;
    std::vector<Real3> facetN;
    MeshTools::CalculateTriangleAreasAndFaceNormals(mesh, triangleA, facetN);

    std::ofstream ofs;
    ofs.open(filename);

    for (int i=0;i<triangleA.size();i++){
        ofs << triangleA[i] << "\n";
    }

    ofs.close();
}

void MeshTools::writeCuttedMesh(std::string filename, Mesh& mesh, Real3& volume)
{
    std::vector<INT3> faces;
    std::vector<Real3> verts;
    CutMesh(mesh, faces, verts, volume);
    writePLY(filename, verts, faces);
}



void MeshTools::CheckDegenerateTriangle(Mesh& mesh, \
                                                    std::vector<int>& MergeFaces, \
                                                    std::vector<INT2>& MergeVertices)
{
    const auto& f = mesh.gettriangles();
    const auto& v = mesh.getvertices();
    MergeFaces.clear();
    Real epsilon=1e-8;
    Real merge_epsilon=1e-4;

    // what we check is if one side is whether or not a + b = c 
    for (int i=0;i<f.size();i++){
        INT3 indices = f[i].triangleindices_;

        // get all three vertices of the triangle
        Real3 A = v[indices[0]].position_;
        Real3 B = v[indices[1]].position_;
        Real3 C = v[indices[2]].position_;

        // obtain the side length
        Real3 vec;
        Real3 triangle_length;
        Real AB, BC, CA;
        mesh.getVertexDistance(A, B, vec, AB);
        mesh.getVertexDistance(B, C, vec, BC);
        mesh.getVertexDistance(C, A, vec, CA);

        // AB + BC > CA
        // BC + CA > AB
        // CA + AB > BC
        Real diff1 = AB + BC - CA;
        Real diff2 = BC + CA - AB;
        Real diff3 = CA + AB - BC;
 
        //ASSERT((diff1 >= 0 && diff2 >= 0 && diff3 >= 0), "Something weird with triangles.");

        if ((AB < merge_epsilon) || (BC < merge_epsilon) || (CA < merge_epsilon)){
            MergeFaces.push_back(i);
            INT2 verts_ind1, verts_ind2,  verts_ind3;

            verts_ind1 = {{indices[0], indices[1]}};
            verts_ind2 = {{indices[1], indices[2]}};
            verts_ind3 = {{indices[2], indices[0]}};
            Algorithm::sort(verts_ind1); Algorithm::sort(verts_ind2); Algorithm::sort(verts_ind3);

            if ((! Algorithm::contain(MergeVertices, verts_ind1)) && (AB<merge_epsilon)){MergeVertices.push_back(verts_ind1);}
            if ((! Algorithm::contain(MergeVertices, verts_ind2)) && (BC<merge_epsilon)){MergeVertices.push_back(verts_ind2);}
            if ((! Algorithm::contain(MergeVertices, verts_ind3)) && (CA<merge_epsilon)){MergeVertices.push_back(verts_ind3);}
        }
    }
}

bool MeshTools::decimateDegenerateTriangle(Mesh& mesh)
{
    std::vector<int> MergeFaceIndices;
    std::vector<INT2> MergeVerts;
    CheckDegenerateTriangle(mesh, MergeFaceIndices, MergeVerts);

    if (MergeFaceIndices.size() != 0){
        std::cout << "MergeFace Indices = " << MergeFaceIndices << "\n";
        auto& triangles = mesh.accesstriangles();
        auto& verts     = mesh.accessvertices();

        // declare new triangles and vertices
        std::vector<triangle> newT;
        std::vector<vertex>   newV;

        int nv = verts.size();
        int nf = triangles.size();

        std::vector<int> MapOldIndicesToNew = Algorithm::arange(0,nv,1);
        std::vector<int> Min(MergeVerts.size());
        std::vector<int> Max(MergeVerts.size());

        for (int i=0;i<MergeVerts.size();i++){
            Min[i] = MergeVerts[i][0];
            Max[i] = MergeVerts[i][1];
        }

        int index=0;
        for (int i=0;i<nv ;i++){
            int ind;
            if (Algorithm::contain(Max, i, ind)){
                MapOldIndicesToNew[i] = MapOldIndicesToNew[Min[ind]];
            }
            else{
                MapOldIndicesToNew[i] = index;
                newV.push_back(verts[i]);
                index ++;
            }
        }

        // obtain the new triangle indices 
        for (int i=0;i<triangles.size();i++){
            if (! Algorithm::contain(MergeFaceIndices,i)){
                triangle t;
                auto ind  = triangles[i].triangleindices_;
                t.triangleindices_ = {{MapOldIndicesToNew[ind[0]], MapOldIndicesToNew[ind[1]], MapOldIndicesToNew[ind[2]]}};
                if (Algorithm::is_unique(t.triangleindices_))
                {
                    newT.push_back(t);
                }
            }
        }

        triangles.clear();
        triangles = newT;
        verts.clear();
        verts = newV;

        return true;
    }

    return false;
}

void MeshTools::CorrectMesh(Mesh& mesh, std::vector<int>& FaceIndices)
{
    auto& triangles = mesh.accesstriangles();
    auto& vertices  = mesh.accessvertices();

    // declare the new triangles 
    std::vector<triangle> newT;
}

void MeshTools::CalculateBoundaryBarycenter(Mesh& mesh, std::vector<bool>& boundary_indicator, Real3& barycenter)
{
    // shift all boundary vertices wrt to the first vertex
    const auto& verts = mesh.getvertices();
    barycenter=verts[0].position_;

    for (int i=1;i<verts.size();i++){
        Real3 newpos = mesh.getShiftedVertexPosition(verts[i], verts[0]);
        barycenter = barycenter + newpos;
    }

    // average
    barycenter = barycenter * 1.0/verts.size();

    // shift into box
    Real3 shift = mesh.getShiftIntoBox(barycenter);

    barycenter = barycenter + shift;
}

bool MeshTools::IsIsolatedFace(Mesh& mesh, int faceIndex, const std::map<INT2, std::vector<int>>& mapEdgeToFace){
    const auto& faces = mesh.gettriangles();
    INT3 vInd  = faces[faceIndex].triangleindices_;

    int numBoundary=0;
    for (int i=0;i<3;i++){
        int next = (i+1) % 3;
        INT2 edge = MeshTools::makeEdge(vInd[i],vInd[next]);
        std::vector<int> faces;
        bool found = Algorithm::FindInMap(mapEdgeToFace, edge, faces);
        ASSERT((found), "The edge " << edge << " is not found.");

        if (faces.size() == 1){
            numBoundary += 1;
        }
    }

    if (numBoundary == 3){
        return true;
    }

    return false;
}

bool MeshTools::IsTeethlikeFace(Mesh& mesh, int faceIndex, const std::map<INT2, std::vector<int>>& mapEdgeToFace){
    const auto& faces = mesh.gettriangles();
    INT3 vInd  = faces[faceIndex].triangleindices_;

    int numBoundary=0;
    for (int i=0;i<3;i++){
        int next = (i+1) % 3;
        INT2 edge = MeshTools::makeEdge(vInd[i],vInd[next]);
        std::vector<int> faces;
        bool found = Algorithm::FindInMap(mapEdgeToFace, edge, faces);
        ASSERT((found), "The edge " << edge << " is not found.");

        if (faces.size() == 1){
            numBoundary += 1;
        }
    }

    if (numBoundary == 2){
        return true;
    }

    return false;
}

void MeshTools::ReconstructMeshAfterFaceCut(Mesh& mesh)
{
    std::vector<std::vector<int>> neighborIndex;
    MeshTools::CalculateVertexNeighbors(mesh, neighborIndex);

    const auto& v = mesh.getvertices();
    const auto& f = mesh.gettriangles();

    int index=0;
    std::vector<int> MapOldToNewIndex(v.size(),-1);
    for (int i=0;i<neighborIndex.size();i++){
        if (neighborIndex[i].size() != 0){
            MapOldToNewIndex[i] = index;            
            index++;
        }
    }

    std::vector<triangle> newFace;
    std::vector<vertex> newVerts;
    for (int i=0;i<f.size();i++){
        triangle t;
        for (int j=0;j<3;j++){
            int oldIndex = f[i].triangleindices_[j];
            int newIndex = MapOldToNewIndex[oldIndex];
            t.triangleindices_[j] = newIndex;        
        }
        newFace.push_back(t);
    }

    for (int i=0;i<v.size();i++){
        if (neighborIndex[i].size() != 0){
            newVerts.push_back(v[i]);
        }
    }

    auto& oldT = mesh.accesstriangles();
    auto& oldV = mesh.accessvertices();
    oldT.clear(); oldT.insert(oldT.end(), newFace.begin(), newFace.end());
    oldV.clear(); oldV.insert(oldV.end(), newVerts.begin(), newVerts.end());
}
