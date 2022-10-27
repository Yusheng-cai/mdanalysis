/*
    Tables and conventions from
    http://paulbourke.net/geometry/polygonise/
*/

#include "marching_cubes.hpp"

int MarchingCubes::calculate_cube_index(GridCell &cell, Real isovalue)
{
    int cubeIndex = 0;
    for (int i = 0; i < 8; i++)
        if (cell.value[i] < isovalue) cubeIndex |= (1 << i);
    return cubeIndex;
}

std::vector<Point> MarchingCubes::get_intersection_coordinates(GridCell &cell, Real isovalue)
{
    std::vector<Point> intersections (12);

    int cubeIndex = calculate_cube_index(cell, isovalue);
    int intersectionsKey = edgeTable[cubeIndex];

    int idx = 0;
    while (intersectionsKey)
    {
        if (intersectionsKey&1)
        {
            int v1 = edgeToVertices[idx].first, v2 = edgeToVertices[idx].second;
            Point intersectionPoint = interpolate(cell.vertex[v1], cell.value[v1],
                                                    cell.vertex[v2], cell.value[v2], isovalue);
            intersections[idx] = intersectionPoint;
        }
        idx++;
        intersectionsKey >>= 1;
    }


    return intersections;
}

MarchingCubes::Real3 MarchingCubes::getPBCDistance(Real3& v1, Real3& v2)
{
    Real3 res;
    for (int i=0;i<3;i++)
    {
        Real dist = v2[i] - v1[i];

        if (dist > end_[i]/2.0) dist -= end_[i];
        else if (dist < -end_[i]/2.0) dist += end_[i];
        else dist=dist;
        res[i] = dist;
    }

    return res;
}


Point MarchingCubes::interpolate(Point& v1, Real val1, Point& v2, Real val2, Real isovalue)
{
    Point interpolated;
    Real mu = (isovalue - val1) / (val2 - val1);

    Real3 v1pos = {{v1.x, v1.y, v1.z}};
    Real3 v2pos = {{v2.x, v2.y, v2.z}};

    // v2 - v1
    Real3 dist;
    dist = getPBCDistance(v1pos, v2pos);
    interpolated.x = mu*dist[0] + v1.x;
    interpolated.y = mu*dist[1] + v1.y;
    interpolated.z = mu*dist[2] + v1.z;

    // let's put the interpolated positions into the box 
    Real3 interpolatedPos = {{interpolated.x, interpolated.y, interpolated.z}};
    for (int i=0;i<3;i++){
        if (interpolatedPos[i] < 0){
            interpolatedPos[i] += N_[i];
        }

        if (interpolatedPos[i] > N_[i]){
            interpolatedPos[i] -= N_[i];
        }

        if (std::abs(interpolatedPos[i] - N_[i]) < tol_[i]){
            interpolatedPos[i] = 0.0; 
        }
    }

    interpolated.x = interpolatedPos[0];
    interpolated.y = interpolatedPos[1];
    interpolated.z = interpolatedPos[2];

    return interpolated;
}

std::vector<std::vector<Point>> MarchingCubes::get_triangles(std::vector<Point>& intersections, int cubeIndex){
    std::vector<std::vector<Point>> triangles;
    for (int i = 0; triangleTable[cubeIndex][i] != -1; i += 3)
    {
        std::vector<Point> triangle (3);
        for (int j = 0; j < 3; j++)
            triangle[j] = intersections[triangleTable[cubeIndex][i + j]];
        triangles.push_back(triangle);
    }

    return triangles;
}

std::vector<std::vector<Point>> MarchingCubes::triangulate_cell(GridCell &cell, Real isovalue){
    int cubeIndex = calculate_cube_index(cell, isovalue);
    std::vector<Point> intersections = get_intersection_coordinates(cell, isovalue);
    std::vector<std::vector<Point>> triangles = get_triangles(intersections, cubeIndex);

    return triangles;
}


void MarchingCubes::VerticesForGridCell(INT3& index, std::vector<Point>& initialV, std::vector<Point>& neighborV){
    int intIndex = GridIndexToIndex(index);

    // these are the points from the self grid cell
    for (auto tri : triangles_[intIndex]){
        for (auto vert : tri){
            initialV.push_back(vert);
        }
    }

    // if we are performing with periodic boundary condition
    for (INT3 off : offsets_){
        // calculate the offset index
        INT3 updated = index + off;

        // check if updated is a valid point --> only change if non pbc
        bool valid = true;

        // * * * *(end)
        // for a system with 4 points, our grid cell only go to the 3rd one
        if (! pbc_){
            for (int i=0;i<3;i++){
                if (updated[i] >= (N_[i]-1)){
                    valid = false;
                }
            }
        }

        // if this is a valid system
        if (valid){
            int ind = GridIndexToIndex(updated);

            // ind = -1 happens when we are converting cell grid with no pbc but ended up outside of grid 
            for (auto& tri : triangles_[ind]){
                for (auto vert : tri){
                    neighborV.push_back(vert);
                }
            }
        }
    }

    return; 
}

void MarchingCubes::initializeGridSearch()
{
    for (int i=-1;i<=1;i++){
        for (int j=-1;j<=1;j++){
            for (int k=-1;k<=1;k++){
                INT3 ind = {{i,j,k}};
                offsets_.push_back(ind);
            }
        }
    }
}

int MarchingCubes::GridIndexToIndex(INT3& index){
    int x,y,z;
    if (index[0] < 0){x += N_[0];}else{x %= N_[0];}
    if (index[1] < 0){y += N_[1];}else{y %= N_[1];}
    if (index[2] < 0){z += N_[2];}else{z %= N_[2];}

    int idx = z * N_[0] * N_[1] + y * N_[0] + x;

    return idx;
}

MarchingCubes::INT3 MarchingCubes::IndexToGridIndex(int index){
    INT3 idx;
    idx[2] = std::round(index / (N_[0] * N_[1]));
    index = index - idx[2] * N_[0] * N_[1];
    idx[1] = std::round(index / N_[0]); 
    index = index - idx[1] * N_[0];
    idx[0] = index;

    return idx;
}

void MarchingCubes::CorrectPBCposition(Point& p)
{
    Real3 pos = {{p.x, p.y, p.z}};

    for (int i=0;i<3;i++){
        if (pos[i] > N_[i]) { pos[i] -= N_[i];}
        if (pos[i] < 0) { pos[i] += N_[i];}
        if (std::abs(pos[i] - N_[i]) < tol_[i]) { pos[i] = 0.0;}
    }

    p.x = pos[0];
    p.y = pos[1];
    p.z = pos[2];
}

void MarchingCubes::update(){
    triangles_.clear();
    MapFromCellGridIndexToIndex_.clear();
    MapFromIndexToCellGridIndex_.clear();
    offsets_.clear();
}

void MarchingCubes::triangulate_field(Lattice<Real>& field, Mesh& mesh, Real3 spacing, INT3 N, Real isovalue, bool pbc)
{
    // give an update --> clear the dictionaries etc.
    update();

    // field does the following thing * * * | --> 3 lattice points and 3 segments
    // whether or not we are performing periodic mesh
    pbc_ = pbc;

    // spacing of lattice points in the field 
    spacing_ = spacing;
    N_ = N;
    for (int i =0;i<3;i++){end_[i] = N_[i];}

    // Find where the end points are  --> that constructs the pbc box too
    Real spaceX = spacing_[0];
    Real spaceY = spacing_[1];
    Real spaceZ = spacing_[2];

    /// Whether we are performing pbc or not
    if (pbc) {inc_ = 0;}else {inc_=1;}
    
    // Perform the marching cubes 
    int Size = (N_[0] - inc_) * (N_[1] - inc_) * (N_[2] - inc_);
    triangles_.resize(Size);

    // max 
    INT3 max = N_ - inc_;
    std::vector<int> NumVertsPerCell(Size,0);

    #pragma omp parallel 
    {
        #pragma omp for collapse(3)
        for (int i = 0; i < max[0]; i++){
            for (int j = 0; j < max[1]; j++){
                for (int k = 0; k < max[2]; k++){
                    // define the x y z of this particular cell grid 
                    int x = i, y = j, z = k;
                    INT3 index = {{i,j,k}};
                    Real xplus, yplus, zplus;

                    // Correct for the positions 
                    xplus = (x+1) % N_[0];
                    yplus = (y+1) % N_[1];
                    zplus = (z+1) % N_[2];

                    // find the index 
                    int idx = GridIndexToIndex(index);
                    
                    // cell ordered according to convention in referenced website
                    GridCell cell = 
                    {
                        {
                            {x, y, z}, {xplus, y, z}, 
                            {xplus, y, zplus}, {x, y, zplus},
                            {x, yplus, z}, {xplus, yplus, z}, 
                            {xplus, yplus, zplus}, {x, yplus, zplus}
                        },

                        // indexing into field already takes care of pbc
                        {
                            field(i,j,k), field(i+1, j, k), 
                            field(i+1, j, k + 1), field(i, j, k+1), 
                            field(i, j+1, k), field(i+1, j+1, k), 
                            field(i+1, j+1, k+1), field(i, j+1, k+1)
                        }
                    };
                    
                    // obtain the triangles in cell grid
                    std::vector<std::vector<Point>> cellTriangles = triangulate_cell(cell, isovalue);

                    // keep track of the triangles 
                    triangles_[idx] = cellTriangles; 

                    // keep track of number of vertices per cell
                    NumVertsPerCell[idx] = cellTriangles.size() * 3;
                }
            }
        }
    }

    // sum the num verts per cell
    std::vector<int> partial_sum_num_verts(NumVertsPerCell.size(),0);
    std::partial_sum(NumVertsPerCell.begin(), NumVertsPerCell.end(), partial_sum_num_verts.begin(), std::plus<Real>());

    // index all vertices 
    std::vector<int> MapVertexIndexToCellGridIndex;
    #pragma omp parallel for
    for (int i=0;i<Size;i++){
        int offset;
        if (i == 0){offset = 0;}
        else{offset = partial_sum_num_verts[i-1];}

        int index=0;
        for (auto& tri : triangles_[i]){
            for (auto& v : tri){
                v.index = offset + index;
                index++;
            }
        }
    }

    // Now make a map of all the indices to some initial INT value
    std::vector<int> MapFromVertexIndexToNewIndex(partial_sum_num_verts[Size-1]+1, INITIAL_);

    // start cell Grid 
    // make a list of the offsets
    initializeGridSearch();

    // We go through this voxel by voxel
    std::vector<Real3> AllVertices;
    for (int i=0;i<max[0];i++){
        for (int j=0;j<max[1];j++){
            for (int k=0;k<max[2];k++){
                INT3 initialIndex={{i,j,k}};

                // stores neighbor vertices and self Vertices
                std::vector<Point> NeighborV;
                std::vector<Point> initialV;

                // find vertices for grid cell containing neighbor vertices and self vertices
                VerticesForGridCell(initialIndex, initialV, NeighborV); 

                //  
                int selfnum = initialV.size();
                int num = NeighborV.size();

                for (int m=0;m<selfnum;m++){
                    Real3 mpos = {{initialV[m].x, initialV[m].y, initialV[m].z}};
                    int GridIndex = initialV[m].index; 

                    std::vector<int> NeighborNumber;
                    std::vector<int> NeighborIndex;

                    for (int n=0;n<num;n++){
                        int indexN = NeighborV[n].index;
                        // avoid self calculation
                        if (indexN != GridIndex){
                            Real3 v2 = {{NeighborV[n].x, NeighborV[n].y, NeighborV[n].z}};

                            Real3 pbcdist = getPBCDistance(mpos, v2);

                            if (std::abs(pbcdist[0]) <= tol_[0] && std::abs(pbcdist[1]) <= tol_[1] && std::abs(pbcdist[2]) <= tol_[2]){
                                NeighborIndex.push_back(indexN);
                                int number = MapFromVertexIndexToNewIndex[indexN];
                                if (number != INITIAL_){
                                    NeighborNumber.push_back(number);
                                }
                            }
                        }
                    }

                    // First case scenario, all neighbors have not been added to the list yet 
                    if (NeighborNumber.size() == 0){
                        if (MapFromVertexIndexToNewIndex[GridIndex] == INITIAL_){
                            MapFromVertexIndexToNewIndex[GridIndex] = AllVertices.size();
                            CorrectPBCposition(initialV[m]);

                            Real3 mposCorrected = {{initialV[m].x * spaceX, initialV[m].y * spaceY , initialV[m].z * spaceZ}};
                            for (int index : NeighborIndex){
                                MapFromVertexIndexToNewIndex[index] = AllVertices.size();
                            }
                            AllVertices.push_back(mposCorrected);
                        }
                    }
                    else{
                        auto it = std::adjacent_find(NeighborNumber.begin(), NeighborNumber.end(), std::not_equal_to<>());
                        ASSERT((it == NeighborNumber.end()), "The elements in the list are not all identical.");

                        int IDX = NeighborNumber[0];

                        MapFromVertexIndexToNewIndex[GridIndex] = IDX;
                    }
                }
            }
        }
    }


    // vertices
    auto& verts = mesh.accessvertices();
    verts.clear();
    verts.resize(AllVertices.size());
    for (int i=0;i<verts.size();i++){verts[i].position_ = AllVertices[i];}

    auto& triMesh = mesh.accesstriangles();
    triMesh.clear();

    // fill the triangles 
    for (auto& block : triangles_){
        for (auto& tri : block){
            ASSERT((tri.size() == 3), "the triangle size must be 3.");
            std::vector<int> triIndex(3,0);
            for (int i=0;i<3;i++){
                triIndex[i] = MapFromVertexIndexToNewIndex[tri[i].index];
            }

            std::sort(triIndex.begin(), triIndex.end());
            auto it = std::unique(triIndex.begin(), triIndex.end());
            bool isUnique = (it == triIndex.end());

            if (isUnique){
                triangle t;
                for (int i=0;i<3;i++){
                    int index = MapFromVertexIndexToNewIndex[tri[i].index];
                    t.triangleindices_[i] = index;
                    t.vertices_[i] = verts[index];
                }

                triMesh.push_back(t);
            }
        }
    }

    mesh.CalcVertexNormals();
}