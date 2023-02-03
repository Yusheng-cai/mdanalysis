#include "Calculation.h"
#include "CellGrid.h"
#include "tools/CommonOperations.h"
#include "LinAlgTools.h"
#include "Lattice.hpp"
#include "parallel/OpenMP_buffer.h"

class QtensorVoxel : public Calculation
{
    public:
        QtensorVoxel(const CalculationInput& input);

        virtual void calculate() override;

        virtual void update() override;

        virtual void finishCalculate() override;

        INT3 findVoxel(Real3& position);

        void reset_buffers();

    private:
        INT3 voxel_size_;
        Real3 voxel_separation_;

        std::string residue_group_;
        int head_index_=1, tail_index_=2;

        std::vector<Real3> uij_;

        Matrix Global_Qtensor_;

        Lattice<Matrix> Qtensor_Voxels_;
        Lattice<Real> NumAtoms_Voxels_;
        Lattice<Real> Order_Voxels_;
        Lattice<Real3> Director_Voxels_;

        // qtensor, numatoms buffer
        OpenMP::OpenMP_buffer<Lattice<Matrix>> Qtensor_voxels_buffer_;
        OpenMP::OpenMP_buffer<Lattice<Real>> NumAtoms_voxels_buffer_;

        Matrix zero_matrix_;
};