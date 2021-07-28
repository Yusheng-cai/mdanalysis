#pragma once

#include "tools/CommonTypes.h"
#include "Assert.h"
#include "tools/GenericFactory.h"

#include <array>
#include <vector>

class Frame
{
    public:
        using Real  = float;
        using Real3 = std::array<Real,3>; 
        using VectorReal3 = std::vector<Real3>;
        using VectorInt = std::vector<int>;
        using Matrix = std::array<Real3,3>;

        Frame(){}; 
        ~Frame(){};
        VectorReal3& accessPositions(){return positions_;};
        VectorReal3& accessVelocities(){return velocities_;}
        VectorReal3& accessForces(){return forces_;}
        VectorInt& accessGlobalIndex(){return global_index_;};
        Matrix& accessBoxMatrix(){return box_;}

        const VectorReal3& getPositions() const{return positions_;}
        const VectorInt& getGlobalIndex() const{return global_index_;}
        const VectorReal3& getVelocities()const {return velocities_;}
        const VectorReal3& getForces()const{return forces_;}
        const Matrix& getBoxMatrix()const{return box_;}

        Real getTime(){return time_;}
        void setTime(Real time){time_ = time;}
        int getStep(){return step_;}
        void setStep(int step){step_ = step;}
        void setBox(const Matrix& box)
        {
            for(int i=0;i<3;i++)
            {
                for(int j=0;j<3;j++)
                {
                    box_[i][j] = box[i][j];
                }
            }
        }

        void setNumAtoms(int num_atoms);

    private:
        VectorReal3 positions_, velocities_, forces_;
        VectorInt global_index_;
        Real time_=0;
        Matrix box_;

        int step_ = 0;
        int num_atoms_;
};
