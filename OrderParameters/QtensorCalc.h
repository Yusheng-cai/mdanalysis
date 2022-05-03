#pragma once
#include "Calculation.h"
#include "Bin.h"
#include "LinAlgTools.h"
#include "SimulationState.h"

#include <string>
#include <numeric>
#include <algorithm>
#include <vector>
#include <array>
#include <memory>
#include <map>
#include <iostream>
#include <iomanip>

// calculate distribution of cosine theta of the angles in a probe volume

class QtensorCalc : public Calculation
{
    public:
        using Matrix = CommonTypes::Matrix;
        using Real3  = CommonTypes::Real3;
        using Real   = CommonTypes::Real;

        QtensorCalc(const CalculationInput& input);

        virtual void calculate() override;
        virtual void finishCalculate() override;

        // dots the vectors uij with respect to the largest eigenvector
        void printcos2thetaPerIter(std::ofstream& ofs);
        void printevPerIter(std::ofstream& ofs);
        void printp2PerIter(std::ofstream& ofs);
        void printQtensorPerIter(std::ofstream& ofs);
        void printcos2PerIter(std::ofstream& ofs);
        void printcosPerIter(std::ofstream& ofs);
        void printaverageQ(std::string name);

        Real getBiaxiality() const {return biaxiality_;}
        Real getP2() const {return p2_;}
        Real getQxx() const {return Qtensor_[0][0];}
        Real getQxy() const {return Qtensor_[0][1];}
        Real getQxz() const {return Qtensor_[0][2];}
        Real getQyy() const {return Qtensor_[1][1];}
        Real getQyz() const {return Qtensor_[1][2];}
        Real getv1x() const {return eigenvector_[0][0];}
        Real getv1y() const {return eigenvector_[1][0];}
        Real getv1z() const {return eigenvector_[2][0];}
    
    private:
        Matrix Qtensor_;
        // The average Qtensor over time
        Matrix QtensorTot_;

        // The director of each of the molecules
        std::vector<Real3> uij_;
        std::vector<Real> norm_;

        std::vector<Real3> COM_;

        // The eigenvector directions 
        Matrix eigenvector_;
        Real3  eigenval_;
        Real3 v0_;
        Real3 v1_;
        Real3 v2_;
        Real biaxiality_;
        Real p2_;

        // P(costheta) wrt to each of the director

        int head_index_;
        int tail_index_;

        std::string residue_name_;

        int atomSize_;
        int size_;

        Real3 arr_ = {{0,0,1}};
};