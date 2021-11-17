#pragma once
#include "Calculation.h"
#include "Bin.h"
#include "CalculationTools.h"
#include "Qtensor.h"
#include "LinAlgTools.h"

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
        virtual void finishCalculate() override {};

        // dots the vectors uij with respect to the largest eigenvector
        void printcos2thetaPerIter(std::ofstream& ofs);
        void printevPerIter(std::ofstream& ofs);
        void printp2PerIter(std::ofstream& ofs);
        void printQtensorPerIter(std::ofstream& ofs);
        void printcos2PerIter(std::ofstream& ofs);

        Real getBiaxiality() const {return biaxiality_;}
    
    private:
        Matrix Qtensor_;

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

        // P(costheta) wrt to each of the director

        int head_index_;
        int tail_index_;

        std::string residue_name_;

        int atomSize_;
        int size_;

        std::string pvName_;

        Real3 arr_ = {{0,0,1}};
};