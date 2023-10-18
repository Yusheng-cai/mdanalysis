#pragma once

#include "Calculation.h"
#include "CellGrid.h"
#include "Eigen/Core"
#include <Eigen/SVD>

#include <string>
#include <vector>
#include <memory>

class CageFinder : public Calculation
{
    public:
        using cellptr = std::unique_ptr<CellGrid>;
        CageFinder(const CalculationInput& input);
        virtual void calculate() override;
        virtual void update();
        virtual void finishCalculate() {};
        void RemoveDuplicateFaces(std::vector<std::vector<int>>& Faces);
        void RemoveDuplicate(std::vector<std::vector<int>>& Faces);
        void CheckCoplanar(std::vector<std::vector<int>>& Faces);
        void CheckConvex(std::vector<std::vector<int>>& Faces);

        void FindOccupancy512();
        void FindOccupancy62512();

        void PrintCage512(std::ofstream& ofs);
        void PrintCage62512(std::ofstream& ofs);
        void PrintPerIterOccupation(std::ofstream& ofs);
        void PrintNonOccupied512(std::ofstream& ofs);
        void PrintNonOccupied62512(std::ofstream& ofs);
    
    private:
        std::string atomgroup_, solute_group_;
        std::vector<std::vector<int>> neighbor_list_;

        std::vector<std::vector<int>> solvent_neighbor_indices_;
        Real cut_off_=0.35;
        Real cut_off_sq_;
        Real plane_cutoff_=0.1;
        Real solute_512_cutoff_=0.6, solute_512_cutoffsq_;
        Real solute_62512_cutoff_=0.8, solute_62512_cutoffsq_;
        Real internal_angle_cutoff_=10;

        cellptr cell_;

        std::vector<std::vector<int>> ring5_;
        std::vector<std::vector<int>> ring6_;
        std::vector<std::vector<int>> Cage_512_rings_;
        std::vector<std::vector<int>> Cage_62512_rings_;
        std::vector<int> Non_occupied_cage62512_;

        int occupied_512_;
        int occupied_62512_;
};