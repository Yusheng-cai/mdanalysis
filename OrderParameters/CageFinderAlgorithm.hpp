//
//                  ***GRADE_v1.00***
//  MyFunctions.hpp
//
//  Created by Farbod Mahmoudinobar and Cristiano L. Dias on 12/6/18.
//  Copyright © 2018 Farbod Mahmoudinobar and Cristiano L. Dias.
//
//  This file is part of GRADE.
//
//  GRADE is a free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  GRADE is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with GRADE.  If not, see <https://www.gnu.org/licenses/>.

#pragma once
#include <stdio.h>
#include <vector>
#include <string>
#include "tools/Algorithm.h"
#include "tools/Assert.h"
#include "tools/CommonOperations.h"

namespace CageFinderAlgorithm
{
    void find_shared_edges_ring5(int count_ring5, std::vector<std::vector<int>>& ring5, std::vector<std::vector<int>>& My_neigh_ring5, std::vector<unsigned long int>& N_ring5_neigh);

    void find_shared_edges_ring6(int count_ring6, std::vector<std::vector<int>>& ring6, std::vector<std::vector<int>>& My_neigh_ring6, std::vector<unsigned long int>& N_ring6_neigh);

    void find_shared_edges_ring6_ring5(std::vector<std::vector<int>>& ring5, std::vector<std::vector<int>>& ring6, std::vector<std::vector<int>>& My_neigh_ring6_ring5, std::vector<unsigned long int>& N_ring6_ring5_neigh );

    bool  compare(std::vector<int>& ringA, std::vector<int>& ringB, int ringA_N, int ringB_N );

    bool  compare_adjacant(std::vector<int>& ringA, std::vector<int>& ringB, int ringA_N, int ringB_N , std::vector<int>& base);

    int remove_duplicates_map(std::vector<std::vector<int>>& cup512);

    void cup_512_Finder(std::vector<std::vector<int>>& ring5,int count_ring5, std::vector<unsigned long int>& N_ring5_neigh, std::vector<std::vector<int>>& My_neigh_ring5, std::vector<std::vector<int>>& cup512 );

    void cup_62512_Finder(std::vector<std::vector<int>>& ring6, int count_ring6, std::vector<std::vector<int>>& My_neigh_ring6_ring5, std::vector<unsigned long int> N_ring6_ring5_neigh, std::vector<std::vector<int>>& ring5, std::vector<unsigned long int> N_ring5_neigh, std::vector<std::vector<int>>& My_neigh_ring5, std::vector<std::vector<int>>& cup62512);

    int cage_Finder(std::vector<std::vector<int>> cups, unsigned long count_cups, std::vector<std::vector<int>>& neighour_rings, std::vector<std::vector<int>>& cage, std::vector<std::vector<int>>& cage_rings, std::string time);

    int cage_Finder_64512(std::vector<std::vector<int>> cup62512, int count_62512_cups, std::vector<std::vector<int>>& cage_64512_rings);
}