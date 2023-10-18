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


#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include "CageFinderAlgorithm.hpp"

// New Function--------------------------------------------------------------------------------------

//using namespace std::chrono;

// New Function--------------------------------------------------------------------------------------
void CageFinderAlgorithm::find_shared_edges_ring5(int count_ring5, std::vector<std::vector<int>>& ring5, std::vector<std::vector<int>>& My_neigh_ring5, std::vector<unsigned long int>& N_ring5_neigh){
    
    // Finding shared edges for 5-membered rings
   std::vector<int> temp_vecn;
    temp_vecn.clear();
    int counter=0;
    int neigh_ring_counter=0;
    //std::vector<int> N_ring6_ring5_neigh;
    std::vector<std::vector<int>> ring6_ring5_neigh;
    std::vector<int> temp_vec4 = {0,0,0,0};

    // Map edge to ring5 
    std::map<std::array<int,2>, std::vector<int>> MapEdgeToRing5;

    // find all the edges in ring5
    for (int i=0;i<ring5.size();i++){
        for (int j=0;j<ring5[i].size()-1;j++){
            std::array<int,2> edge = {ring5[i][j], ring5[i][j+1]};
            std::sort(edge.begin(), edge.end());
            std::vector<int> temp;

            if (Algorithm::IsInMap(MapEdgeToRing5, edge, temp)){
                Algorithm::InsertInVectorMap(MapEdgeToRing5, edge, i);
            }
            else{
                std::vector<int> one = {i};
                Algorithm::InsertInMap(MapEdgeToRing5, edge, one);
            }
        }
    }

    for (int i=0;i<ring5.size();i++){
        std::vector<int> temp;
        for (int j=0;j<ring5[i].size()-1;j++){
            std::array<int,2> edge = {ring5[i][j], ring5[i][j+1]};
            std::sort(edge.begin(), edge.end());

            std::vector<int> Ring5_indices;
            ASSERT(Algorithm::IsInMap(MapEdgeToRing5, edge, Ring5_indices), "Ring 5 is wrong.");
            temp.insert(temp.end(), Ring5_indices.begin(), Ring5_indices.end());
        }
        std::sort(temp.begin(), temp.end());
        temp.erase(unique(temp.begin(), temp.end()), temp.end());
        My_neigh_ring5.push_back(temp);
    }
    
    N_ring5_neigh.clear();
    for (int i = 0 ; i < ring5.size() ; i++){
        N_ring5_neigh.push_back(My_neigh_ring5[i].size()) ;
    }
    
    
    
}

// New Function--------------------------------------------------------------------------------------

void CageFinderAlgorithm::find_shared_edges_ring6(int count_ring6, std::vector<std::vector<int>>& ring6, std::vector<std::vector<int>>& My_neigh_ring6, std::vector<unsigned long int>& N_ring6_neigh){
    
    std::vector<int> temp_vecn;
    temp_vecn.clear();
    int counter=0;
    int neigh_ring_counter=0;
    //std::vector<int> N_ring6_neigh;
    std::vector<std::vector<int>> ring6_neigh;
    std::vector<int> temp_vec4 = {0,0,0,0};
    
    
    for (int i = 0 ; i < count_ring6 ; i++)
    {
        for (int j = 0 ; j < ring6[i].size() ; j++)
        {
            for (int k = 0 ; k < count_ring6 ; k++)
            {
                for (int l=0 ; l < ring6[k].size() ; l++)
                {
                    if (i != k && j+1 < ring6[i].size() && ( ( ring6[i][j] == ring6[k][l] && ring6[i][j+1] == ring6[k][l+1]) || (ring6[i][j] == ring6[k][l+1] && ring6[i][j+1] == ring6[k][l]) ))
                    {
                        counter++;
                        neigh_ring_counter++;
                        temp_vec4[0] = i+1;
                        temp_vec4[1] = k+1;
                        temp_vec4[2] = ring6[i][j];
                        temp_vec4[3] = ring6[i][j+1];
                        
                        ring6_neigh.push_back(temp_vec4);
                        
                        temp_vecn.push_back(k);
                    }
                }
            }
        }
        N_ring6_neigh.push_back(neigh_ring_counter);
        neigh_ring_counter = 0;
        
        sort(temp_vecn.begin(), temp_vecn.end());
        temp_vecn.erase(unique(temp_vecn.begin(), temp_vecn.end()), temp_vecn.end());
        My_neigh_ring6.push_back(temp_vecn);
        temp_vecn.clear();
        
    }
    
    N_ring6_neigh.clear();
    for (int i = 0 ; i < count_ring6 ; i++){
        N_ring6_neigh.push_back(My_neigh_ring6[i].size()) ;
    }
    
}

// New Function--------------------------------------------------------------------------------------

void CageFinderAlgorithm::find_shared_edges_ring6_ring5(std::vector<std::vector<int>>& ring5, std::vector<std::vector<int>>& ring6, std::vector<std::vector<int>>& My_neigh_ring6_ring5, std::vector<unsigned long int>& N_ring6_ring5_neigh ){
    
    std::vector<int> temp_vecn;
    temp_vecn.clear();
    int counter=0;
    int neigh_ring_counter=0;
    //std::vector<int> N_ring6_ring5_neigh;
    std::vector<std::vector<int>> ring6_ring5_neigh;
    std::vector<int> temp_vec4 = {0,0,0,0};

    // Map edge to ring5 
    std::map<std::array<int,2>, std::vector<int>> MapEdgeToRing5;

    // find all the edges in ring5
    for (int i=0;i<ring5.size();i++){
        for (int j=0;j<ring5[i].size()-1;j++){
            std::array<int,2> edge = {ring5[i][j], ring5[i][j+1]};
            std::sort(edge.begin(), edge.end());
            std::vector<int> temp;

            if (Algorithm::IsInMap(MapEdgeToRing5, edge, temp)){
                Algorithm::InsertInVectorMap(MapEdgeToRing5, edge, i);
            }
            else{
                std::vector<int> one = {i};
                Algorithm::InsertInMap(MapEdgeToRing5, edge, one);
            }
        }
    }

    for (int i=0;i<ring6.size();i++){
        std::vector<int> temp;
        for (int j=0;j<ring6[i].size()-1;j++){
            std::array<int,2> edge = {ring6[i][j], ring6[i][j+1]};
            std::sort(edge.begin(), edge.end());

            std::vector<int> Ring5_indices;
            if (Algorithm::IsInMap(MapEdgeToRing5,edge, Ring5_indices)){
                temp.insert(temp.end(), Ring5_indices.begin(), Ring5_indices.end());
            }
        }
        std::sort(temp.begin(), temp.end());
        temp.erase(unique(temp.begin(), temp.end()), temp.end());
        My_neigh_ring6_ring5.push_back(temp);
    }
    
    N_ring6_ring5_neigh.clear();
    for (int i = 0 ; i < ring6.size() ; i++){
        N_ring6_ring5_neigh.push_back(My_neigh_ring6_ring5[i].size()) ;
    }
    
}


// New Function--------------------------------------------------------------------------------------

// New Function--------------------------------------------------------------------------------------

//This function compares two rings for shared sides. Each side is made up of two consecutive atoms.
//Every two consecutive atoms from ringA are compared with every two consecutive atoms of ringB.
//For example, if ringA={1,3,5,7,9} and ringB={2,4,6,8,10}, {1,3} will be compared with both {2,4} and {4,2} since both cases should be considered.

bool CageFinderAlgorithm::compare(std::vector<int>& ringA, std::vector<int>& ringB, int ringA_N, int ringB_N )
{
    int neigh_ring_counter=0;
    std::map<std::array<int,2>,bool> m;

    for (int i=0;i<ringA_N;i++){
        std::array<int,2> edge = {ringA[i], ringA[i+1]};
        std::sort(edge.begin(), edge.end());
        bool a;

        if (! Algorithm::IsInMap(m, edge, a)){
            Algorithm::InsertInMap(m, edge, true);
        }
    }

    for (int i=0;i<ringB_N;i++){
        std::array<int,2> edge = {ringB[i], ringB[i+1]};
        std::sort(edge.begin(), edge.end());
        bool a;
        if (Algorithm::IsInMap(m, edge,a)){
            neigh_ring_counter++;
        }
    }
    
    // for (int i = 0 ; i < ringA_N  ; i++)
    // {
    //     for (int j = 0 ; j < ringB_N  ; j++)
    //     {
            
    //         if  ( (ringA[i] == ringB[j] && ringA[i+1] == ringB[j+1]) || (ringA[i] == ringB[j+1] && ringA[i+1] == ringB[j]) )
    //         {
    //             neigh_ring_counter++;
    //             break;
    //         }
    //     }
    // }
    if (neigh_ring_counter == 0)
    {
        return false;
    }
    else return true;
}       //close compare function

// New Function---------------------------------------------------------------------------------------

// This functions compares 3 rings. If two 5-rings have a common edge/side, then it compares the common side with the sides of base ring(6-ring).
// If the common side b/w two 5-rings are shared with the base ring, these two 5-rings are not neighbors and function returns "false", otherwise
// it returns "true".
bool CageFinderAlgorithm::compare_adjacant(std::vector<int>& ringA, std::vector<int>& ringB, int ringA_N, int ringB_N , std::vector<int>& base)
{
    int neigh_ring_counter=0;
    std::map<std::array<int,2>,bool> m;
    std::map<std::array<int,2>,bool> shared_edges;

    for (int i=0;i<ringA_N;i++){
        std::array<int,2> edge = {ringA[i], ringA[i+1]};
        std::sort(edge.begin(), edge.end());
        bool a;

        if (! Algorithm::IsInMap(m, edge, a)){
            Algorithm::InsertInMap(m, edge, true);
        }
    }

    for (int i=0;i<ringB_N;i++){
        std::array<int,2> edge = {ringB[i], ringB[i+1]};
        std::sort(edge.begin(), edge.end());
        bool a;
        if (Algorithm::IsInMap(m, edge,a)){
            neigh_ring_counter++;
            Algorithm::InsertInMap(shared_edges, edge, true);
        }
    }

    for (int i=0;i<base.size()-1;i++){
        std::array<int,2> edge = {base[i], base[i+1]};
        std::sort(edge.begin(), edge.end());
        bool a;

        if (Algorithm::IsInMap(shared_edges, edge,a)){
            return false;
        }
    }

    if (neigh_ring_counter == 0){
        return false;
    }
    else{
        return true;
    }
    
    // for (int i = 0 ; i < ringA_N  ; i++)
    // {
    //     for (int j = 0 ; j < ringB_N ; j++)
    //     {
    //         if  ( (ringA[i] == ringB[j] && ringA[i+1] == ringB[j+1]) || (ringA[i] == ringB[j+1] && ringA[i+1] == ringB[j] ) )       //compare sides of the two 5-rings. If they have a common side, compare the common side with base ring.
    //         {
    //             for (int k6 = 0 ; k6 < 6; k6++)     // Loop to go over all the sides of base ring.
    //             {
    //                 if( (ringA[i] == base[k6] && ringA[i+1] == base[k6+1])  || (ringA[i] == base[k6+1] && ringA[i+1] == base[k6]) )
    //                 {
    //                     return false;
    //                 }
    //             }
    //             neigh_ring_counter++;
    //             break;
    //         }
    //     }
    // }

    // if (neigh_ring_counter == 0)
    // {
    //     return false;
    // }
    // else return true;
}       //close compare function

// New Function---------------------------------------------------------------------------------------


/* This Functions reads in a std::vector<std::vector<int>> which holds all the cages/cups/rings and removes the duplicates from them. Then it
 writes the new list without duplicates inside the original input file. This function returns an "int" which is the count number of
 the cups/rings without duplicates. Pay attention that the after this function is executed, the original input variable still holds the same number of rows/columns as before, but the first "Count" number of rows hold the results without any duplicate. So it is crucial to use the first "Count" rows in future functions. 
 */

int CageFinderAlgorithm::remove_duplicates_map(std::vector<std::vector<int>>& cups)
{
    unsigned long int cup_size;
    if(cups.size() != 0){cup_size = cups[0].size();}         //Number of rings in each cup, including base ring.
    else cup_size = 0;
    
    unsigned long int N_cups = cups.size();              //Number of cups before removing the duplicates.
    int count=0;                                                //Number of cups.
    std::vector<int> current;                                        //Cup which is being read and compared at each loop.
    std::map<std::vector<int>,int> mymap;                                 //A map which counts the frequencty of each cup( [cup,frequency] )
    int prev = 0 ;
    
    for (int i = 0 ; i < N_cups ; i++)      //Loop over all the cups.
    {
        for(int j = 0 ; j < cup_size ; j++){current.push_back(cups[i][j]);}      //Loop to put the read cup into "current".
        sort(current.begin(), current.end());
        int duplicate_count= 0;
        for(int j = 0 ; j < current.size()-1 ; j++ )
        {
            if(current[j] == current[j+1]) {duplicate_count++;}
        }
        
        if(duplicate_count == 0)
        {
            if(!mymap[current])                  //For each read cup, check if it exists in the map. If it does NOT, write it in the map.
                
            {
                
                mymap[current] = 1;
                count++;
                for (int k = 0 ; k < cup_size ; k++)
                {
                    cups[prev][k] = cups[i][k];   //If the cup does not exist in map, write it to first line of original cup variable.
                }
                prev++;
            }
            else {mymap[current]++;}        //If the cup already exists in the map, jump the writing step and add one to frequencty of that cup.
        }
        else duplicate_count = 0;
            
            current.clear();
    }
    
    return count;
}       // close remove_dupclicates_maps


/*
 This Function finds cups (half of cage) of 5(12) cages. Sometimes cups and half_cups are used interchangeably("half" reminds me it is not a full cage and distinguishes them).
 It gets ring5 and all related info about ring5 and gives back cup512.
 */
void CageFinderAlgorithm::cup_512_Finder(std::vector<std::vector<int>>& ring5,int count_ring5, std::vector<unsigned long int>& N_ring5_neigh, std::vector<std::vector<int>>& My_neigh_ring5, std::vector<std::vector<int>>& cup512 ){
    
    // finding half of 5(12) cages
    
    std::vector<int> temp_vecn, temp_vec5={0,0,0,0,0,0,};
    temp_vecn.clear();
    int first_ring, second_ring, third_ring, fourth_ring, fifth_ring, sixth_ring;       //These refer to first, second, third, ... ring of cup. First ring is the base of the cup.
    std::vector<int> temp_ring1, temp_ring2, temp_ring3, temp_ring4, temp_ring5, temp_ring6, temp_ring7;
    std::vector<int> repeated_atoms;
    
    for (int index1 = 0 ; index1 < count_ring5 ; index1++)
    {
        first_ring = index1;
        if (N_ring5_neigh[first_ring] > 4)
        {            
            for (int index2 = 0 ; index2 < N_ring5_neigh[first_ring] ; index2++)
            {
                second_ring = My_neigh_ring5[first_ring][index2];
                if (N_ring5_neigh[second_ring] > 2)
                {
                    if(second_ring != first_ring)
                    {
                        for (int index3 = 0 ; index3 < N_ring5_neigh[second_ring] ; index3++)
                        {
                            if (My_neigh_ring5[second_ring][index3] != first_ring &&
                                My_neigh_ring5[second_ring][index3] != second_ring
                                )
                            {
                                if(compare(ring5[first_ring], ring5[My_neigh_ring5[second_ring][index3]], 5, 5))
                                {
                                    third_ring = My_neigh_ring5[second_ring][index3];
                                    
                                    if (N_ring5_neigh[third_ring] > 2)
                                    {
                                        for (int index4 = 0 ; index4 < N_ring5_neigh[third_ring]; index4++)
                                        {
                                            if (My_neigh_ring5[third_ring][index4] != first_ring &&
                                                My_neigh_ring5[third_ring][index4] != second_ring &&
                                                My_neigh_ring5[third_ring][index4] != third_ring &&
                                                compare(ring5[first_ring], ring5[My_neigh_ring5[third_ring][index4]], 5, 5) &&
                                                compare_adjacant(ring5[third_ring], ring5[My_neigh_ring5[third_ring][index4]], 5, 5, ring5[first_ring])
                                                )
                                            {
                                                fourth_ring = My_neigh_ring5[third_ring][index4];
                                                if (N_ring5_neigh[fourth_ring] > 2)
                                                {
                                                    for (int index5 = 0 ; index5 < N_ring5_neigh[fourth_ring]; index5++)
                                                    {
                                                        if (My_neigh_ring5[fourth_ring][index5] != first_ring &&
                                                            My_neigh_ring5[fourth_ring][index5] != second_ring &&
                                                            My_neigh_ring5[fourth_ring][index5] != third_ring  &&
                                                            My_neigh_ring5[fourth_ring][index5] != fourth_ring &&
                                                            compare(ring5[first_ring], ring5[My_neigh_ring5[fourth_ring][index5]], 5, 5) &&
                                                            compare_adjacant(ring5[fourth_ring], ring5[My_neigh_ring5[fourth_ring][index5]], 5, 5, ring5[first_ring])
                                                            )
                                                        {
                                                            fifth_ring = My_neigh_ring5[fourth_ring][index5];
                                                            
                                                            if (N_ring5_neigh[fifth_ring] != 0 && N_ring5_neigh[fifth_ring] > 2)
                                                            {
                                                                for (int index6 = 0 ; index6 < N_ring5_neigh[fifth_ring]; index6++)
                                                                {
                                                                    if (My_neigh_ring5[fifth_ring][index6] != first_ring  &&
                                                                        My_neigh_ring5[fifth_ring][index6] != second_ring &&
                                                                        My_neigh_ring5[fifth_ring][index6] != third_ring  &&
                                                                        My_neigh_ring5[fifth_ring][index6] != fourth_ring &&
                                                                        My_neigh_ring5[fifth_ring][index6] != fifth_ring  &&
                                                                        compare(ring5[second_ring], ring5[My_neigh_ring5[fifth_ring][index6]], 5, 5) &&
                                                                        compare_adjacant(ring5[fifth_ring], ring5[My_neigh_ring5[fifth_ring][index6]], 5, 5, ring5[first_ring])
                                                                        )
                                                                    {
                                                                        sixth_ring = My_neigh_ring5[fifth_ring][index6];
                                                                        
                                                                        if (compare(ring5[first_ring], ring5[sixth_ring], 5, 5))
                                                                        {
                                                                            
                                                                            temp_ring1.clear();
                                                                            temp_ring2.clear();
                                                                            temp_ring3.clear();
                                                                            temp_ring4.clear();
                                                                            temp_ring5.clear();
                                                                            temp_ring6.clear();
                                                                            
                                                                            
                                                                            for (int i_1 = 0; i_1 < 5; i_1++)
                                                                            {
                                                                                temp_ring1.push_back(ring5[first_ring ][i_1]);
                                                                                temp_ring2.push_back(ring5[second_ring][i_1]);
                                                                                temp_ring3.push_back(ring5[third_ring ][i_1]);
                                                                                temp_ring4.push_back(ring5[fourth_ring][i_1]);
                                                                                temp_ring5.push_back(ring5[fifth_ring ][i_1]);
                                                                                temp_ring6.push_back(ring5[sixth_ring ][i_1]);
                                                                                
                                                                            }
                                                                            
                                                                            repeated_atoms.clear();
                                                                            set_intersection(temp_ring1.begin(),temp_ring1.end(),temp_ring2.begin(),temp_ring2.end(),back_inserter(repeated_atoms));    //check ring1 and ring2 for overlaps
                                                                            
                                                                            if (repeated_atoms.size() <= 2)
                                                                            {
                                                                                repeated_atoms.clear();
                                                                                set_intersection(temp_ring2.begin(),temp_ring2.end(),temp_ring3.begin(),temp_ring3.end(),back_inserter(repeated_atoms));
                                                                                if (repeated_atoms.size() <= 2)
                                                                                {
                                                                                    repeated_atoms.clear();
                                                                                    set_intersection(temp_ring3.begin(),temp_ring3.end(),temp_ring4.begin(),temp_ring4.end(),back_inserter(repeated_atoms));
                                                                                    if (repeated_atoms.size() <= 2)
                                                                                    {
                                                                                        repeated_atoms.clear();
                                                                                        set_intersection(temp_ring4.begin(),temp_ring4.end(),temp_ring5.begin(),temp_ring5.end(),back_inserter(repeated_atoms));
                                                                                        if (repeated_atoms.size() <= 2)
                                                                                        {
                                                                                            repeated_atoms.clear();
                                                                                            set_intersection(temp_ring5.begin(),temp_ring5.end(),temp_ring6.begin(),temp_ring6.end(),back_inserter(repeated_atoms));
                                                                                            if (repeated_atoms.size() <= 2)
                                                                                            {
                                                                                                repeated_atoms.clear();
                                                                                                set_intersection(temp_ring6.begin(),temp_ring6.end(),temp_ring1.begin(),temp_ring1.end(),back_inserter(repeated_atoms));
                                                                                                
                                                                                                temp_vec5[0] = first_ring ;
                                                                                                temp_vec5[1] = second_ring;
                                                                                                temp_vec5[2] = third_ring ;
                                                                                                temp_vec5[3] = fourth_ring;
                                                                                                temp_vec5[4] = fifth_ring ;
                                                                                                temp_vec5[5] = sixth_ring ;
                                                                                                
                                                                                                cup512.push_back(temp_vec5);
                                                                                                
                                                                                                temp_ring1.clear();
                                                                                                temp_ring2.clear();
                                                                                                temp_ring3.clear();
                                                                                                temp_ring4.clear();
                                                                                                temp_ring5.clear();
                                                                                                temp_ring6.clear();
                                                                                                
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }       //end of finding half of 5(12) cages
    
    
}

// New Function--------------------------------------------------------------------------------------

// New Function--------------------------------------------------------------------------------------

void CageFinderAlgorithm::cup_62512_Finder(std::vector<std::vector<int>>& ring6, int count_ring6, std::vector<std::vector<int>>& My_neigh_ring6_ring5, std::vector<unsigned long int> N_ring6_ring5_neigh, std::vector<std::vector<int>>& ring5, std::vector<unsigned long int> N_ring5_neigh, std::vector<std::vector<int>>& My_neigh_ring5, std::vector<std::vector<int>>& cup62512){
    
    
    // finding half of 6(2)5(12) cages
    
    
    int first_ring, second_ring, third_ring, fourth_ring, fifth_ring, sixth_ring, seventh_ring;
    std::vector<int> temp_ring1, temp_ring2, temp_ring3, temp_ring4, temp_ring5, temp_ring6, temp_ring7;
    std::vector<int> repeated_atoms, temp_vec6={0,0,0,0,0,0,0};
    
    for (int index1 = 0 ; index1 < count_ring6 ; index1++)
    {
        first_ring = index1;
        if (N_ring6_ring5_neigh[first_ring] > 5)
        {
            for (int index2 = 0 ; index2 < N_ring6_ring5_neigh[first_ring] ; index2++)
            {
                second_ring = My_neigh_ring6_ring5[first_ring][index2];
                if (N_ring5_neigh[second_ring] > 3)
                {
                    if(second_ring != first_ring)
                    {
                        for (int index3 = 0 ; index3 < N_ring5_neigh[second_ring] ; index3++)
                        {
                            if (My_neigh_ring5[second_ring][index3] != first_ring &&
                                My_neigh_ring5[second_ring][index3] != second_ring
                                )
                            {
                                if(compare(ring6[first_ring], ring5[My_neigh_ring5[second_ring][index3]], 6, 5) &&
                                   compare_adjacant(ring5[second_ring], ring5[My_neigh_ring5[second_ring][index3]], 5, 5, ring6[first_ring]) )
                                {
                                    third_ring = My_neigh_ring5[second_ring][index3];
                                    if (N_ring5_neigh[third_ring] > 3)
                                    {
                                        for (int index4 = 0 ; index4 < N_ring5_neigh[third_ring]; index4++)
                                        {
                                            if (My_neigh_ring5[third_ring][index4] != first_ring &&
                                                My_neigh_ring5[third_ring][index4] != second_ring &&
                                                My_neigh_ring5[third_ring][index4] != third_ring &&
                                                compare(ring6[first_ring], ring5[My_neigh_ring5[third_ring][index4]], 6, 5) &&
                                                compare_adjacant(ring5[third_ring], ring5[My_neigh_ring5[third_ring][index4]], 5, 5, ring6[first_ring])
                                                )
                                            {
                                                fourth_ring = My_neigh_ring5[third_ring][index4];
                                                if (N_ring5_neigh[fourth_ring] > 3)
                                                {
                                                    for (int index5 = 0 ; index5 < N_ring5_neigh[fourth_ring]; index5++)
                                                    {
                                                        if (My_neigh_ring5[fourth_ring][index5] != first_ring &&
                                                            My_neigh_ring5[fourth_ring][index5] != second_ring &&
                                                            My_neigh_ring5[fourth_ring][index5] != third_ring  &&
                                                            My_neigh_ring5[fourth_ring][index5] != fourth_ring &&
                                                            compare(ring6[first_ring], ring5[My_neigh_ring5[fourth_ring][index5]], 6, 5) &&
                                                            compare_adjacant(ring5[fourth_ring], ring5[My_neigh_ring5[fourth_ring][index5]], 5, 5, ring6[first_ring] )
                                                            )
                                                        {
                                                            fifth_ring = My_neigh_ring5[fourth_ring][index5];
                                                            if (N_ring5_neigh[fifth_ring] != 0 && N_ring5_neigh[fifth_ring] > 3)
                                                            {
                                                                for (int index6 = 0 ; index6 < N_ring5_neigh[fifth_ring]; index6++)
                                                                {
                                                                    if (My_neigh_ring5[fifth_ring][index6] != first_ring  &&
                                                                        My_neigh_ring5[fifth_ring][index6] != second_ring &&
                                                                        My_neigh_ring5[fifth_ring][index6] != third_ring  &&
                                                                        My_neigh_ring5[fifth_ring][index6] != fourth_ring &&
                                                                        My_neigh_ring5[fifth_ring][index6] != fifth_ring  &&
                                                                        compare(ring6[first_ring], ring5[My_neigh_ring5[fifth_ring][index6]], 6, 5) &&
                                                                        compare_adjacant(ring5[fifth_ring], ring5[My_neigh_ring5[fifth_ring][index6]], 5, 5, ring6[first_ring])
                                                                        )
                                                                    {
                                                                        sixth_ring = My_neigh_ring5[fifth_ring][index6];
                                                                        
                                                                        if (N_ring5_neigh[sixth_ring] != 0 && N_ring5_neigh[sixth_ring] > 3)
                                                                        {
                                                                            for (int index7 = 0 ; index7 < N_ring5_neigh[sixth_ring]; index7++)
                                                                            {
                                                                                if (My_neigh_ring5[sixth_ring][index7] != first_ring  &&
                                                                                    My_neigh_ring5[sixth_ring][index7] != second_ring &&
                                                                                    My_neigh_ring5[sixth_ring][index7] != third_ring  &&
                                                                                    My_neigh_ring5[sixth_ring][index7] != fourth_ring &&
                                                                                    My_neigh_ring5[sixth_ring][index7] != fifth_ring  &&
                                                                                    My_neigh_ring5[sixth_ring][index7] != sixth_ring  &&
                                                                                    compare(ring5[second_ring], ring5[My_neigh_ring5[sixth_ring][index7]], 5, 5) &&
                                                                                    compare_adjacant(ring5[sixth_ring], ring5[My_neigh_ring5[sixth_ring][index7]], 5, 5, ring6[first_ring])
                                                                                    )
                                                                                {
                                                                                    seventh_ring = My_neigh_ring5[sixth_ring][index7];
 
                                                                                    for (int index8 = 0 ; index8 < N_ring6_ring5_neigh[first_ring]; index8++)
                                                                                    {
                                                                                        if(My_neigh_ring6_ring5[first_ring][index8] == seventh_ring)
                                                                                        {
                                                                                            if (compare(ring6[first_ring], ring5[seventh_ring], 6, 5))
                                                                                            {
                                                                                                
                                                                                                temp_ring1.clear();
                                                                                                temp_ring2.clear();
                                                                                                temp_ring3.clear();
                                                                                                temp_ring4.clear();
                                                                                                temp_ring5.clear();
                                                                                                temp_ring6.clear();
                                                                                                temp_ring7.clear();
                                                                                                
                                                                                                //put all the rings in temp_ring std::vectors.
                                                                                                for (int i_1 = 0; i_1 < 6; i_1++)
                                                                                                {
                                                                                                    temp_ring1.push_back(ring6[first_ring ][i_1]);
                                                                                                }
                                                                                                for (int i_1 = 0; i_1 < 5; i_1++)
                                                                                                {
                                                                                                    temp_ring2.push_back(ring5[second_ring][i_1]);
                                                                                                    temp_ring3.push_back(ring5[third_ring ][i_1]);
                                                                                                    temp_ring4.push_back(ring5[fourth_ring][i_1]);
                                                                                                    temp_ring5.push_back(ring5[fifth_ring ][i_1]);
                                                                                                    temp_ring6.push_back(ring5[sixth_ring ][i_1]);
                                                                                                    temp_ring7.push_back(ring5[seventh_ring][i_1]);
                                                                                                }
                                                                                                
                                                                                                repeated_atoms.clear();
                                                                                                set_intersection(temp_ring1.begin(),temp_ring1.end(),temp_ring2.begin(),temp_ring2.end(),back_inserter(repeated_atoms));    //check ring1 and ring2 for overlaps
                                                                                                
                                                                                                if (repeated_atoms.size() <= 2)
                                                                                                {
                                                                                                    repeated_atoms.clear();
                                                                                                    set_intersection(temp_ring2.begin(),temp_ring2.end(),temp_ring3.begin(),temp_ring3.end(),back_inserter(repeated_atoms));
                                                                                                    if (repeated_atoms.size() <= 2)
                                                                                                    {
                                                                                                        repeated_atoms.clear();
                                                                                                        set_intersection(temp_ring3.begin(),temp_ring3.end(),temp_ring4.begin(),temp_ring4.end(),back_inserter(repeated_atoms));
                                                                                                        if (repeated_atoms.size() <= 2)
                                                                                                        {
                                                                                                            
                                                                                                            repeated_atoms.clear();
                                                                                                            set_intersection(temp_ring4.begin(),temp_ring4.end(),temp_ring5.begin(),temp_ring5.end(),back_inserter(repeated_atoms));
                                                                                                            if (repeated_atoms.size() <= 2)
                                                                                                            {
                                                                                                                repeated_atoms.clear();
                                                                                                                set_intersection(temp_ring5.begin(),temp_ring5.end(),temp_ring6.begin(),temp_ring6.end(),back_inserter(repeated_atoms));
                                                                                                                if (repeated_atoms.size() <= 2)
                                                                                                                {
                                                                                                                    repeated_atoms.clear();
                                                                                                                    set_intersection(temp_ring6.begin(),temp_ring6.end(),temp_ring7.begin(),temp_ring7.end(),back_inserter(repeated_atoms));
                                                                                                                    
                                                                                                                    if (repeated_atoms.size() <= 2)
                                                                                                                    {
                                                                                                                        repeated_atoms.clear();
                                                                                                                        set_intersection(temp_ring7.begin(),temp_ring7.end(),temp_ring1.begin(),temp_ring1.end(),back_inserter(repeated_atoms));
                                                                                                                        if (repeated_atoms.size() <= 2)
                                                                                                                        {
                                                                                                                            
                                                                                                                            
                                                                                                                            temp_vec6[0] = first_ring ;
                                                                                                                            temp_vec6[1] = second_ring;
                                                                                                                            temp_vec6[2] = third_ring ;
                                                                                                                            temp_vec6[3] = fourth_ring;
                                                                                                                            temp_vec6[4] = fifth_ring ;
                                                                                                                            temp_vec6[5] = sixth_ring ;
                                                                                                                            temp_vec6[6] = seventh_ring ;
                                                                                                                            
                                                                                                                            
                                                                                            cup62512.push_back(temp_vec6);
                                                                                                                            
                                                                                                                            temp_ring1.clear();
                                                                                                                            temp_ring2.clear();
                                                                                                                            temp_ring3.clear();
                                                                                                                            temp_ring4.clear();
                                                                                                                            temp_ring5.clear();
                                                                                                                            temp_ring6.clear();
                                                                                                                            temp_ring7.clear();
                                                                                                                        }
                                                                                                                    }
                                                                                                                }
                                                                                                            }
                                                                                                        }
                                                                                                    }
                                                                                                }
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }           //end of finding half of 6(2)5(12) cages---------------------------------------------------------------------------
    
}

// New Function--------------------------------------------------------------------------------------

// New Function--------------------------------------------------------------------------------------

//Find Cages of 5(12) from cups of 5(6).

int CageFinderAlgorithm::cage_Finder(std::vector<std::vector<int>> cups, unsigned long count_cups, std::vector<std::vector<int>>& My_neigh_ring5, std::vector<std::vector<int>>& cage, std::vector<std::vector<int>>& cage_rings, std::string time )
{
    
    int count = 0;      //Counter for number of cages.
    std::vector<int> current_cup;
    int current_ring = 0 ;
    int N = 0;      //counter for number of common rings between side rings of current cup and side rings of rest of the cups.
    int M = 0;      //counter for number of side rings of current cup which have 2 common rings with side rings of rest of the cups.
    unsigned long int cup_size=0;
    if( cups.size() != 0 ){cup_size = cups[0].size();}
    std::vector<int> temp_vec, temp_vec_rings;
    
    
    for (int i = 0  ; i < count_cups ; i++)     //Loop over all the cups.
    {
        for (int j = 0 ; j < cups[i].size() ; j++ ) {current_cup.push_back(cups[i][j]);}        //Write current cup into a new variable.
        for (int k = i+1 ; k < count_cups ; k++)        //Loop over the rest of the cups to compare with current cup.
        {
            for(int l = 1 ; l < current_cup.size() ; l++)                           //Loop over side rings of current cup.
                                                                                    //"l" is the counter for side rings of current cup. It starts from 1
                                                                                    //because first ring in each cup is the base ring.
            {
                current_ring = current_cup[l];
                for(int m = 1 ; m < current_cup.size() ; m++)                          //Loop over side rings of the rest of the cups.
                {
                    for(int p = 0 ; p < My_neigh_ring5[current_ring].size() ; p++)      //Loop over neighbour rings of current ring.
                                                                                        //"p" is counter for all the neighbour rings of
                                                                                        //side rings of current cup.
                    {
                        if( My_neigh_ring5[current_ring][p] == cups[k][m] )
                        {N++;}
                    }
                }
                if(N==2){ M++; N = 0;}      //If current side ring has 2 neighbors in other cups, move on
                                            //to next side ring of the current cup.
                else {N = 0 ; break;}       //If current side ring has less than 2 neighbors, move on to the next cup.
            }
            if(M == cup_size -1)
            {
                
                for(int h = 0 ; h < cup_size  ; h++){temp_vec_rings.push_back(cups[i][h]);}
                for(int h = 0 ; h < cup_size  ; h++){temp_vec_rings.push_back(cups[k][h]);}
                cage_rings.push_back(temp_vec_rings);
                temp_vec_rings.clear();
                
                temp_vec.push_back(i);
                temp_vec.push_back(k);
                cage.push_back(temp_vec);
                temp_vec.clear();
                count++;
            }
            M = 0;
        }
        current_cup.clear();
    }
    return count;
    
}


int CageFinderAlgorithm::cage_Finder_64512(std::vector<std::vector<int>> cup62512, int count_62512_cups, std::vector<std::vector<int>>& cage_64512_rings)
{
    //cage_64512_rings:   This parameter holds the index of rings which make the 64512 cage in the following way:
    /*Each 64512 cage has 4 cups of 6156. So each cage will have 4 6-rings + 4*6 5-rings in total. These rings are written in cage_64512_rings this way: index 0,1,2,3 hold base rings which are hexagons and index 4-27 hold lateral rings which are pentagons. */
    
    int count_cage_64512 = 0;
    int N = 0;      //counter for number of shared lateral rings between any two cups.
    std::vector<std::vector<int>> cage_64512;     //cage_64512 holds the index of 62512 (or 6156) cups which make the 64512 cage.
    cage_64512.clear();
    int p;      //This holds the index of the 62512 cups. It is for simplicity of the code only.
    std::vector<int> temp,temp2;
    temp.clear();
    
    
    for (int i = 0 ; i < count_62512_cups - 1 ; i++)
    {
        temp.push_back(i);
        
        for (int j = i+1 ; j < count_62512_cups ; j++)
        {
            for (int k = 1 ; k < 7 ; k++)
            {
                for (int l = 1 ; l < 7 ; l++)
                {
                    if (cup62512[i][k] == cup62512[j][l])
                    {
                        N++;
                        if( N >= 2 ){temp.push_back(j);}
                    }
                }
            }
            N = 0;
        }
        if(temp.size() == 4)
        {
            cage_64512.push_back(temp);
            count_cage_64512++;
        }
        temp.clear();
    }
    
   /*To print 4 cups that form a 64512 cage uncomment this part.*/
//    for (int i = 0 ; i < cage_64512.size() ; i++)
//    {
//        if(cage_64512[i].size() >= 3 )
//        {
//            cout << "cup " << i << "  neighbors: " ;
//            for (int j = 0 ; j < cage_64512[i].size() ; j++) {cout << cage_64512[i][j] << " ";}
//            cout  << "\n";
//        }
//    }
    
    
    for (int i = 0 ; i < count_cage_64512 ; i++)
    {
        for (int j = 0 ; j < 4 ; j++ )
        {
            p = cage_64512[i][j];               //p is for simplicity of code (re-naming basically).
            temp.push_back(cup62512[p][0]);     //here temp holds the index of base rings (hexagons) of each of the four, 62512 cups.
        }
        for (int j = 0 ; j < 4 ; j++)           //this loops gets the pentagons of the four cups of cage_64512 into cage_64512_rings.
        {
            p = cage_64512[i][j];
            for(int k=1;k<=6;k++){temp.push_back(cup62512[p][k]);}
        }
        cage_64512_rings.push_back(temp);       //this line puts the 24 lateral rings of the 4 cups in cage_64512_rings to fill indices 4-27.
        temp.clear();
    }
    
    
    
    return count_cage_64512;
}

