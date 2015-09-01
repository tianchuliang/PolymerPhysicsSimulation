//
//  polymer.h
//  PolSim
//
//  Created by Tianchu Liang on 10/31/14.
//  Copyright (c) 2014 Tianchu Liang. All rights reserved.
//

#ifndef PolSim_polymer_h
#define PolSim_polymer_h

#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include "vertex.h"

using namespace std;

const float hb_potential_depth=1.;//potential depth for the hb bonding well
const float exv_potential_depth=1.;//depth for the exv energy well.

const float r_0 = 1.0; //radius at which both slj and tlj achieve minima.

const float r_co_slj=r_0; // cutoff for exv slj
const float sigma=0.89089;
const float r_co_tlj=2.5*sigma; // cutoff for hb

const float du_tlj=hb_potential_depth*0.0163; // offset for the truncated lennard jones
const float stepSize=1.;
const float temperature=1.;//really is k_B*T
const float springConst=0.5;//For linking energy: we have it as: U=3*k_B*T/(2*b^2)*|ri-ri-1|^2=3/(2b^2)*|r-r|^2

/*
 Description of polymer class:
 
 The polymer class has one protected attribute, that is a vector of vertices, called monomers. 
 Each vertex stands for a monomer, containing monomer's spatial location and its hydrogen bond status. 
 
 The polymer class has a default constructor and a normal constructor. The normal constructor takes an int,N, as the only parameter, which tells the how many monomers this polymer is going to have. 
 
 There are 7 operations, as shown below. Among 7 operations, dU(), step(), getMonomer() all require careful treatment to the hb attribute of each vertex. Detailed implementation can be found at polymer.cpp file. */


class polymer
{
    
protected:
    vector<vertex> monomers; // r_i
    
public:
    polymer();//Default Constructor
    
    polymer(int N, bool link, bool exv, bool hbond);//Constructor; three boolean variables are switches for which
    //energy components are to be taken into account. 
    
    //Calculate the total energy of the entire polymer:
    float U_tot();
    
    //Calculate the change in energy for each single vertex movement:
    //The change in energy comes from:
    //Change in U_link
    //Change in U_exv
    float dU (vertex *beforeMove, vertex *afterMove, int itsIndex);
    
    //move one monomer from its previous position to another new position. This movement guarantees that
    //the monomer will not collide with other monomers.
    vertex make_a_legal_move(vertex *beforeMove, int itsIndex);
    
    //Make a change to the monomer at (index), that obeys Boltzmann distribution.
    void step(int index);
    
    //Gibbs sampling:
    void Gibbs_step(int start, int end);
    
    //Return a specific monomer's information in the format of a vertex.
    vertex getMonomer(int index);
    
    bool check_exv(vertex toBeChecked, int avoidIndex);
    
    //Return the address of the entire polymer.
    vector<vertex>* getMonomers();
    
    float delta_hbond(vertex *beforeMove, vertex *afterMove);
    
    void breakBond(int i);
    
    void formBond (int i, int j);
    
    int itsAt(vertex *ThisOne);
    
    int findOne(vector<int> indices, vertex *ThisVertex);
    
    vector<int> formOrNot(int i);
    
};







#endif
