//
//  polymer.cpp
//  PolSim
//
//  Created by Tianchu Liang on 10/31/14.
//  Copyright (c) 2014 Tianchu Liang. All rights reserved.
//

#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include "polymer.h"
#include "vertex.h"

using namespace std;

/*
 There are some auxiliary functions before class member functions of polymer.h.
 These auxiliary functions are (in order):
 
 isInShell(), help generating spatially 'not so loose and not so tight' coordinates.
 getPos(), based on the previous vertex, generate the next possible random coordinate.
 getPos_chain(), generate vertex coordinates in a linear, chain fashion.
 U_link(), calculate linking energy between two monomers.
 U_exv(), calculate exclusion volume energy between two monomers.
 U_hb, calculate hydrogen bond energy between two bonded monomers.
 isValid_Candidate(), For initializing monomers: Check if the proposed new monomoer member is in anybody's way. Return true if the new vertex is in nobody's way.
 Accpt_or_Rejt, accept or reject the movement of one monomer, based on whether this movement obeys Boltzmann distribution.
 delta_hbond(), calculate change in hb for a specific pair of monomers.
 breakOrNot(), check whether a certain deltaU_hb value means we need to break the bond.
 formOrNot(), check if one vertex is in the viscinity of any other monomer(s) such that a hydrogen bond can form.
 ...
 ...
 ...
 ...
 ...
 UNFINISHED
 ...
 ...
 ...
 ...
 ...
 
 
 */

bool link_E=true;
bool exv_E=true;
bool hbond_E=true;



//Check, whether the randomly generated
//vertex is in the acceptable spatial shell:
bool isInShell(float a, float b, float c, float rin, float rout)
{
    return (rin<sqrtf(a*a+b*b+c*c)<rout);
}

//Get new vertex by generating random 3-D coordinates:
vertex getPos(vertex *prevVert)
{
    float rout;
    rout=2.5;
    float rexv;
    
    //If exclusion volume is taken into account, don't allow phantom monomers.
    //Else, allow phantom monomers.
    if (exv_E) {
        rexv=r_0;
    }else{
        rexv=0.;
    }
    
    // change to use RAND_MAX
    float xmove;
    float ymove;
    float zmove;
 
    do{
        xmove=(2.*rout*rand())/RAND_MAX-rout;
        ymove=(2.*rout*rand())/RAND_MAX-rout;
        zmove=(2.*rout*rand())/RAND_MAX-rout;
    }while(!isInShell(xmove,ymove,zmove,rexv,rout));
    
    vertex newVert=vertex();
    
    newVert.x=prevVert->x+xmove;
    newVert.y=prevVert->y+ymove;
    newVert.z=prevVert->z+zmove;
    
    return newVert;
    
}

//Get new vertex like a chain:
vertex getPos_chain(vertex *prevVert)
{
    float rout;
    rout=2.5;
    
    float xmove;
    float ymove;
    float zmove;
    
    xmove=rout;
    ymove=0.;
    zmove=0.;
    
    //cerr<<prevVert->x<<" "<<prevVert->y<<" "<<prevVert->z<<endl;
    vertex newVert = vertex();
    newVert.x=prevVert->x+xmove;
    newVert.y=prevVert->y+ymove;
    newVert.z=prevVert->z+zmove;
    //cerr<<newVert.x<<" "<<newVert.y<<" "<<newVert.z<<endl;
    return newVert;
    
}

//Calculate energies of monomers:
float U_link(vertex *A, vertex *B)
{
    float result;
    float r;
    r=((*A)-(*B)).dist();
    result=(1./2.)*springConst*r*r;
    return result;
}

float U_exv(vertex *A, vertex *B)
{
    //U(r)=piecewise{{U_0*((r_0/r)^12-2*(r_0/r)^6)+U_0,{r<r_0}},{0,{r>=r_0}}}
    float U_0=exv_potential_depth;
    float r;
    r=((*A)-(*B)).dist();
    float result;
    
    //r_co_slj=r_0
    if (r>=r_co_slj) {
        result=0.;
    }
    else{
        float x = r_0/r;
        x = x*x*x*x*x*x;
        result=U_0*(x*x-2*x)+U_0;
    }
    return result;
}

float U_hb(vertex *A, vertex *B)
{
    float U_0=hb_potential_depth;//depth of the well
    float r;
    r=((*A)-(*B)).dist();
    float result;
    
    //r_co_tlj=2.5*sigma
    if (r>=r_co_tlj) {
        result=0.;
    }
    else{
        float x = r_0/r;
        x = x*x*x*x*x*x;
        result=U_0*((x*x-2*x)+du_tlj);
    }
    return result;
}

//For initializing monomers:
//Check if the proposed new monomoer member is in anybody's way.
//Return true if the new vertex is in nobody's way.
bool isValid_Candidate(vector<vertex> * Current, vertex * candidate)
{
    
    //If exclusion volume is taken into account, don't allow phantom monomers.
    //Else, allow phantom monomers.
    float rexv;
    if (exv_E) {
        rexv=r_0;
    }else{
        rexv=0.;
    }
    
    for (int i = 0; i<Current->size(); i++) {
        
        
        if (((Current->at(i))-(*candidate)).dist()<=rexv) {
            //cerr<<"Unacceptable"<<endl;
            return false;
        }
        
    }
    //cerr<<"Acceptable"<<endl;
    return true;
}

//Accept or reject the move of one monomer based on change in Energy:
bool Accpt_or_Rejt(float dU )
{
    
    float ratio;
    ratio=exp(-(dU)/temperature);
    
    //cerr<<"value of exp(-(dU)/temperature): "<<(ratio)<<endl;
    
    //take the ratio as the percent for accepting
    
    float random;

    random=rand();
    
    float percentage;
    percentage=random/RAND_MAX;
    
    if (ratio>1)
    {
//        cerr<<"1";
        return true;
        
    }
    else{
        float percentage;
        percentage=(1.*rand())/RAND_MAX;
        
        //cerr<<"initially rejected."<<endl;
        //cerr<<"percentage to be compared with exp(-(dU)/temperature): "<<percentage<<endl;
        if (percentage<ratio)
        {
  //          cerr<<"2";
            return true;
        }
    }
  //  cerr<<"0";
    return false;
}

//Calculate change in hb for a specific pair of monomers:
float polymer:: delta_hbond(vertex *beforeMove, vertex *afterMove)
{
    float result;

        if (beforeMove->hb!=-1) {
            
            float hb1;
            float hb2;
            
            hb1=U_hb(beforeMove, &monomers.at(beforeMove->hb));
            hb2=U_hb(afterMove, &monomers.at(afterMove->hb));
            
            result=hb1-hb2;
            
        }else
        {
            result=0.;
        }
  
    
    return result;
}

//Check whether a certain deltaU_hb value means we need to break the bond.
bool breakOrNot(float deltaU_hb){
    
    float P_hb;
    P_hb=exp(-deltaU_hb/(temperature));
    
    float random;
    random=rand();
    random=random/RAND_MAX;
    
    if (random<P_hb) {
        return true;
    }
    
    return false;
}

//
void polymer:: breakBond(int ToBeBroken)
{
    //cerr<<"breaking bond between "<<ToBeBroken<<"th monomer and"<<monomers.at(ToBeBroken).hb<<"th monomer."<<endl;
    //Cut the other end first:
    monomers.at(monomers.at(ToBeBroken).hb).hb=-1;
    //Then we cut this end:
    monomers.at(ToBeBroken).hb=-1;
}

//
void polymer:: formBond (int thisOne, int thatOne){
    
    monomers.at(thisOne).hb= thatOne;
    monomers.at(thatOne).hb= thisOne;
    
}

//Given a current monomer, find a possible monomer that is the closest for hb to happen.
int polymer:: findOne(vector<int> indices, vertex *ThisVertex)
{
    float minDist=INFINITY;
    float minIndex=0;
    
    vertex temp=vertex();
    
    for (int a =0; a<indices.size(); a++) {
        temp=*ThisVertex - monomers.at(indices.at(a));
        
        if (temp.dist()<minDist) {
            minDist=temp.dist();
            minIndex=indices.at(a);
        }
    }
    
    return minIndex;
}

//Check if one vertex is in the viscinity of any other monomer(s) such that a hydrogen bond can form.
vector<int> polymer:: formOrNot(int thisOne)
{
    vector<int> PossibleBonding;
    vertex temp=vertex();
    
    for (int i =0; i< monomers.size(); i++) {
        
        //Check if this monomer already has a bond!!!!!!
        if (monomers.at(i).hb==-1 && !(monomers.at(i)==monomers.at(thisOne))) {
            
            temp=monomers.at(i)-monomers.at(thisOne);
            float dist=temp.dist();
            //            cerr<<"temp.dist()="<<dist<<endl;
            //            cerr<<"r-cutoff is"<<r_co_tlj<<endl;
            //If we do form a bond, return the possible index(es) of the monomer, to which ThisVertex is going
            //to bond with.
            
            if (dist<r_co_tlj) {
                PossibleBonding.push_back(i);
            }
        }
    }
    
    //If we don't form bond, return the empty PossibleBonding
    return PossibleBonding;
}



//Polymer constructor:
polymer::polymer(int N, bool link, bool exv, bool hbond)
{
    
    link_E=link;
    exv_E=exv;
    hbond_E=hbond;
    cerr<<"==="<<link<<"==="<<exv<<"==="<<hbond<<"==="<<endl;
    
    
    //cerr<<"Constructing Polymer..."<<endl;
    //Initialize first monomer, with (0,0,0) as coordinate, and NULL for hb.
    vertex initial=vertex(0.,0.,0., -1);
    monomers.push_back(initial);
    
    //Constructing other monomer's coordinate based on the first monomer, using getPos().
    for (int i =0; i<N-1; i++) {
        
        //cerr<<"current monomer's end is:"<<(*(--monomers.end())).x<<" "<<(*(--monomers.end())).y<<" "<<(*(--monomers.end())).z<<" "<<endl;
        vertex candidate=vertex();
        
        candidate=getPos_chain(&(*(--monomers.end())));
        
        while (!isValid_Candidate(&monomers, &candidate)) {
            
            candidate=getPos(&(*(--monomers.end())));
            
        }
        
        //Each monomer is initially NULL in the hb attribute.
        candidate.hb=-1;
        
        //cerr<<"============"<<endl;
        //cerr<<candidate.x<<" "<<candidate.y<<" "<<candidate.z<<" "<<endl;
        monomers.push_back(candidate);
        
    }
    
}

vertex polymer::getMonomer(int index){
    vertex result=monomers.at(index);
    return result;
}


////Calculate the total energy of the entire polymer:
//float polymer:: U_tot()
//{
//    float link=0;
//    float exv=0;
//    
//    //Sum up all the U_link from all monomers:
//    if (link_E) {
//        
//        for (int i =0; i<monomers.size()-1; i++)
//        {
//            link=link+U_link(&(monomers.at(i)), &(monomers.at(i+1)));
//        }
//    }
//
//    if (exv_E) {
//        //Sum up all the U_exv from all monomers:
//        
//        for (int j=0; j<monomers.size()-1; j++) {
//            for (int k=j+1; k<monomers.size(); k++) {
//                exv=U_exv( &(monomers.at(j)), &(monomers.at(k)));
//            }
//        }
//    }
//    
//    
//    /*Following code for calculating total hb energy may not work correctly.
//     Currently, the program doesn't call U_tot().
//     Sum up all the U_hb from all monomers:*/
//    
//    float result;
//    result=link+exv;
//    return result;
//    
//}




//Calculate the change in energy for each single vertex movement:
//The change in energy comes from:
//Change in U_link
//Change in U_exv
//Change in U_hb
float polymer:: dU (vertex *beforeMove, vertex *afterMove, int itsIndex)
{
    
    //Calculate U_link energy change:
    float dU_link1=0;//Energy change in the left link.
    float dU_link2=0;//Energy change in the right link.
    
    if (link_E) {
        //cerr<<"linkE is considered."<<endl;
        if (itsIndex==0) {
            dU_link1=0;//When the first vertex gets moved, its energy change in left link is 0.
            dU_link2=U_link(&(monomers.at(itsIndex+1)), afterMove)-U_link(&(monomers.at(itsIndex+1)), beforeMove);
        }
        else if(itsIndex==monomers.size()-1)
        {
            dU_link1=U_link(&(monomers.at(itsIndex-1)), afterMove)-U_link(&(monomers.at(itsIndex-1)), beforeMove);
            dU_link2=0;//When the last vertex gets moved, its energy change in the right link is 0.
        }
        else
        {
            //Energy change in the left link:
            dU_link1=U_link(&(monomers.at(itsIndex-1)), afterMove)-U_link(&(monomers.at(itsIndex-1)), beforeMove);
            
            //Energy change in the right link:
            dU_link2=U_link(&(monomers.at(itsIndex+1)), afterMove)-U_link(&(monomers.at(itsIndex+1)), beforeMove);
        }
    }
    
    float dU_link;
    dU_link=dU_link1+dU_link2;
    
    //Calculate change in U_exv energy:
    float dU_exv;
    dU_exv=0;
    float U_exv_previous=0.;
    float U_exv_after=0.;
    if (exv_E) {
        //cerr<<"exvE is considered."<<endl;
        for (int i =0; i<monomers.size(); i++) {
            
            if (i!=itsIndex) {
                
                U_exv_previous=U_exv_previous+U_exv(beforeMove, &(monomers.at(i)));
                U_exv_after=U_exv_after+U_exv(afterMove, &(monomers.at(i)));

//                float add=0.;
//                add=U_exv(afterMove, &(monomers.at(i)))-U_exv(beforeMove, &(monomers.at(i)));
//                dU_exv=dU_exv+add;
            }
            else
            {
                //Do nothing
                //Skip comparing U_exv with it self.
            }
        }
    }
    dU_exv=U_exv_after-U_exv_previous;
    
    
    //Calculate change in U_hb:
    float dU_hb=0;
    
    if (hbond_E) {
        //cerr<<"hbondE is considered."<<endl;
        dU_hb=delta_hbond(beforeMove, afterMove);
    }
    
    //cerr<<"===="<<endl;
    //cerr<<"dU_exv, dU_link, and dU_hb are:"<<endl;
    //cerr<<dU_exv<<" "<<dU_link<<" "<<dU_hb<<endl;
    float result;
    result=dU_link+dU_hb+dU_exv;
    return result;
}
/*The following may not be useful; because when making movements, all the constraints should be 
handled by the thermal energy constraint.
////During making new movement:
////Check if exclusion volume is violated:
//bool polymer:: check_exv(vertex toBeChecked,int avoidIndex){
//    
//    //If exclusion volume is taken into account, don't allow phantom monomers.
//    //Else, allow phantom monomers.
//    float rexv;
//    if (exv_E) {
//        rexv=r_0;
//    }else{
//        rexv=0.;
//    }
//    
//    for (int i = 0; i<monomers.size(); i++) {
//        if (i!=avoidIndex) {
//            if (((monomers.at(i))-(toBeChecked)).dist()<=rexv) {
//                
//                //cerr<<((monomers.at(i))-(toBeChecked)).dist()<<"-------"<<rexv<<endl;
//                //cerr<<"Unacceptable"<<endl;
//                return false;
//            }
//            //cerr<<((monomers.at(i))-(toBeChecked)).dist()<<"========"<<endl;
//        }else{
//        //Do nothing, skip comparing with itself.
//        }
//    }
//    //cerr<<"***"<<rexv<<"***"<<endl;
//    //cerr<<"Acceptable"<<endl;
//    return true;
//}
 */

//Make a move such that the new position does not collide with any other monomers.
vertex polymer:: make_a_legal_move(vertex *beforeMove, int itsIndex)
{
    vertex result=vertex();
    //vertex temp=vertex();
    
    float xmove;
    float ymove;
    float zmove;
    
    do{
        //temp=*beforeMove;
        xmove=(2.*stepSize*rand())/RAND_MAX-stepSize;
        
        ymove=(2.*stepSize*rand())/RAND_MAX-stepSize;
        
        zmove=(2.*stepSize*rand())/RAND_MAX-stepSize;
        
        result = vertex( xmove, ymove, zmove,-1);
        
        //temp=result+temp;
        
    }while( (result.dist() < stepSize));

    //cerr<<xmove<<" "<<ymove<<" "<<zmove<<endl;
    result=result+*beforeMove;
     
    //The result inherits the hb attribute from the beforeMove.
    result.hb=beforeMove->hb;
    return result;
}

//Step
void polymer:: step(int index)
{
    //cerr<<"=================================="<<endl;
    //cerr<<"start stepping!"<<endl;
    //Move vertex to a new position, which does not
    //collide with any other existing vertices:
    
    vertex beforeMove=vertex();
    beforeMove=monomers.at(index);
    
    //cerr<<"before steping the "<< index <<"th monomer"<<", "<<"its original coordinates are: "<<endl;
    //cerr<<beforeMove.x<<" "<<beforeMove.y<<" "<<beforeMove.z<<endl;
    
    
    //Give aftermove the new position;
    //NOTE: aftermove will inherit the hb attribute from beforeMove.
    vertex afterMove=vertex();
    afterMove=make_a_legal_move(&beforeMove, index);
    
    //cerr<<"A possible new position is:"<<endl;
    //cerr<<afterMove.x<<" "<<afterMove.y<<" "<<afterMove.z<<endl;
    
    //calculate dU.
    float deltaU;
    deltaU=dU(&beforeMove, &afterMove, index);
    
    //cerr<<"Change in Energy is:"<<endl;
    //cerr<<deltaU<<endl;
    
    //Accept or reject
    //cerr<<Accpt_or_Rejt(deltaU)<<endl;
    if (Accpt_or_Rejt(deltaU)) {
        //cerr<<"Change is Accepted!"<<endl;
        monomers.at(index)=afterMove;
    }
    
    //cerr<<"afterMove end up: "<<endl;
    //cerr<<monomers.at(index).x<<" "<<monomers.at(index).y<<" "<<monomers.at(index).z<<endl;
    
}

//Auxillary method that returns the index of a given monomer.
int polymer:: itsAt(vertex *ThisOne)
{
    int result = -5;
    
    for (int i =0; i<monomers.size(); i++) {
        if (*ThisOne==monomers.at(i)) {
            result=i;
        }
    }
    
    return result;
}

void polymer::Gibbs_step(int start, int end)
{
    vertex previous=vertex();
    vertex after=vertex();
    
    //cerr<<"Gibbs stepping "<<end<<"monomers"<<endl;
    for (int i=start; i<end; i++) {
        previous=monomers.at(i);
        step(i);
        after=monomers.at(i);
        
        //If hbond interaction is considered, hb breaking/formation:
        if (hbond_E) {
            if (monomers.at(i).hb!=-1) {
                
                //cerr<<i<<"th monomer has an existing hbond with"<<monomers.at(i).hb<<"th monomer."<<endl;
                
                float deltaU_hb;
                deltaU_hb=delta_hbond(&previous, &after);
                
                if (breakOrNot(deltaU_hb)) {
                    cerr<<"breaking bond"<<endl;
                    breakBond(i);
                }
            }else{
                //cerr<<i<<"th monomer does not have an existing hbond"<<endl;
                
                vector<int> PossibleCandidates=formOrNot(i);
                
                if (PossibleCandidates.size()>0) {
                    int whichOne=findOne(PossibleCandidates,&monomers.at(i));
                    //cerr<<i<<"th monomer forms bond with"<<whichOne<<"th monomer."<<endl;
                    cerr<<"forming bond"<<endl;
                    formBond(i, whichOne);
                }
            }
        }
    }
}


vector<vertex> * polymer::getMonomers()
{
    return &monomers;
}
