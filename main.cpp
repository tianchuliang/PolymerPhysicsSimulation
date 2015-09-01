//  main.cpp
//  PolSim
//
//  Created by Tianchu Liang on 10/31/14.
//  I hereby reaffirm the Lawrence University Honor Code.
//  Copyright (c) 2014 Tianchu Liang. All rights reserved.
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include "vertex.h"
#include "polymer.h"
/*
The main.cpp file achieves the following tasks:
1. Construct the running_avg object, which is responsible for calculating running averages.
*/
const int numbers[]={16,32,64,128,256,512};
const int limit=6;
const int NumOfRuns=30;
const bool OutPutPosits=false;
const bool OutPutRgAvg=true;
const bool link=true;
const bool exv=true;
const bool hbond=false;
const int numOfGibbs=100;
const int numOfRounds=100;
using namespace std;

//Running Average of Rg class:
class running_avg
{
    float avg;
    int count;
public:
    running_avg(){
        avg=0.;
        count=0;
    }
    float operator()(float next_avg){   
        avg=(next_avg+count*avg)/(count+1);
        count++;
        return avg;
    }
};

//Calculate the mean position of all monomers:
vertex meanPosition(vector<vertex> *monomers)
{
    float meanx=0.;
    float meany=0.;
    float meanz=0.;
    //One for loop to iteratively calculate mean values of x,y,z:
    for (int i=0; i<monomers->size(); i++) {
        //        meanx += monomers->at(i).x;
        //        meany += monomers->at(i).y;
        //        meanz += monomers->at(i).z;
        meanx=(meanx*i+monomers->at(i).x)/(i+1);
        meany=(meany*i+monomers->at(i).y)/(i+1);
        meanz=(meanz*i+monomers->at(i).z)/(i+1);
    }
    //    meanx /= static_cast<float>(monomers->size());
    //    meany /= static_cast<float>(monomers->size());
    //    meanz /= static_cast<float>(monomers->size());
    vertex mean=vertex(meanx,meany,meanz,NULL);
    //cerr<<meanx<<" "<<meany<<" "<<meanz<<endl;
    return mean;
}

//Calculate Running Avg R_g using a functor:
//Calculation of Radius of Gyration:
float R_g(vector<vertex> *currentMonomers)
{
    //With each Gibbs sampling, the monomers positions are all changed,
    //So, before calculating running average on R_g, we need to update the mean position:
    vertex meanPos=vertex();
    meanPos=meanPosition(currentMonomers);
    
    //Calculate the summation of (monomer_i_position - monomer_mean_position)^2 of i=0 to total.
    float NewSumR_g=0.;
    for (int i =0; i<currentMonomers->size(); i++) {
        vertex DeltaPos=vertex();
        DeltaPos=currentMonomers->at(i)- meanPos;
        NewSumR_g=NewSumR_g+(DeltaPos.dist())*(DeltaPos.dist());
    }
    
    //newRg is NewSumR_g divided by N:
    float newRg;
    newRg=NewSumR_g/currentMonomers->size();
    
    //cerr << "--" << newRg << endl;
    //meanPos = currentMonomers->at(9);
    //cerr<< meanPos.x << " " << meanPos.y << " " << meanPos.z <<endl;
    return sqrt(newRg);
}


void writeOutposits(string fileName, int NumOfGibbs, int Rounds)
{
    ofstream out;
    out.open(fileName.c_str(),ios::out);
        
        for (int i=0; i<limit; i++) {
            
            cerr<<"Processing polymer with "<<numbers[i]<<"monomers."<<endl;
            //Set up the polymer:
            polymer pol=polymer(numbers[i],link,exv,hbond);
            
            //Write out initial positions before burn in:
            for (int j =0; j<numbers[i];j++) {
                out<<pol.getMonomer(j).x<<" "<<pol.getMonomer(j).y<<" "<<pol.getMonomer(j).z<<endl;
            }
            
            //Do burn in:
            for (int j =0; j<1000; j++) {
                
                //cerr<<"doing burning in!"<<endl;
                pol.Gibbs_step(0, numbers[i]);
                
                
                //Write-out data during burn in:
                out<<"==="<<endl;
                for (int k =0; k<numbers[i];k++) {
                    out<<pol.getMonomer(k).x<<" "<<pol.getMonomer(k).y<<" "<<pol.getMonomer(k).z<<endl;
                }
                
            }
            
            cerr<<"..."<<endl;
            int index=0;
            
            while (index!=Rounds)
            {
                for (int k=0; k<NumOfGibbs; k++) {
                    pol.Gibbs_step(0, numbers[i]);
                }
                //cerr<<"writing output"<<endl;
                out<<"==="<<endl;
                for (int a =0; a<numbers[i];a++) {
                    out<<pol.getMonomer(a).x<<" "<<pol.getMonomer(a).y<<" "<<pol.getMonomer(a).z<<endl;
                }
            index++;
            }
        }
}

void writeOutRgAvg(string fileName, int NumOfGibbs, int Rounds)
{
    ofstream out;
    out.open(fileName.c_str(),ios::out);
    
    for (int b=0; b<NumOfRuns; b++) {
        
        cerr<<b+1<<"th run"<<endl;
        
        float avgs[limit];
        for (int a=0; a<limit; a++) {
            avgs[a]=0.;
        }
        
        for (int i=0; i<limit; i++) {
            cerr<<"Processing polymer with "<<numbers[i]<<"monomers."<<endl;
            
            //Set up the polymer:
            polymer pol=polymer(numbers[i],link,exv,hbond);
            
            //Burn in
            float avgRg=0.;
            //cerr<<"about to burn in"<<endl;
            for (int j =0; j<1000; j++) {
                //cerr<<"doing burning in!"<<endl;
                pol.Gibbs_step(0, numbers[i]);
            }
            
            cerr<<"..."<<endl;
            //Starts 100 gibbs step sampling:
            running_avg rg_ravg;
            int index=0;
            while (index!=Rounds) {
                
                for (int k=0; k<NumOfGibbs; k++) {
                    pol.Gibbs_step(0, numbers[i]);
                }

                avgRg=rg_ravg(R_g(pol.getMonomers()));
                index++;
            }
            avgs[i]=avgRg;
        }
        
            for (int i=0; i<limit; i++) {
                 out<<numbers[i]<<" "<<avgs[i]<<endl;
            }
    }

}



//Main method:
int main(int argc, const char * argv[]) {
    srand ((unsigned int)time(NULL));
    const string OUTPUT_PROMPT = "Output file name:  ";
    
    string out_name;
    cerr<<OUTPUT_PROMPT<<endl;
    cin>>out_name;
    
    if(OutPutPosits){
        writeOutposits(out_name,numOfGibbs,numOfRounds);
    }
    
    if(OutPutRgAvg){
        writeOutRgAvg(out_name,numOfGibbs,numOfRounds);
    }
    
    return 0;
}




//The following main method specfically tests the hydrogen bond attributes:
//int main(int argc, const char * argv[])
//{
//    polymer pol=polymer(8);
//
//    cerr<<"Attributes of 5 monomers are:"<<endl;
//    for (int i =0; i<8; i++) {
//            cerr<<pol.getMonomer(i).x<<" "<<pol.getMonomer(i).y<<" "<<pol.getMonomer(i).z<<" "<<pol.getMonomer(i).hb<<endl;
//    }
//
//    for (int i =0; i<5; i++) {
//        pol.Gibbs_step(0, 8);
//
//        cerr<<"After the "<<i<<"th Gibbs step, the coordinates are:"<<endl;
//        for (int i =0; i<8; i++) {
//            cerr<<pol.getMonomer(i).x<<" "<<pol.getMonomer(i).y<<" "<<pol.getMonomer(i).z<<" "<<pol.getMonomer(i).hb<<endl;
//        }
//    }
//
//
//
//
//    return 0;
//}




