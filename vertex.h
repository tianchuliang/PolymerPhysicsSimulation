//
//  vertex.h
//  PolSim
//
//  Created by Tianchu Liang on 10/31/14.
//  Copyright (c) 2014 Tianchu Liang. All rights reserved.
//

#ifndef PolSim_vertex_h
#define PolSim_vertex_h

#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
using namespace std;

/*
 Description of vertex class:
 
 The vertex class is an object that has four attributes:
 float x, y, z as the 3-D coordinates, and 
 vertex *hb, which is a pointer to another vertex.
 
 The hb attribute points to another vertex that is bonded with 'this' vertex, through hydrogen bond. 
 
 The vertex class has aa Default Constructor as well as a normal constructor.
 
 There are 12 operations in the vertex class. 
 The first 11 operations all deal with the coordinates of vertices, either overloading certain operations (+,-,*,/,etc) or computing the distance between two vertices.
 
 Notice that the first 11 operations don't do anything to the hb attribute. This is because we don't want to mix coordinate-related calculations together with pointer manipulations. Most hb (pointer) manipulations are done in the polymer.cpp file, where many operations (such as breaking or forming hb bonds) require changing the hb attribute.
 
 Sometimes, it is very helpful to get the entire vertex information including the hb attribute. The overloaded "=" operator does not quite finish the job, since it only equates a vertex's coordinate with another vertex's coordinate. Therefore, the last operation is copy(), which returns a vertex's full information. 
 */

class vertex
{
public:
    float x;
    float y;
    float z;
    int hb;
    
    
    vertex(){}//Default Constructor
    
    vertex(float a, float b, float c, int hb);

    
    //Overloading the "+" operator, does not do any change to *hb
    vertex operator+ (const vertex& other);
    
    //Overloading the "-" operator, does not do any change to *hb
    vertex operator- (const vertex& other);
    
    //Overloading the "-" (unary), does not do any change to *hb
    vertex operator- ();
    
    //Overloading * as scalar product, does not do any change to *hb
    vertex scalarProd (int scalar, vertex *other);
    
    //Overloading * as inner product, does not do any change to *hb
    float inner_Prod(vertex *A, vertex *B);
    
    //Overloading / as scalar division, does not do any change to *hb
    vertex scalarDiv(int scalar, vertex *other);
    
    //Unit_Vec, gives you the unit vector (just the parallel vector of unit length), does not do any change to *hb
    vertex Unit_Vec(vertex *A);
    
    //dist(), gives the length between two vertices, does not do any change to *hb
    //For example, the distance between vertices A and B is (A-B).dist().
    float dist();
    
    //Cross Product, does not do any change to *hb
    vertex Cross_Prod(vertex *A, vertex *B);
    
    //Overloading the "=" operator, does not do any change to *hb
    void operator= (const vertex& other);
    
    //Overloading the "==" operator, does not do any change to *hb
    bool operator== (const vertex& other);
    
    
};






#endif

