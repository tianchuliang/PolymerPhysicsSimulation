//
//  vertex.cpp
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
#include "vertex.h"
using namespace std;


//Normal Constructor:
vertex:: vertex (float a, float b, float c, int d)
{
    x=a;
    y=b;
    z=c;
    hb=d;
}

//Overloading the "+" operator, it does not do anything to the hb attribute.
vertex vertex::operator+ (const vertex& other)
{
    vertex result=vertex();
    
    result.x=this->x+other.x;
    result.y=this->y+other.y;
    result.z=this->z+other.z;
    
    return result;
}

//Overloading the "-" operator, it does not do anything to the hb attribute.
vertex vertex:: operator- (const vertex& other)
{
    vertex result=vertex();
    
    result.x=this->x-other.x;
    result.y=this->y-other.y;
    result.z=this->z-other.z;
    
    return result;
}

//Overloading the "-" (unary), it does not do anything to the hb attribute.
vertex vertex:: operator- ()
{
    vertex result=vertex(-x,-y,-z,-1);
    return result;
}

//Overloading * as scalar product, it does not do anything to the hb attribute.
vertex vertex:: scalarProd (int scalar, vertex *other)
{
    vertex result = vertex (scalar*other->x, scalar* other->y, scalar* other->z,-1);
    return result;
}

//Overloading * as inner product, it does not do anything to the hb attribute.
float vertex:: inner_Prod(vertex *A, vertex *B)
{
    float result=(A->x)*(B->x)+(A->y)*(B->y)+(A->z)*(B->z);
    return result;
}

//Overloading / as scalar division, it does not do anything to the hb attribute.
vertex vertex:: scalarDiv(int scalar, vertex *other)
{
    vertex result= vertex((other->x)/scalar, (other->y)/scalar,(other->z)/scalar, -1);
    return result;
}

//Unit_Vec, gives you the unit vector (just the parallel vector of unit length), it does not do anything to the hb attribute.
vertex vertex:: Unit_Vec(vertex *A)
{
    float length;
    length=sqrt((A->x)*(A->x)+(A->y)*(A->y)+(A->z)*(A->z));
    
    vertex result=vertex(A->x/length,A->y/length,A->z/length,-1);
    return result;
}



//dist, gives the length between two verterx, it does not do anything to the hb attribute.
float vertex:: dist(){
    return sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
}

//Cross Product, it does not do anything to the hb attribute.
vertex vertex:: Cross_Prod(vertex *A, vertex *B)
{
    vertex result=vertex((A->y)*(B->z)-(A->z)*(B->y),(A->z)*(B->x)-(A->x)*(B->z),(A->x)*(B->y)-(A->y)*(B->x),-1);
    return result;
}

//Overloading the "=" operator, it does not do anything to the hb attribute.
void vertex:: operator= (const vertex& other)
{
    this->x=other.x;
    this->y=other.y;
    this->z=other.z;
    this->hb=other.hb;
    
}

//Overloading the "==" operator, it does not do anything to the hb attribute.
bool vertex:: operator== (const vertex& other)
{
    
    if (this->x==other.x && this->y==other.y &&this->z==other.z)
    {
        return true;
    }
    
    return false;
}

