//
// Created by s2533081 on 2/21/19.
//

/**
#ifndef SMOKE_QUAD_H
#define SMOKE_QUAD_H


#include <GL/gl.h>
#include <rfftw.h>              //the numerical simulation FFTW library
#include <stdio.h>              //for printing the help text
#include <math.h>               //for various math functions

class Quad
{
public:
    Quad(){
        glBegin(GL_QUADS);
    }
    void addPoint(float x, float y, float z){
        glVertex3f(x,y,z);
    }
    void draw(){
        glEnd();
    }

} quad;

#endif //SMOKE_QUAD_H
**/

#include <GL/gl.h>
#pragma once



class Quad										//A simple class for passing a quad polygon to OpenGL
{
public:

    Quad()									//Constructor: Starts drawing the quad
    {  glBegin(GL_QUADS); }
    ~Quad()
    {}

    void	addPoint(float x, float y, float z)		//Adds a new vertex to the quad
    {  glVertex3f(x,y,z); }

    void	addNormal(float* n)						//Adds a new normal to the quad
    {  glNormal3f(n[0],n[1],n[2]); }

    void    draw()									//Drawing function: ends the quad definition. This also draws the quad.
    { glEnd(); }
};




