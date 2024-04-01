// Copyright (c) 2014  
// University of Illinois at Urbana-Champaign, All rights reserved. 
//
// This file is a CGAL version to calculate the smallest enclosing ellipse 
// using CGAL library (https://www.cgal.org/). The function prototype is:
//     [major,minor,angle,x,y]=CGAL_EllipseSize(c);     
//     Input:   
//       c - the binary matrix of the single particle image
//     Output:  
//       major - the lenght of major axis of the ellipse
//       minor - the lenght of minor axis of the ellipse
//       angle - the angle of the ellipse in radian
//       x     - the x position of the ellipse center
//       y     - the y position of the ellipse center
//
// Author(s)     : Wei Wu   weiwu3@illinois.edu
//               : Lee Dunnavan (edwin.dunnavan@noaa.gov)
//               : Joe Finlon (joseph.a.finlon@nasa.gov)
//
// Compile command: mex CGAL_EllipseSize.cpp -I/opt/local/include -L/opt/local/lib -lCGAL -lgmp
// 
// Code History:
//   * First version created by Wei Wu, July 24, 2014
//   * Added center x/y position, Lee Dunnavan and Joe Finlon, 29 Jan 2024
//
#include "mex.h"
#include <iostream>
#include <cassert>
#include <stdio.h>
#include <stdlib.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Min_ellipse_2.h>
#include <CGAL/Min_ellipse_2_traits_2.h>
#include <CGAL/Gmpq.h>
typedef  CGAL::Gmpq                       NT;
typedef  CGAL::Cartesian<NT>              K;
typedef  CGAL::Point_2<K>                 Point;
typedef  CGAL::Min_ellipse_2_traits_2<K>  Traits;
typedef  CGAL::Min_ellipse_2<Traits>      Min_ellipse;
#define PI 3.1415926

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{   
      
    #define B_OUT1 plhs[0]
    #define B_OUT2 plhs[1]
    #define B_OUT3 plhs[2]
    #define B_OUT4 plhs[3] // center x
    #define B_OUT5 plhs[4] // center y

    #define A_IN prhs[0]

    double a,b,c,d,f,g;
    double *B0,*B1,*B2,*B3,*B4;
    int m,n;
    mxChar *ICEPIC;
    

    B_OUT1=mxCreateDoubleScalar(0.0);
    B_OUT2=mxCreateDoubleScalar(0.0);
    B_OUT3=mxCreateDoubleScalar(0.0);
    B_OUT4=mxCreateDoubleScalar(0.0);
    B_OUT5=mxCreateDoubleScalar(0.0);

    B0=mxGetPr(B_OUT1);
    B1=mxGetPr(B_OUT2);
    B2=mxGetPr(B_OUT3);
    B3=mxGetPr(B_OUT4);
    B4=mxGetPr(B_OUT5);

    m=mxGetM(A_IN);
    n=mxGetN(A_IN);    
    ICEPIC=mxGetChars(A_IN);

    int i,j,k=0,mm=0,nn=0;
    //double L,W,anglex1,angley1,anglex2,angley2;
    std::vector<Point> points;
    for( i = 0; i < n; i++ )
    {
        for (j=0; j<m; j++)
        {
            //cout<< ICEPIC[i*m+j] << ", ";
            if(ICEPIC[i*m+j]==48) //48
            {
                points.push_back(Point(i,j));
            }
        }
        //cout<<"\n"<<endl;
        
    }
    //cout<< k << endl;

    if (points.size()==0)
    {
        //cout << "No Illuminated Doide!\n" << endl;
        *B0 = 0;
        *B1 = 0;
        *B2 = 0;
        *B3 = 0;
        *B4 = 0;
        return;
    }
    
    Min_ellipse  me2( points.begin(), points.end(), true);     // fast
    me2.ellipse().double_coefficients( a, c, b, d, f, g);
    // Weird fix needed in case minimum_ellipse_2 function gives negative parameter for a (and others).
    if (a<0)
    {
       a = -a;
       c = -c;
       b = -b;
       d = -d;
       f = -f;
       g = -g;

    }
//     std::cout << "ellipse has the equation " <<
//       r << " x^2 + " <<
//       s << " y^2 + " <<
//       t << " xy + " <<
//       u << " x + " <<
//       v << " y + " <<
//       w << " = 0." << std::endl;
    
    double major,minor;
    major=2*sqrt((2*(a*pow(f/2,2)+c*pow(d/2,2)+g*pow(b/2,2)-b*d*f/4-a*c*g))/((pow(b/2,2)-a*c)*(sqrt(pow((a-c),2)+pow(b,2))-(a+c))));
    minor=2*sqrt((2*(a*pow(f/2,2)+c*pow(d/2,2)+g*pow(b/2,2)-b*d*f/4-a*c*g))/((pow(b/2,2)-a*c)*(-sqrt(pow((a-c),2)+pow(b,2))-(a+c))));
    
    *B0=major>minor? major:minor;
    *B1=major>minor? minor:major;
    
    if (0==b && a<c)
        *B2=0; 
    else if (0==b && a>c)
        *B2=PI/2;
    else if (0 != b && a<c)
        *B2=(atan(b/(a-c)))/2;
    else if (0 != b && a>c)
        *B2=(atan(b/(a-c)))/2+PI/2;
    else if (a==c && b>0) // this line and rest of block for refined angles
        *B2 = -PI/4;
    else if (a==c && b<0)
        *B2 = PI/4;
    else
        *B2=0;
    
    if (*B2 > PI/2)
        *B2 = *B2 - PI;
    
    if (*B2 < -PI/2)
        *B2 = *B2 + PI;
    
    *B3 = 1+(2*c*d-b*f)/(b*b-4*a*c);
    *B4 = 1+(2*a*f-b*d)/(b*b-4*a*c);
    
    return;   
}
