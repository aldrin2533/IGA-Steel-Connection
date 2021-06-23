#include<string>
#include<iostream>
#include<vector>
#include<fstream>
#include<cmath>
#include<iomanip>

#include "functions.h"

using namespace std;

//----------------------------------------------------------------------------------------------------------------------------------
//                                                          GAUSS - JORDAN ELIMINATION
//----------------------------------------------------------------------------------------------------------------------------------
void perform_gauss_jordan(vector<double> Matrix_A, int n,vector<double> load_vector,vector<double>& unknowns)
{
    //n Size of a single column
    fill(unknowns.begin(),unknowns.end(),0);

    double b=0;

    for(int j=0; j<n; j++) {
      for(int i=0; i<n; i++) {
         if(i!=j) {
            b=Matrix_A[n*i+j]/Matrix_A[n*j+j];
            for(int k=0; k<=n; k++) {
               if(k<n){
               Matrix_A[n*i+k]=Matrix_A[n*i+k]-b*Matrix_A[n*j+k];
               }
               else if(k==n){
               load_vector[i]=load_vector[i]-b*load_vector[j];
               }
            }
         }
      }
   }

   for(int i=0; i<n; i++) {
      unknowns[i]=load_vector[i]/Matrix_A[n*i+i];
   }



}


//----------------------------------------------------------------------------------------------------------------------------------
//                                                                                                            ** verified 04/08/2021
//
//----------------------------------------------------------------------------------------------------------------------------------
