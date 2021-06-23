#include<string>
#include<iostream>
#include<vector>
#include<fstream>
#include<cmath>
#include<iomanip>

#include "functions.h"

using namespace std;

//----------------------------------------------------------------------------------------------------------------------------------
//                                                                  GAUSS-SEIDEL METHOD
//----------------------------------------------------------------------------------------------------------------------------------
void perform_gauss_seidel(vector<double> Matrix_A, int Matrix_A_size,vector<double> load_vector,double TOL,int max_iter,vector<double>& unknowns)
{
    fill(unknowns.begin(),unknowns.end(),0);
    double sum_1=0;
    double sum_2=0;

    vector<double> unknowns_old;
    unknowns_old.resize(Matrix_A_size,0);
    int iter=0;

    while(iter<10)
    {
        iter++;
        cout<<"iter"<<iter<<endl;

        for(int i=0;i<Matrix_A_size;i++)
            {
            sum_1=0;
            sum_2=0;

                if(i==0){
                    sum_1=0;}

                else{
                    for(int j=0;j<i;j++)
                        {
                        //if(i==0){sum_1=0;}
                        //else if(i==1){sum_1+=Matrix_A[Matrix_A_size*i+j]*unknowns[j]/Matrix_A[Matrix_A_size*i+i];}
                        sum_1+=(Matrix_A[Matrix_A_size*i+j]*unknowns[j]/Matrix_A[Matrix_A_size*i+i]);
                                        //cout<<"i j"<<i<<" "<<j<<endl;
                                        //cout<<Matrix_A[Matrix_A_size*i+j]<<endl;
                                        //cout<<unknowns[j]<<endl;}
                    }}

                        cout<<"sum1: "<<sum_1<<endl;

                if(i==(Matrix_A_size-1)){
                    sum_2=0;}

                else{

                    for(int j=i+1;j<Matrix_A_size;j++)
                        {
                            sum_2+=Matrix_A[Matrix_A_size*i+j]*unknowns_old[j]/Matrix_A[Matrix_A_size*i+i];
                    }}

                        cout<<"sum2: "<<sum_2<<endl;

            unknowns[i]=load_vector[i]/Matrix_A[Matrix_A_size*i+i]-sum_1-sum_2;

            }

        for(int i=0;i<Matrix_A_size;i++)
            {
                cout<<unknowns[i]<<endl;
                unknowns_old[i]=unknowns[i];
            }
    }

        //for(int i=0;i<Matrix_A_size;i++)
          //  {
            //    cout<<unknowns[i]<<endl;
            //}


}


//----------------------------------------------------------------------------------------------------------------------------------
//                                                                                                            ** verified 04/08/2021
//
//----------------------------------------------------------------------------------------------------------------------------------
