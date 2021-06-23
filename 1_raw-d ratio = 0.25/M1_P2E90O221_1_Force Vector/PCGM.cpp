#include<string>
#include<iostream>
#include<vector>
#include<fstream>
#include<cmath>
#include<iomanip>

#include "functions.h"

using namespace std;

//----------------------------------------------------------------------------------------------------------------------------------
//                                             CONJUGATE GRADIENT METHOD W/ PRE-CONDITIONER
//----------------------------------------------------------------------------------------------------------------------------------
void perform_PCGM(vector<double> Matrix_A, int Matrix_A_size,vector<double> load_vector,double TOL,int max_iter,vector<double>& unknowns)
{
    int iter = 1;

    // Initialize values of unknown
    //vector<double> unknowns;
        for(int i =0;i<Matrix_A_size;i++){unknowns.push_back(0);}

    vector<double> vect_r;
    double sum;
    vect_r.resize(Matrix_A_size);
        for(int i=0;i<Matrix_A_size;i++){
            sum=0;
        for(int j=0;j<Matrix_A_size;j++){
            sum+=(Matrix_A[Matrix_A_size*i+j]*unknowns[j]);
        }
            vect_r[i]=load_vector[i]-sum;
            //cout<<vect_r[i]<<endl;
        }

    vector<double> Matrix_A_diagonal;
    for(int i = 0;i<Matrix_A_size;i++){
            Matrix_A_diagonal.push_back(Matrix_A[Matrix_A_size*i+i]);
            //cout<<fixed<<setprecision(8)<<Matrix_A_diagonal[i]<<endl;
        }

    vector<double> vect_d;
    for(int i=0;i<Matrix_A_size;i++){
        vect_d.push_back(vect_r[i]/Matrix_A_diagonal[i]);
        //cout<<fixed<<setprecision(8)<<vect_d[i]<<endl;
        }

    double c_new=0;
        for(int i=0;i<Matrix_A_size;i++){
            c_new+=vect_r[i]*vect_d[i];
        }

    double c_0=c_new;
    double c_old=0;

    vector<double> vect_q;
    vect_q.resize(Matrix_A_size);
    double alpha=0,denom=0,beta=0;
    vector<double> vect_s;
    vect_s.resize(Matrix_A_size);

    while(iter<max_iter && c_new>pow(TOL,2)*c_0)
    //while(c_new>pow(TOL,2)*c_0)
    {
        vect_q.clear();
        double sum=0;
        for(int i=0;i<Matrix_A_size;i++){
        for(int j=0;j<Matrix_A_size;j++){
            sum+=Matrix_A[Matrix_A_size*i+j]*vect_d[j];
        }
        vect_q.push_back(sum);
        sum=0;
        }

        denom=0;
        for(int i=0;i<Matrix_A_size;i++){
            denom+=vect_d[i]*vect_q[i];
            }

        alpha=c_new/(denom);

        for(int i=0;i<Matrix_A_size;i++){
            unknowns[i]+=alpha*vect_d[i];
            vect_r[i]-=alpha*vect_q[i];
            vect_s[i]=vect_r[i]/Matrix_A_diagonal[i];
            //cout<<fixed<<setprecision(8)<<unknowns[i]<<endl;
            //cout<<fixed<<setprecision(8)<<vect_r[i]<<endl;
            }

        c_old=c_new;

        c_new=0;
        for(int i=0;i<Matrix_A_size;i++){
            c_new+=vect_r[i]*vect_s[i];
            }
            //cout<<c_new;

        beta=c_new/c_old;

        for(int i=0;i<Matrix_A_size;i++){
            vect_d[i]=vect_s[i]+beta*vect_d[i];
            //cout<<fixed<<setprecision(8)<<vect_d[i]<<endl;
            }

        iter++;


    //printf("Iterations:%d \n",iter);
    }

//printf("\t\t Performed Conjugate Gradient Method with %d Iterations \n",iter);


}


//----------------------------------------------------------------------------------------------------------------------------------
//                                                                                                            ** verified 04/08/2021
//
//----------------------------------------------------------------------------------------------------------------------------------








