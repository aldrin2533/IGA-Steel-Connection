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
void perform_PCGM(vector<double> Matrix_A, int Matrix_A_size,vector<double> load_vector,double TOL,int max_iter,vector<double>& unknowns,vector<double> initial_values)
{
    int iter = 1;

    // Initialize values of unknown
    //vector<double> unknowns;
        for(int i =0;i<Matrix_A_size;i++){unknowns[i]=initial_values[i];}

    vector<double> vect_r;
    double sum;
    vect_r.resize(Matrix_A_size,0);
        for(int i=0;i<Matrix_A_size;i++){
            sum=0;
        for(int j=0;j<Matrix_A_size;j++){
            sum+=(Matrix_A[Matrix_A_size*i+j]*unknowns[j]);
        }
            vect_r[i]=load_vector[i]-sum;
            //cout<<vect_r[i]<<endl;
        }

    vector<double> Matrix_A_diagonal;
    Matrix_A_diagonal.resize(Matrix_A_size,0);
    for(int i = 0;i<Matrix_A_size;i++){
            Matrix_A_diagonal[i]=Matrix_A[Matrix_A_size*i+i];
        }

    vector<double> vect_d;
    vect_d.resize(Matrix_A_size,0);
        for(int i=0;i<Matrix_A_size;i++){
        vect_d[i]=(vect_r[i]/Matrix_A_diagonal[i]);
        //cout<<fixed<<setprecision(8)<<vect_d[i]<<endl;
        }

    double c_new=0;
        for(int i=0;i<Matrix_A_size;i++){
            c_new+=vect_r[i]*vect_d[i];
        }
        //cout<<"cnew"<<c_new<<endl;

    double c_0=c_new;
    cout<<"c_0"<<c_0<<endl;
    double c_old=0;

    vector<double> vect_q;
    vect_q.resize(Matrix_A_size,0);
    double alpha=0,denom=0,beta=0;
    vector<double> vect_s;
    vect_s.resize(Matrix_A_size,0);

   do
    //while(c_new>pow(TOL,2)*c_0)
    {
        //cout<<"******************ITER :"<<iter<<endl;
        //vect_q.clear();
        double sum=0;
        for(int i=0;i<Matrix_A_size;i++){
        for(int j=0;j<Matrix_A_size;j++){
            sum+=Matrix_A[Matrix_A_size*i+j]*vect_d[j];
        }
        vect_q[i]=sum;
        //cout<<vect_q[i]<<endl;
        sum=0;
        }

        denom=0;
        for(int i=0;i<Matrix_A_size;i++){
            denom+=vect_d[i]*vect_q[i];
            }

        alpha=c_new/(denom);

        for(int i=0;i<Matrix_A_size;i++){
            unknowns[i]+=(alpha*vect_d[i]);
            vect_r[i]=vect_r[i]-alpha*vect_q[i];
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
    //cout<<fixed<<setprecision(16)<<c_new/c_0<<endl;
    }
while(iter<max_iter && c_new>pow(TOL,1)*c_0);
//cout<<c_new<<" "<<c_0<<endl;
printf("\t\t Performed Conjugate Gradient Method with %d Iterations, TOL: %.11f \n",iter,c_new/c_0);


}


//----------------------------------------------------------------------------------------------------------------------------------
//                                                                                                            ** verified 04/08/2021
//
//----------------------------------------------------------------------------------------------------------------------------------
