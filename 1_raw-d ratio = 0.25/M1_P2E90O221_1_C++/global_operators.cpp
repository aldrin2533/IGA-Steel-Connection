#include<string>
#include<iostream>
#include<vector>
#include<fstream>
#include<cmath>

#include "functions.h"


using namespace std;


int calc_num_lines_textfile(string file_name)
{
    string line;
    int num_lines=0;
    std::ifstream ifile(file_name, std::ios::in);
    if (ifile.is_open())
    {
        while(!ifile.eof())
            {
            getline(ifile,line);
            num_lines++;
            }
            ifile.close();
    }
    return num_lines-1;
}


//----------------------------------------------------------------------------------------------------------------------------------
//                                                      LOAD ELASTIC TANGENT MODULUS
//----------------------------------------------------------------------------------------------------------------------------------

void load_D_e(vector<double>& D_e,double bulk_mod,vector<int> rank_4_unit,double shear_mod,vector<double> dev_tens)
    {
       for(int i=0;i<36;i++){
            D_e.push_back(0);}

        for(int i=0;i<6;i++){
        for(int j=0;j<6;j++){
            D_e[6*i+j]=bulk_mod*rank_4_unit[6*i+j]+2*shear_mod*dev_tens[6*i+j];
        }}
    }

//----------------------------------------------------------------------------------------------------------------------------------
//                                                          BUILD IEN ARRAY
//----------------------------------------------------------------------------------------------------------------------------------

void build_IEN(int nen,int nel,int n,int m,int l,int p,int q,int r,vector<int>& IEN)
{


    //double INN[nnp][2];                                 // NURBS coordinate system
    //for(int aa=0;aa<nnp;aa++)                           // Initialize Values of INN
    //{for(int bb=0;bb<2;bb++)
        //{INN[aa][bb]=0;}}

    //double IEN[nen][nel];                               // Connectivity array
    for(int aa=0;aa<nen*nel;aa++)                           // Initialize Values of IEN
        {IEN.push_back(0);}

    int A=0, B=0,el=0,b=0;                              // Initialize local variables

    for(int k=1;k<l+1;k++)
    {for(int j=1;j<m+1;j++)
    {for(int i=1;i<n+1;i++)
        {A=A+1;                                         // Global Function Number
        //INN[A-1][0]=i;                                  // Assigning NURBS coordinates
        //INN[A-1][1]=j;

            if(i>=(p+1)&&j>=(q+1)&&k>=(r+1))
            {el=el+1;
                for(int kloc=0;kloc<r+1;kloc++)
                {for(int jloc=0;jloc<q+1;jloc++)
                {for(int iloc=0;iloc<p+1;iloc++)
                    {
                        B=A-kloc*n*m-jloc*n-iloc;
                        b=kloc*(p+1)*(q+1)+jloc*(p+1)+iloc+1;            // Local Function Number
                        IEN[nen*(b-1)+(el-1)]=B;               // Assigning Connectivity Array
                    }
                }}}}}}

}


//----------------------------------------------------------------------------------------------------------------------------------
//                                                  LOAD AND INITIALIZE GLOBAL PARAMETERS
//----------------------------------------------------------------------------------------------------------------------------------
void load_global_param(vector<double>& u_global_t,vector<double>& u_g_t1,vector<double>& Fext_global_t,vector<double>& Fext_g_t,
                       vector<double>& K_global_t,int DOF,int noctrlpt,int no_unrestrained_dof,vector<double>& u_global_t1,vector<double>& u_g_t)
{
    u_global_t.resize(DOF*noctrlpt,0);
    u_global_t1.resize(DOF*noctrlpt,0);
    u_g_t.resize(no_unrestrained_dof,0);
    u_g_t1.resize(no_unrestrained_dof,0);
    Fext_global_t.resize(DOF*noctrlpt,0);
    Fext_g_t.resize(no_unrestrained_dof,0);
    K_global_t.resize(pow(DOF*noctrlpt,2),0);
}


//----------------------------------------------------------------------------------------------------------------------------------
//                                                     INITIALIZE LOCAL PARAMETERS
//----------------------------------------------------------------------------------------------------------------------------------

void initialize_knot_point_par(double& xi,double& eta, double& zeta,int& xi_index,int& eta_index,int& zeta_index)
{
    xi=0;eta=0;zeta=0;
    xi_index=0;eta_index=0;zeta_index=0;
}


//----------------------------------------------------------------------------------------------------------------------------------
//                                              BUILD LOCAL STIFFNESS MATRIX AT A CONTROL POINT
//----------------------------------------------------------------------------------------------------------------------------------

//  ********************************************************
//  COMPUTE CONTRIBUTION OF A GAUSS POINT TO THE GLOBAL
//  STIFFNESS MATRIX
//  ********************************************************

void build_k_local(double bulk_modulus,double shear_modulus,vector<double> dR_dx,double J_mod,int deg_of_freedom,vector<double>& K_local_mod,
                   int index_i,int index_j)
{
     //diagonals
     K_local_mod[0]=((bulk_modulus+4*shear_modulus/3)*dR_dx[deg_of_freedom*index_i+0]*dR_dx[deg_of_freedom*index_j+0]+
            shear_modulus*(dR_dx[deg_of_freedom*index_i+1]*dR_dx[deg_of_freedom*index_j+1]+dR_dx[deg_of_freedom*index_i+2]*dR_dx[deg_of_freedom*index_j+2]))*J_mod;
     K_local_mod[4]=((bulk_modulus+4*shear_modulus/3)*dR_dx[deg_of_freedom*index_i+1]*dR_dx[deg_of_freedom*index_j+1]+
            shear_modulus*(dR_dx[deg_of_freedom*index_i+0]*dR_dx[deg_of_freedom*index_j+0]+dR_dx[deg_of_freedom*index_i+2]*dR_dx[deg_of_freedom*index_j+2]))*J_mod;
     K_local_mod[8]=((bulk_modulus+4*shear_modulus/3)*dR_dx[deg_of_freedom*index_i+2]*dR_dx[deg_of_freedom*index_j+2]+
            shear_modulus*(dR_dx[deg_of_freedom*index_i+0]*dR_dx[deg_of_freedom*index_j+0]+dR_dx[deg_of_freedom*index_i+1]*dR_dx[deg_of_freedom*index_j+1]))*J_mod;

     //off - diagonal entries
     K_local_mod[1]=((bulk_modulus-2*shear_modulus/3)*dR_dx[deg_of_freedom*index_i+0]*dR_dx[deg_of_freedom*index_j+1]+
            shear_modulus*dR_dx[deg_of_freedom*index_i+1]*dR_dx[deg_of_freedom*index_j+0])*J_mod;
     K_local_mod[2]=((bulk_modulus-2*shear_modulus/3)*dR_dx[deg_of_freedom*index_i+0]*dR_dx[deg_of_freedom*index_j+2]+
            shear_modulus*dR_dx[deg_of_freedom*index_i+2]*dR_dx[deg_of_freedom*index_j+0])*J_mod;

     K_local_mod[3]=((bulk_modulus-2*shear_modulus/3)*dR_dx[deg_of_freedom*index_i+1]*dR_dx[deg_of_freedom*index_j+0]+
            shear_modulus*dR_dx[deg_of_freedom*index_i+0]*dR_dx[deg_of_freedom*index_j+1])*J_mod;
     K_local_mod[5]=((bulk_modulus-2*shear_modulus/3)*dR_dx[deg_of_freedom*index_i+1]*dR_dx[deg_of_freedom*index_j+2]+
            shear_modulus*dR_dx[deg_of_freedom*index_i+2]*dR_dx[deg_of_freedom*index_j+1])*J_mod;

     K_local_mod[6]=((bulk_modulus-2*shear_modulus/3)*dR_dx[deg_of_freedom*index_i+2]*dR_dx[deg_of_freedom*index_j+0]+
            shear_modulus*dR_dx[deg_of_freedom*index_i+0]*dR_dx[deg_of_freedom*index_j+2])*J_mod;
     K_local_mod[7]=((bulk_modulus-2*shear_modulus/3)*dR_dx[deg_of_freedom*index_i+2]*dR_dx[deg_of_freedom*index_j+1]+
            shear_modulus*dR_dx[deg_of_freedom*index_i+1]*dR_dx[deg_of_freedom*index_j+2])*J_mod;
}

//----------------------------------------------------------------------------------------------------------------------------------
//                                                      COMPUTE DETERMINANT OF A MATRIX
//----------------------------------------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------------------------------------
//                                         EXTRACT K_gg_t AND delta_F_g_t FROM GLOBAL VECTORS
//----------------------------------------------------------------------------------------------------------------------------------
void extr_reduced_vector(vector<double>& K_gg_t,vector<double>& delta_F_g_t1,vector<double> P_loads,vector<int> bound_condition,
                         vector<double> K_global_t,int DOF, int noctrlpt, int no_unrestrained_dof)
{
K_gg_t.resize(no_unrestrained_dof*no_unrestrained_dof,0);
    delta_F_g_t1.clear();
    int counter=0;
    for(int j=0;j<noctrlpt;j++){
    for(int jj=0;jj<DOF;jj++){
        if(bound_condition[DOF*j+jj]==0){
            for(int i=0;i<noctrlpt;i++){
            for(int ii=0;ii<DOF;ii++){
                if(bound_condition[DOF*i+ii]==0){
                    //K_gg_t.push_back(K_global_t[noctrlpt*DOF*(DOF*j+jj)+DOF*i+ii]);
                    K_gg_t[counter]=K_global_t[noctrlpt*DOF*(DOF*j+jj)+DOF*i+ii];
                    counter++;


                    //if(j==367 && jj==2 && i==367 && ii==2){

                    //printf("j: %d jj:%d i:%d ii:%d K:%f counter:%d \n",j,jj,i,ii,K_gg_t[1149183],counter);}
                }
            }}
        delta_F_g_t1.push_back(P_loads[DOF*j+jj]);
        }
    }}
}

//----------------------------------------------------------------------------------------------------------------------------------
//                                          EXTRACT K_gg_t ONLY, FROM GLOBAL VECTORS
//----------------------------------------------------------------------------------------------------------------------------------
void extr_reduced_vector_Kgg(vector<double>& K_gg_t,vector<int> bound_condition,
                         vector<double> K_global_t,int DOF, int noctrlpt)
{
    for(int j=0;j<noctrlpt;j++){
    for(int jj=0;jj<DOF;jj++){
        if(bound_condition[DOF*j+jj]==0){
            for(int i=0;i<noctrlpt;i++){
            for(int ii=0;ii<DOF;ii++){
                if(bound_condition[DOF*i+ii]==0){
                    K_gg_t.push_back(K_global_t[noctrlpt*DOF*(DOF*j+jj)+DOF*i+ii]);
                }
            }}
        }
    }}
}

//----------------------------------------------------------------------------------------------------------------------------------
//                                      CALCULATE INITIAL RESIDUAL FOR NEWTON RAPHSON METHOD
//----------------------------------------------------------------------------------------------------------------------------------

void calc_residual(vector<double> K_gg_t,vector<double>& r_initial,vector<double> Fext_g_t,vector<double> u_g_t,
                           double& norm_initial,int no_unrestrained_dof)
{
    norm_initial=0;
    r_initial.clear();
    r_initial.resize(no_unrestrained_dof);
        for(int i = 0;i<no_unrestrained_dof;i++){
        for(int j=0;j<no_unrestrained_dof;j++){
            r_initial[i]+=(-K_gg_t[no_unrestrained_dof*i+j]*u_g_t[j]);
        }
            r_initial[i]+=Fext_g_t[i];
            norm_initial+=pow(r_initial[i],2);
        }
    norm_initial=sqrt(norm_initial);
}

//----------------------------------------------------------------------------------------------------------------------------------
//                                          SET CURRENT RESIDUAL = TO INCREMENTAL LOAD
//----------------------------------------------------------------------------------------------------------------------------------
void set_curr_residual(vector<double>& r_current,int no_unrestrained_dof,vector<double> delta_F_g_t1,double& norm_current)
{
    norm_current=0;
    r_current.clear();
    r_current.resize(no_unrestrained_dof);
        for(int i=0;i<no_unrestrained_dof;i++){
        r_current[i]=delta_F_g_t1[i];
        }
}

//----------------------------------------------------------------------------------------------------------------------------------
//                                                      UPDATE CONTROL POINTS
//----------------------------------------------------------------------------------------------------------------------------------
void update_ctrl_pt(vector<double>& ctrlpt_coordinates,vector<double> delta_u_g_t1,int no_unrestrained_dof,vector<int> unrestrained_dof_data,double alpha,vector<double> ctrlpt_coordinates_prev)
{
    for(int i=0;i<no_unrestrained_dof;i++){
        int subscript=unrestrained_dof_data[2*i]-1;
        int subscript_dof=unrestrained_dof_data[2*i+1];
        ctrlpt_coordinates[4*subscript+subscript_dof]=ctrlpt_coordinates_prev[4*subscript+subscript_dof]+alpha*delta_u_g_t1[i];
    }
}

//----------------------------------------------------------------------------------------------------------------------------------
//                                                  LOAD AND INITIALIZE POST PROCESSING PARAMETERS
//----------------------------------------------------------------------------------------------------------------------------------
void load_postproc_param(vector<double>& strain_t,int noctrlpt)
{
    strain_t.resize(6*noctrlpt,0);
}

//----------------------------------------------------------------------------------------------------------------------------------
//                                                          CALCULATE DEVIATORIC STRESS
//----------------------------------------------------------------------------------------------------------------------------------

void calc_dev_stress(vector<double> stress_point, vector<double>& dev_stress_point)
{
    fill(dev_stress_point.begin(),dev_stress_point.end(),0);
    double ave=0;
    double sum=0;

    for(int i=0;i<3;i++){
        sum+=stress_point[i];
    }

    ave=sum/3;

    for(int i=0;i<6;i++){
        if(i<3){dev_stress_point[i]=stress_point[i]-ave;}
        else{dev_stress_point[i]=stress_point[i];}
    }
}

//----------------------------------------------------------------------------------------------------------------------------------
//                                                          EVALUATE YIELD CRITERION
//----------------------------------------------------------------------------------------------------------------------------------

void evaluate_yield_func(vector<double> dev_stress_point,double fy)
{
    double norm=sqrt(pow(dev_stress_point[0],2)+pow(dev_stress_point[1],2)+pow(dev_stress_point[2],2)+
                2*pow(dev_stress_point[3],2)+2*pow(dev_stress_point[4],2)+2*pow(dev_stress_point[5],2));
    double diff = norm-fy;

    //if(diff<0){printf("Elastic Gauss Point, %f \n",norm);}
    //else{printf("Plastic Gauss Point \n");}
}



