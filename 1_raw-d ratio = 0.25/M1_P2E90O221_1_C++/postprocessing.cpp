#include<string>
#include<iostream>
#include<vector>
#include<fstream>
#include<cmath>
#include<iomanip>

#include "functions.h"

using namespace std;

//----------------------------------------------------------------------------------------------------------------------------------
//                                                    CALCULATE PHYSICAL POINT
//----------------------------------------------------------------------------------------------------------------------------------
void calc_phy_pt(vector<int> ctrlpt_incidence,vector<double> ctrlpt_coord,vector<double> R,int element,int nen,vector<double>& physical_point)
{
    fill(physical_point.begin(),physical_point.end(),0);
    for(int i=0;i<nen;i++)
    {
        int cc=ctrlpt_incidence[nen*(element-1)+i];
        //cout<<cc<<" "<<u_global_t[3*(cc-1)+2]<<endl;
        physical_point[0]+=R[i]*ctrlpt_coord[4*(cc-1)];
        physical_point[1]+=R[i]*ctrlpt_coord[4*(cc-1)+1];
        physical_point[2]+=R[i]*ctrlpt_coord[4*(cc-1)+2];
    }
}


//----------------------------------------------------------------------------------------------------------------------------------
//                                                    CALCULATE DEFORMATION AT A POINT
//----------------------------------------------------------------------------------------------------------------------------------
void calc_displacement(vector<int> ctrlpt_incidence,vector<double> u_global_t,vector<double> R,int element,int nen,vector<double>& disp_point)
{
    fill(disp_point.begin(),disp_point.end(),0);
    for(int i=0;i<nen;i++)
    {
        int cc=ctrlpt_incidence[nen*(element-1)+i];
        //cout<<cc<<" "<<u_global_t[3*(cc-1)+2]<<endl;
        disp_point[0]+=R[i]*u_global_t[3*(cc-1)];
        disp_point[1]+=R[i]*u_global_t[3*(cc-1)+1];
        disp_point[2]+=R[i]*u_global_t[3*(cc-1)+2];
    }
}

//----------------------------------------------------------------------------------------------------------------------------------
//                                                    CALCULATE STRAIN AT A POINT
//----------------------------------------------------------------------------------------------------------------------------------
void calc_strain(vector<double>& du_dx_point, vector<double> dR_dx,vector<double> u_global_t,int deg_of_freedom, int nen, vector<int> ctrlpt_incidence,int element, vector<double>& strain_point)
{
    fill(du_dx_point.begin(),du_dx_point.end(),0);
    for(int i=0;i<deg_of_freedom;i++){
        for(int j=0;j<deg_of_freedom;j++){
            for(int aa=0;aa<nen;aa++){
                int cc=ctrlpt_incidence[nen*(element-1)+aa];
                //cout<<u_global_t[deg_of_freedom*(cc-1)+i]<<endl;
                du_dx_point[deg_of_freedom*i+j]+=(dR_dx[deg_of_freedom*aa+j]*u_global_t[deg_of_freedom*(cc-1)+i]);
                //du_dx_point[deg_of_freedom*i+j]+=(dR_dx[deg_of_freedom*aa+j]);
        }}}

    fill(strain_point.begin(),strain_point.end(),0);
    strain_point[0]=du_dx_point[0]; //  dux_dx
    strain_point[1]=du_dx_point[4]; //  duy_dy
    strain_point[2]=du_dx_point[8]; //  duz_dz
    // Engineering strains
    strain_point[3]=(du_dx_point[5]+du_dx_point[7]);    //  duy_dz & duz_dy
    strain_point[4]=(du_dx_point[6]+du_dx_point[2]);    //  dux_dz & duz_dx
    strain_point[5]=(du_dx_point[1]+du_dx_point[3]);    //  dux_dy & duy_dx

}


//----------------------------------------------------------------------------------------------------------------------------------
//
//----------------------------------------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------------------------------------
//                                                    CALCULATE STRESS AT A POINT
//----------------------------------------------------------------------------------------------------------------------------------
void calc_stress(vector<double> strain_point,vector<double>& stress_point,vector<double> D_e)
{
    fill(stress_point.begin(),stress_point.end(),0);
    for(int i=0;i<6;i++){
        for(int j=0;j<6;j++){
            stress_point[i]+=(D_e[6*i+j]*strain_point[j]);
        }
    }

}

//----------------------------------------------------------------------------------------------------------------------------------
//                                                    EQUIVALENT STRESS/STRAIN
//----------------------------------------------------------------------------------------------------------------------------------

double calc_equiv(vector<double> stress_or_strain)
{
    double value=0;
    vector<int> operator_1{1,1,1,2,2,2};
    for(int i=0;i<6;i++)
    {
        value+=operator_1[i]*pow(stress_or_strain[i],2);
    }
    value=sqrt(value);

    return value;
}


