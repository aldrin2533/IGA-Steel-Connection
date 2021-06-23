#include<iostream>
#include<cmath>
#include<iomanip>
#include <fstream>
#include<string>
#include<vector>
#include<sstream>
#include<iterator>
#include <cstdio>
#include "data.h"
#include "functions.h"

using namespace std;

int main()
{
//----------------------------------------------------------------------------------------------------------------------------------
//                                          READ DATA FROM TEXT FILES AND STORE TO VECTORS
//----------------------------------------------------------------------------------------------------------------------------------
printf("Reading Data...\n");
int noctrlpt_actual=286;
    //  *************************************************
    //  READ KNOT VECTORS DATA
    //  Xi{} for Xi Direction
    //  Eta{} for Eta Direction
    //  Zeta{} for Zeta Direction
    //  *************************************************
    printf("\t Knot Vectors \n");
    vector<double> Xi,Eta;
    Xi.clear();
    Eta.clear();
    text_file DATA;
    DATA.input_string_file_name();
        int size_KV_xi = DATA.determine_size(DATA.KV_xi_string);
        int size_KV_eta = DATA.determine_size(DATA.KV_eta_string);
    for(int i=0;i<size_KV_xi;i++){
        Xi.push_back(DATA.read_textfile_double(DATA.KV_xi_string,i));}
    for(int i=0;i<size_KV_eta;i++){
        Eta.push_back(DATA.read_textfile_double(DATA.KV_eta_string,i));}

    //  *************************************************
    //  PROPERTIES OF B-SPLINE BASIS FUNCTIONS
    //  *************************************************
        printf("\t Knot Vectors Properties \n");
        calc_bspline_prop(Xi,p,n);
        calc_bspline_prop(Eta,q,m);

    //  *************************************************
    //  MODEL PROPERTIES
    //  *************************************************
        calc_nurbs_mod_prop(p,q,n,m,nel,nnp,nen,noctrlpt,noxi,noeta,nophyspts,Xi,Eta);

    //  *************************************************
    //  READ CONTROL POINTS DATA
    //  ctrlpt_coordinates_t[][4]     Control points coordinates.
    //  Memb_Inc[][27]              Member Incidences.
    //  Supp_cond[][]               Support conditions
    //                              (1 = Free to move,
    //                               0 = Restrained ).
    //  *************************************************
    printf("\t Control Points Data \n");

    vector<double> ctrlpt_coordinates_t;
    vector<double> ctrlpt_coordinates;
    ctrlpt_coordinates_t.clear();
        for(int j=0;j<noctrlpt_actual;j++){
        for(int i=0;i<4;i++){
            ctrlpt_coordinates_t.push_back(DATA.read_textfile_double(DATA.ctrl_pt_coords_string,4*j+i));
            ctrlpt_coordinates.push_back(DATA.read_textfile_double(DATA.ctrl_pt_coords_string,4*j+i));
            }}
    vector<int> ctrlpt_incidence;
    ctrlpt_incidence.clear();
        for(int j=0;j<nel;j++){
        for(int i=0;i<nen;i++){
            ctrlpt_incidence.push_back(DATA.read_textfile_int(DATA.ctrl_pt_inc_string,nen*j+i));
            }}

    printf("Successfully Read Data...\n");


//----------------------------------------------------------------------------------------------------------------------------------
//                                                  PRINT MODEL PROPERTIES
//----------------------------------------------------------------------------------------------------------------------------------

printf("NURBS Model Properties \n");
printf("\t Number of Elements: %d \n",nel);
printf("\t Number of Local Basis Functions with support to an element: %d \n",nen);
printf("\t Number of Control Points: %d \n",noctrlpt);
printf("\t Order of NURBS basis functions in Xi,Eta and Zeta direction: %d,%d \n",p,q);
printf("\t Number of Basis Functions in Xi, Eta and Zeta direction: %d,%d \n",n,m);


//----------------------------------------------------------------------------------------------------------------------------------
//                                                  PARAMETERS AT TIME t=0
//----------------------------------------------------------------------------------------------------------------------------------

    //  *************************************************
    //  LOAD AND INITIALIZE GLOBAL PARAMETERS AT TIME, t=0
    //  *************************************************


        load_global_param(u_global_t,u_g_t,Fext_global_t,Fext_g_t,K_global_t,DOF,noctrlpt_actual,no_unrestrained_dof);

        load_postproc_param(strain_t,noctrlpt_actual);

        load_D_e(D_e,bulk_mod,rank_4_unit,shear_mod,dev_tens);

        build_NURBS_coords(Xi,Eta,NURBS_coords);

        //for(int i=0;i<36;i++){cout<<D_e[i]<<endl;}


//----------------------------------------------------------------------------------------------------------------------------------
//                                                  GLOBAL STIFFNESS MATRIX AT TIME t
//----------------------------------------------------------------------------------------------------------------------------------

printf("Calculate Global Stiffness Matrix using converged parameters at time t...\n");
    //  *************************************************
    //  1. INITIALIZE ELEMENT AND KNOT POINT PARAMETERS
    //  *************************************************

        int element=0;
        xi=0;
        eta=0;
        xi_index=0;
        eta_index=0;

    //  ********************************************************
    //  2. LOOPING THROUGH ELEMENTS, KNOT VECTORS AND GAUSS POINTS
    //  TO CALCULATE CONTRIBUTION OF AN ELEMENT TO THE GLOBAL
    //  STIFFNESS MATRIX
    //  ********************************************************

        //  ********************************************************
        //  2.1 LOOP THROUGH ELEMENTS - (START)
        //  ********************************************************
            //printf("\t 1.1 Looping through elements \n");
            //for(int zeta_i=0;zeta_i<l-r;zeta_i++){
            //for(int eta_i=0;eta_i<m-q;eta_i++){
            //for(int xi_i=0;xi_i<n-p;xi_i++){
            vector<double> R_final;
            R_final.resize(nen,0);
            vector<double> Force_vector;
            Force_vector.resize(DOF*noctrlpt_actual,0);

            for(int el=0;el<nel;el++){
            element++;
            fill(R_final.begin(),R_final.end(),0);
            //printf("\t \t Element %d \n",element);

                //  ********************************************************
                //  2.1.1 DETERMINE NURBS COORDINATES
                //  ********************************************************
                    //printf("\t\t 1.2 Looping through Knot Vectors \n");
                    //xi_index=xi_i+p;
                    //eta_index=eta_i+q;
                    //zeta_index=zeta_i+r;
                    xi_index=NURBS_coords[2*(element-1)];
                    eta_index=NURBS_coords[2*(element-1)+1];
                    //cout<<xi_index<<" "<<eta_index<<endl;
        //  ********************************************************
        //  2.2 LOOP THROUGH GAUSS POINTS - (START)
        //  ********************************************************
            //printf("\t\t 1.2 Looping through Gauss Points \n");
            int gauss_pt_counter=0;
            //for(int gp_z=0;gp_z<NGP;gp_z++){
            for(int gp_y=0;gp_y<NGP;gp_y++){
            for(int gp_x=0;gp_x<NGP;gp_x++){
            gauss_pt_counter++;

                //  ********************************************************
                //  2.2.1 DETERMINE PARAMETRIC COORDINATE
                //  ********************************************************
                    xi=calc_parametric_coord(Xi,xi_index,gauss_pts[gp_x]);
                    eta=calc_parametric_coord(Eta,eta_index,gauss_pts[gp_y]);
                    //printf("*** %f %f \n",xi,eta);
        //  ********************************************************
        //  2.3 CALCULATE B-SPLINE BASIS FUNCTIONS AND THEIR DERIVATIVES
        //  ********************************************************
            form_B_splines(N,N_der,xi_index,p,xi,Xi,n);
            //printf("1,2,3 %f %f %f \n",N[0],N[1],N[2]);
            form_B_splines(M,M_der,eta_index,q,eta,Eta,m);

            //cout<<N_der[0]<<" "<<N_der[1]<<" "<<N_der[2]<<endl;
            //cout<<M_der[0]<<" "<<M_der[1]<<endl;

            //cout<<L[0]<<" "<<L[1]<<" "<<L[2]<<endl;
        //  ********************************************************
        //  2.4 DETERMINE NURBS BASIS FUNCTIONS AND THEIR DERIVATIVES
        //      w.r.t PHYSICAL COORDINATES
        //  ********************************************************
            Shape_func_1(ctrlpt_incidence,ctrlpt_coordinates_t,R,dR_dxi,N,M,N_der,M_der,p,q,element,nen,dR_dx,dx_dxi);
            //cout<<N_der[0]<<" "<<N_der[1]<<" "<<N_der[2]<<endl;
            //for(int i=0;i<nen;i++){printf("%f %f \n",dR_dxi[2*i],dR_dxi[2*i+1]);}



        //  ********************************************************
        //  2.5 CALCULATE DETERMINANT OF THE JACOBIAN MATRIX
        //      (Mapping from parametric to physical coordinates)
        //  ********************************************************
            J=Shape_function_2(xi_index,eta_index,Xi,Eta,dx_dxi);
            //cout<<"Determinant:"<<J<<endl;
            //cout<<"endl J"<<endl;

            for(int i=0;i<nen;i++){
            R_final[i]+=(J*R[i]*gauss_pts_wts[gp_x]*gauss_pts_wts[gp_y]);
            //cout<<"R"<<i<<" "<<J*R[i]*gauss_pts_wts[gp_x]*gauss_pts_wts[gp_y]<<endl;
            }


        //  ********************************************************
        //  2.7 LOOP THROUGH GAUSS POINTS - (END)
        //  ********************************************************
            }}
            //for(int i=0;i<nen;i++){cout<<R_final[i]<<endl;}

        //  ********************************************************
        //  2.8 LOOP THROUGH ELEMENTS - (END)
        //  ********************************************************
             for(int i=0;i<nen;i++){
                int cc=ctrlpt_incidence[(element-1)*nen+i];
                //cout<<cc<<endl;
                Force_vector[3*(cc-1)]+=R_final[i];
                //Force_vector[3*(cc-1)+1]=R_final[i];
                //Force_vector[3*(cc-1)+2]+=R_final[i];
                }
                //cout<<R_final[i]<<endl;}

            }
            //}}}
            for(int i=0;i<DOF*noctrlpt_actual;i++){cout<<Force_vector[i]<<endl;}


    return 0;
}
