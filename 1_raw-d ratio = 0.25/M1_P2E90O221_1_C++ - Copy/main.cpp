#include<iostream>
#include<cmath>
#include<iomanip>
#include <fstream>
#include<string>
#include<vector>
#include<sstream>
#include<iterator>
#include <cstdio>
#include <unistd.h>
#include <ios>

#include "data.h"
#include "functions.h"

using namespace std;



int main()
{

//----------------------------------------------------------------------------------------------------------------------------------
//                                          READ DATA FROM TEXT FILES AND STORE TO VECTORS
//----------------------------------------------------------------------------------------------------------------------------------
printf("Reading Data...\n");
    //  *************************************************
    //  READ KNOT VECTORS DATA
    //  Xi{} for Xi Direction
    //  Eta{} for Eta Direction
    //  Zeta{} for Zeta Direction
    //  *************************************************
    printf("\t Knot Vectors \n");
    vector<double> Xi,Eta,Zeta;
    Xi.clear();
    Eta.clear();
    Zeta.clear();
    text_file DATA;
    DATA.input_string_file_name();

        int size_KV_xi = DATA.determine_size(DATA.KV_xi_string);
        int size_KV_eta = DATA.determine_size(DATA.KV_eta_string);
        int size_KV_zeta = DATA.determine_size(DATA.KV_zeta_string);
    for(int i=0;i<size_KV_xi;i++){
        Xi.push_back(DATA.read_textfile_double(DATA.KV_xi_string,i));}
    for(int i=0;i<size_KV_eta;i++){
        Eta.push_back(DATA.read_textfile_double(DATA.KV_eta_string,i));}
    for(int i=0;i<size_KV_zeta;i++){
        Zeta.push_back(DATA.read_textfile_double(DATA.KV_zeta_string,i));}

    //  *************************************************
    //  PROPERTIES OF B-SPLINE BASIS FUNCTIONS
    //  *************************************************
        printf("\t Knot Vectors Properties \n");
        calc_bspline_prop(Xi,p,n);
        calc_bspline_prop(Eta,q,m);
        calc_bspline_prop(Zeta,r,l);

    //  *************************************************
    //  MODEL PROPERTIES
    //  *************************************************
        calc_nurbs_mod_prop(p,q,r,n,m,l,nel,nnp,nen,noctrlpt_patch,noxi,noeta,nozeta,nophyspts,Xi,Eta,Zeta);
        noctrlpt=calc_num_lines_textfile(DATA.ctrl_pt_coords_string);

    //  *************************************************
    //  READ CONTROL POINTS DATA
    //
    //  *************************************************
    printf("\t Control Points Data \n");

    ctrlpt_coordinates_t1.clear();
    ctrlpt_coordinates_t.clear();
    ctrlpt_coordinates_t0.clear();

        for(int j=0;j<noctrlpt;j++){
        for(int i=0;i<4;i++){
            //ctrlpt_coordinates_t1.push_back(DATA.read_textfile_double(DATA.ctrl_pt_coords_string,4*j+i));
            ctrlpt_coordinates_t1.push_back(DATA.read_textfile_double(DATA.ctrl_pt_coords_string,4*j+i));
            ctrlpt_coordinates_t.push_back(DATA.read_textfile_double(DATA.ctrl_pt_coords_string,4*j+i));
            ctrlpt_coordinates_t0.push_back(DATA.read_textfile_double(DATA.ctrl_pt_coords_string,4*j+i));
            }}

    ctrlpt_incidence.clear();
        for(int j=0;j<nel*nopatch;j++){
        for(int i=0;i<nen;i++){
            ctrlpt_incidence.push_back(DATA.read_textfile_int(DATA.ctrl_pt_inc_string,nen*j+i));
            }}

    bound_condition.clear();
        for(int j=0;j<noctrlpt;j++){
        for(int i=0;i<DOF;i++){
            bound_condition.push_back(DATA.read_textfile_int(DATA.bound_cond_string,DOF*j+i));
            }}

    //  *************************************************
    //  READ PATCH MEMBERS INCIDENCE
    //
    //  *************************************************
    printf("\t Patch Element Members Incidence \n");
    vector<int> patch_incidence;
    patch_incidence.clear();
        for(int i=0;i<nopatch;i++){
        for(int j=0;j<nel;j++){
            patch_incidence.push_back(DATA.read_textfile_int(DATA.patch_inc_string,nel*i+j));
        }}

    //  *************************************************
    //  READ NODAL COORDINATES
    //  *************************************************
    printf("\t Coordinates of Physical Points \n");
   // vector<double> phypt_coordinates;
    //phypt_coordinates.clear();
    //for(int j=0;j<nophyspts;j++){
      //  for(int i=0;i<DOF;i++){
        //    phypt_coordinates.push_back(DATA.read_textfile_double(DATA.phys_pt_coords_string,DOF*j+i));
          //  }}

    printf("Successfully Read Data...\n");

    //  *************************************************
    //  READ POINT LOAD VECTOR
    //
    //  *************************************************
    printf("\t Reading Point Load Vector \n");
    text_file LOAD;
    vector<double> P_loads;
    P_loads.clear();
    for(int i=0;i<DOF*noctrlpt;i++){
        P_loads.push_back(LOAD.read_textfile_double("Point Loads.txt",i));
        }

    //  *************************************************
    //  IDENTIFY UNRESTRAINED DOF.
    //  DETERMINE NUMBER OF RESTRAINED AND UNRESTRAINED DOF.
    //  *************************************************

        det_dof_type(DOF,noctrlpt,unrestrained_dof_data,no_unrestrained_dof,no_restrained_dof,bound_condition);
        //for(int i=0;i<unrestrained_dof_data.size();i++){cout<<unrestrained_dof_data[i]<<endl;}
        //cout<<no_unrestrained_dof<<endl;
        //for(int i=0;i<noctrlpt;i++){printf("%d %d %d \n",bound_condition[3*i],bound_condition[3*i+1],bound_condition[3*i+2]);}

//----------------------------------------------------------------------------------------------------------------------------------
//                                                  PRINT MODEL PROPERTIES
//----------------------------------------------------------------------------------------------------------------------------------

printf("NURBS Model Properties \n");
printf("\t Number of Patches: %d \n",nopatch);
printf("\t Total number of Control Points: %d \n",noctrlpt);
printf("\t Number of Control Points per Patch: %d \n",noctrlpt_patch);
printf("\t Number of Elements per Patch: %d \n",nel);
printf("\t Number of Local Basis Functions with support to an element: %d \n",nen);
printf("\t Order of NURBS basis functions in Xi,Eta and Zeta direction: %d,%d,%d \n",p,q,r);
printf("\t Number of Basis Functions in Xi, Eta and Zeta direction: %d,%d,%d \n",n,m,l);


//----------------------------------------------------------------------------------------------------------------------------------
//                                                  PARAMETERS AT TIME t=0
//----------------------------------------------------------------------------------------------------------------------------------

    //  *************************************************
    //  LOAD AND INITIALIZE GLOBAL PARAMETERS AT TIME, t=0
    //  *************************************************

        load_global_param(u_global_t,u_g_t1,Fext_global_t,Fext_g_t,K_global_t,DOF,noctrlpt,no_unrestrained_dof,u_global_t1,u_g_t);

        //load_postproc_param(strain_t,noctrlpt);

        load_D_e(D_e,bulk_mod,rank_4_unit,shear_mod,dev_tens);

        //for(int i=0;i<36;i++){cout<<D_e[i]<<endl;}

        step_load=0;

    //  *************************************************
    //  Clear txt files values
    //  *************************************************
        std::ofstream step_load_file_res1("Deformation Results.txt",std::ios::trunc);
        step_load_file_res1.close();
        std::ofstream step_load_file_res2("Stress Results.txt",std::ios::trunc);
        step_load_file_res2.close();

//----------------------------------------------------------------------------------------------------------------------------------
//                                                      LOOP THROUGH STEP LOAD
//----------------------------------------------------------------------------------------------------------------------------------
while(step_load<1)
{
    step_load++;
printf("--------------------------------------- STEP LOAD %d ---------------------------------------\n\n",step_load);

//----------------------------------------------------------------------------------------------------------------------------------
//                                                  GLOBAL STIFFNESS MATRIX AT TIME t
//----------------------------------------------------------------------------------------------------------------------------------

printf("Calculate Global Stiffness Matrix using converged parameters at time t...\n");
    //  *************************************************
    //  1. INITIALIZE ELEMENT AND KNOT POINT PARAMETERS
    //  *************************************************


    //  ********************************************************
    //  2. LOOPING THROUGH ELEMENTS, KNOT VECTORS AND GAUSS POINTS
    //  TO CALCULATE CONTRIBUTION OF AN ELEMENT TO THE GLOBAL
    //  STIFFNESS MATRIX
    //  ********************************************************
        fill(K_global_t.begin(),K_global_t.end(),0);

        //  ********************************************************
        //  2.1 LOOP THROUGH PATCHES - (START)
        //  ********************************************************
            for(int patch=0;patch<nopatch;patch++){
            //cout<<"PATCH:"<<patch+1<<endl;
            NURBS_coords.clear();
            build_NURBS_coords(Xi,Eta,Zeta,NURBS_coords);

            int element=0;
            initialize_knot_point_par(xi,eta,zeta,xi_index,eta_index,zeta_index);

        //  ********************************************************
        //  2.1 LOOP THROUGH ELEMENTS - (START)
        //  ********************************************************
            //printf("\t 1.1 Looping through elements \n");
            for(int el=0;el<nel;el++){

            element++;
            //printf("\t \t Element %d \n",element);

            int element_pointer=0;
            element_pointer=patch_incidence[nel*patch+el];
            //cout<<element_pointer<<endl;
            //element=patch_incidence[nel*patch+el];

                //  ********************************************************
                //  2.1.1 DETERMINE NURBS COORDINATES
                //  ********************************************************
                    //printf("\t\t 1.2 Looping through Knot Vectors \n");

                    xi_index=NURBS_coords[3*(element-1)];
                    eta_index=NURBS_coords[3*(element-1)+1];
                    zeta_index=NURBS_coords[3*(element-1)+2];
                    //cout<<xi_index<<" "<<eta_index<<" "<<zeta_index<<endl;
        //  ********************************************************
        //  2.2 LOOP THROUGH GAUSS POINTS - (START)
        //  ********************************************************
            //printf("\t\t 1.2 Looping through Gauss Points \n");
            int gauss_pt_counter=0;
            for(int gp_z=0;gp_z<NGP;gp_z++){
            for(int gp_y=0;gp_y<NGP;gp_y++){
            for(int gp_x=0;gp_x<NGP;gp_x++){
            gauss_pt_counter++;

                //  ********************************************************
                //  2.2.1 DETERMINE PARAMETRIC COORDINATE
                //  ********************************************************
                    xi=calc_parametric_coord(Xi,xi_index,gauss_pts[gp_x]);
                    eta=calc_parametric_coord(Eta,eta_index,gauss_pts[gp_y]);
                    zeta=calc_parametric_coord(Zeta,zeta_index,gauss_pts[gp_z]);
                    //printf("%f %f %f \n",xi,eta,zeta);
        //  ********************************************************
        //  2.3 CALCULATE B-SPLINE BASIS FUNCTIONS AND THEIR DERIVATIVES
        //  ********************************************************
            form_B_splines(N,N_der,xi_index,p,xi,Xi,n);
            form_B_splines(M,M_der,eta_index,q,eta,Eta,m);
            form_B_splines(L,L_der,zeta_index,r,zeta,Zeta,l);

            //cout<<N[0]<<" "<<N[1]<<" "<<N[2]<<endl;
            //cout<<M[0]<<" "<<M[1]<<" "<<M[2]<<endl;
            //cout<<L[0]<<" "<<L[1]<<" "<<L[2]<<endl;
        //  ********************************************************
        //  2.4 DETERMINE NURBS BASIS FUNCTIONS AND THEIR DERIVATIVES
        //      w.r.t PHYSICAL COORDINATES
        //  ********************************************************
            Shape_func_1(ctrlpt_incidence,ctrlpt_coordinates_t,R,dR_dxi,N,M,L,N_der,M_der,L_der,p,q,r,element_pointer,nen,DOF,dR_dx,dx_dxi);
            //cout<<N_der[0]<<" "<<N_der[1]<<" "<<N_der[2]<<endl;
            //for(int i=0;i<nen;i++){printf("%f %f %f \n",dR_dx[3*i],dR_dx[3*i+1],dR_dx[3*i+2]);}

        //  ********************************************************
        //  2.5 CALCULATE DETERMINANT OF THE JACOBIAN MATRIX
        //      (Mapping from parametric to physical coordinates)
        //  ********************************************************
            J=Shape_function_2(DOF,xi_index,eta_index,zeta_index,Xi,Eta,Zeta,dx_dxi);
            //cout<<J<<endl;

        //  ********************************************************
        //  2.6 CALCULATE CONTRIBUTION TO GLOBAL STIFFNESS MATRIX
        //  ********************************************************
            //printf("\t\t Calculating Contribution of Gauss Point %d to Global Stiffness Matrix\n",gauss_pt_counter);
            K_local_mod.clear();
                for(int i=0;i<DOF*DOF;i++){
                    K_local_mod.push_back(0);
                    }

                int ctrl_pt_ii, ctrl_pt_jj;
                    for(int ii=0;ii<nen;ii++){
                        ctrl_pt_ii = ctrlpt_incidence[nen*(element_pointer-1)+ii]-1;
                    for(int jj=0;jj<nen;jj++){
                        ctrl_pt_jj=ctrlpt_incidence[nen*(element_pointer-1)+jj]-1;

                        build_k_local(bulk_mod,shear_mod,dR_dx,J*gauss_pts_wts[gp_x]*gauss_pts_wts[gp_y]*gauss_pts_wts[gp_z],
                                        DOF,K_local_mod,ii,jj);

                        for(int aa=0;aa<DOF;aa++){
                        for(int bb=0;bb<DOF;bb++){
                            K_global_t[DOF*noctrlpt*(DOF*ctrl_pt_ii+aa)+DOF*ctrl_pt_jj+bb]+=K_local_mod[DOF*aa+bb];
                    }}}}

            //cout<<K_global_t[0]<<" "<<K_global_t[1]<<" "<<K_global_t[2]<<endl;

        //  ********************************************************
        //  2.7 LOOP THROUGH GAUSS POINTS - (END)
        //  ********************************************************
            }}}

        //  ********************************************************
        //  2.8 LOOP THROUGH ELEMENTS - (END)
        //  ********************************************************
            }
            //}}}

        //  ********************************************************
        //  2.8 LOOP THROUGH PATCHES - (END)
        //  ********************************************************
            }

 printf("Finished Calculating Global Stiffness Matrix for time t...\n");

//for(int i=0;i<noctrlpt*DOF;i++)
//{for(int j=0;j<noctrlpt*DOF;j++)
//{
  //cout<<K_global_t[noctrlpt*DOF*i+j]<<" ";
//}
//cout<<endl;
//}
            std::ofstream initial_kglobal("Global Stiffness Matrix.txt",std::ios::out);
            //std::ofstream outf13("Fext_g.txt",std::ios::out);
            //outfl2<<"Reduced Matrix"<<"\n";
            for (int i=0; i<noctrlpt*3; i++)
            {
                for (int j=0; j<noctrlpt*3; j++)
                {
                initial_kglobal<<left<<K_global_t[noctrlpt*3*i+j]<<" ";
                }

                initial_kglobal<<"\n";

            }
            initial_kglobal.close();


//cout<<"sad:"<<K_global_t[1218815]<<endl;


//----------------------------------------------------------------------------------------------------------------------------------
//                                          DETERMINE EXTERNAL FORCE VECTOR AT TIME t+1
//----------------------------------------------------------------------------------------------------------------------------------
printf("Calculating Initial Incremental Displacements ( delta_u_global_t1 )...\n");

    printf("\t Forming Initial Reduced Global Stiffness Matrix and Reduced Global Force Vector ( K_gg_t , delta_F_g_t1 ) \n");

    //  *************************************************
    //  1. DETERMINE REDUCED GLOBAL STIFFNESS MATRIX AND
    //     REDUCED GLOBAL LOAD VECTOR
    //  *************************************************

            extr_reduced_vector(K_gg_t,delta_F_g_t1,P_loads,bound_condition,K_global_t,DOF,noctrlpt,no_unrestrained_dof);
        printf("\t Calculating Initial Incremental Displacements at unrestrained control points, ( delta_u_g_t1 ), using Conjugate Gradient Method  \n");

        //cout<<no_unrestrained_dof<<endl;
        //for(int i=0;i<no_unrestrained_dof;i++){
            //for(int j=0;j<no_unrestrained_dof;j++){
                //cout<<K_gg_t[no_unrestrained_dof*i+j]<<" ";
            //}
            //cout<<endl;
        //}

    std::ofstream outfl2("Reduced Stiffness Matrix.txt",std::ios::out);
    std::ofstream outf13("Fext_g.txt",std::ios::out);
    //outfl2<<"Reduced Matrix"<<"\n";
    for (int i=0; i<no_unrestrained_dof; i++)
    {
        for (int j=0; j<no_unrestrained_dof; j++)
        {
        outfl2<<left<<K_gg_t[no_unrestrained_dof*i+j]<<" ";
        }
        outf13<<left<<delta_F_g_t1[i];
        outfl2<<"\n";
        outf13<<"\n";
    }
    outfl2.close();
    outf13.close();



    //  *************************************************
    //  2. UPDATE REDUCED GLOBAL FORCE VECTOR AT TIME t+1
    //  *************************************************
        for(int i = 0;i<no_unrestrained_dof;i++){
        Fext_g_t[i]+=delta_F_g_t1[i];}

        //for(int i = 0;i<no_unrestrained_dof;i++){
        //cout<<Fext_g_t[i]<<endl;}

    printf("\t Finished calculating Initial Incremental Displacements \n");
//----------------------------------------------------------------------------------------------------------------------------------
//                                            CALCULATE INCREMENTAL DISPLACEMENTS
//----------------------------------------------------------------------------------------------------------------------------------
printf("Iterate initial incremental displacements ,( delta_u_g_t1 ), using Newton-Raphson method... \n");
    //  *************************************************
    //  1. DETERMINE INITIAL PARAMETERS FOR NEWTON-RAPHSON METHOD
    //
    //  *************************************************
        printf("\t Setting up initial parameters for Newton-Raphson Method \n");
        //  ********************************************************
        //  1.1 DETERMINE INITIAL RESIDUAL
        //  ********************************************************
            calc_residual(K_gg_t,r_initial,Fext_g_t,u_g_t1,norm_initial,no_unrestrained_dof);

        //  ********************************************************
        //  1.2 SET CURRENT RESIDUAL EQUAL TO INCREMENTAL LOAD
        //  ********************************************************
            set_curr_residual(r_current,no_unrestrained_dof,delta_F_g_t1,norm_current);

        //  ********************************************************
        //  1.3 SET INITIAL ERROR TO 1, ITERATION TO 0 AND MAXIMUM ITERATION
        //  ********************************************************
            iterate=0;
            TOL=1;
            alpha=1;

            norm_previous=0;
            norm_current=norm_initial;
            //cout<<"norm_initial:"<<norm_initial<<endl;
            //cout<<"norm_current:"<<norm_current<<endl;

    //  *************************************************
    //  2. PERFORM NEWTON RAPHSON.
    //  *************************************************

    printf("\t Performing Newton-Raphson Method \n");
        //  ********************************************************
        //  2.1 SET CONDITIONS FOR CONVERGENCE.
        //  ********************************************************
            while(TOL>pow(10,-8)&&iterate<100){
            //cout<<"Iteration : "<<iterate<<endl;

        //  ********************************************************
        //  2.2 DETERMINE REDUCED GLOBAL DISPLACEMENT VECTOR.
        //  ********************************************************

            //  ********************************************************
            //  2.2.1 PERFORM PCGM TO DETERMINE DISPLACEMENTS CORRESPONDING
            //        TO CURRENT RESIDUAL.
            //  ********************************************************
                delta_u_g_t1.clear();
                delta_u_g_t1.resize(no_unrestrained_dof);
                vector<double> initial_value;
                initial_value.resize(no_unrestrained_dof,0);
                //cout<<no_unrestrained_dof<<endl;
                //cout<<"sad"<<K_gg_t[1149184]<<endl;
                    perform_PCGM(K_gg_t,no_unrestrained_dof,r_current,pow(10,-8),150,delta_u_g_t1,initial_value);
                    //perform_gauss_jordan(K_gg_t,no_unrestrained_dof,r_current,delta_u_g_t1);

            //  ********************************************************initial_value
            //  2.2.2 UPDATE REDUCED GLOBAL DISPLACEMENT VECTOR BY ADDING
            //        DELTA DISPLACEMENTS CORRESPONDING TO CURRENT RESIDUAL.
            //  ********************************************************
            norm_previous=norm_current;
            do
                {
                //cout<<"\t LINE SEARCH"<<endl;
                //cout<<"\t alpha :"<<alpha<<endl;

                    for(int i=0;i<no_unrestrained_dof;i++){
                        u_g_t1[i]=u_g_t[i]+alpha*delta_u_g_t1[i];
                        //cout<<u_g_t1[i]<<endl;
                        //u_global_t[DOF*(unrestrained_dof_data[2*i]-1)+unrestrained_dof_data[2*i+1]]=u_g_t1[i];
                        }
                        //cout<<u_g_t1[0]<<" "<<u_g_t1[1]<<" "<<u_g_t1[0]<<endl;

            //  ********************************************************
            //  2.3 UPDATE CONTROL POINTS.
            //  ********************************************************

                update_ctrl_pt(ctrlpt_coordinates_t1,delta_u_g_t1,no_unrestrained_dof,unrestrained_dof_data,alpha,ctrlpt_coordinates_t);
                //for(int i=0;i<4*noctrlpt;i++){cout<<ctrlpt_coordinates_t1[i]<<endl;}
                //cout<<"ctrl pt"<<endl;
                //cout<<ctrlpt_coordinates_t1[48]<<" "<<ctrlpt_coordinates_t1[49]<<" "<<ctrlpt_coordinates_t1[50]<<endl;


            //  ********************************************************
            //  2.4 DETERMINE GLOBAL STIFFNESS MATRIX FOR TIME t+1 USING
            //      UPDATED VALUES OF CONTROL POINTS.
            //  ********************************************************

                //  ********************************************************
                //  2.4.1 LOOPING THROUGH ELEMENTS, KNOT VECTORS AND GAUSS POINTS
                //        TO CALCULATE CONTRIBUTION OF AN ELEMENT TO THE GLOBAL
                //        STIFFNESS MATRIX
                //  ********************************************************
                    fill(K_global_t.begin(),K_global_t.end(),0);

                    //  ********************************************************
                    //  2.1 LOOP THROUGH PATCHES - (START)
                    //  ********************************************************
                    for(int patch=0;patch<nopatch;patch++){
                    NURBS_coords.clear();
                    build_NURBS_coords(Xi,Eta,Zeta,NURBS_coords);
                    xi_index=0,eta_index=0,zeta_index=0;
                    element=0;

                    //  ********************************************************
                    //  2.4.1.1 LOOP THROUGH ELEMENTS - (START)
                    //  ********************************************************
                        //for(int zeta_i=0;zeta_i<l-r;zeta_i++){
                        //for(int eta_i=0;eta_i<m-q;eta_i++){
                        //for(int xi_i=0;xi_i<n-p;xi_i++){
                        for(int el=0;el<nel;el++){
                        element++;
                        int element_pointer=0;
                        element_pointer=patch_incidence[nel*patch+el];

                        //  ********************************************************
                        //  2.4.1.1.1 DETERMINE NURBS COORDINATES
                        //  ********************************************************
                            //xi_index=xi_i+p;
                            //eta_index=eta_i+q;
                            //zeta_index=zeta_i+r;initial_values
                            xi_index=NURBS_coords[3*(element-1)];
                            eta_index=NURBS_coords[3*(element-1)+1];
                            zeta_index=NURBS_coords[3*(element-1)+2];

                    //  ********************************************************
                    //  2.4.2 LOOP THROUGH GAUSS POINTS - (START)
                    //  ********************************************************
                        int gauss_pt_counter=0;
                        for(int gp_z=0;gp_z<NGP;gp_z++){
                        for(int gp_y=0;gp_y<NGP;gp_y++){
                        for(int gp_x=0;gp_x<NGP;gp_x++){
                        gauss_pt_counter++;

                        //  ********************************************************
                        //  2.4.2.1 DETERMINE PARAMETRIC COORDINATE
                        //  ********************************************************
                            xi=calc_parametric_coord(Xi,xi_index,gauss_pts[gp_x]);
                            eta=calc_parametric_coord(Eta,eta_index,gauss_pts[gp_y]);
                            zeta=calc_parametric_coord(Zeta,zeta_index,gauss_pts[gp_z]);

                    //  ********************************************************
                    //  2.4.3 CALCULATE B-SPLINE BASIS FUNCTIONS AND THEIR DERIVATIVES
                    //  ********************************************************
                        form_B_splines(N,N_der,xi_index,p,xi,Xi,n);
                        form_B_splines(M,M_der,eta_index,q,eta,Eta,m);
                        form_B_splines(L,L_der,zeta_index,r,zeta,Zeta,l);

                    //  ********************************************************
                    //  2.4.4 DETERMINE NURBS BASIS FUNCTIONS AND THEIR DERIVATIVES
                    //      w.r.t PHYSICAL COORDINATES
                    //  ********************************************************
                        Shape_func_1(ctrlpt_incidence,ctrlpt_coordinates_t1,R,dR_dxi,N,M,L,N_der,M_der,L_der,p,q,r,element_pointer,nen,DOF,dR_dx,dx_dxi);

                    //  ********************************************************
                    //  2.4.5 EVALUATE IF GAUSS POINT IS ELASTIC OR PLASTIC
                    //  ********************************************************

                        //calc_strain(du_dx_point,dR_dx,u_global_t,DOF,nen,ctrlpt_incidence,element,strain_point);
                        //calc_stress(strain_point,stress_point,D_e);
                        //calc_dev_stress(stress_point,dev_stress_point);
                        //evaluate_yield_func(dev_stress_point,fy);

                    //  ********************************************************
                    //  2.4.6 CALCULATE DETERMINANT OF THE JACOBIAN MATRIX
                    //        (Mapping from parametric to physical coordinates)
                    //  ********************************************************
                        J=Shape_function_2(DOF,xi_index,eta_index,zeta_index,Xi,Eta,Zeta,dx_dxi);

                    //  ********************************************************
                    //  2.4.7 CALCULATE CONTRIBUTION TO GLOBAL STIFFNESS MATRIX
                    //  ********************************************************
                        K_local_mod.clear();
                            for(int i=0;i<DOF*DOF;i++){
                                K_local_mod.push_back(0);
                                }

                            int ctrl_pt_ii, ctrl_pt_jj;
                                for(int ii=0;ii<nen;ii++){
                                    ctrl_pt_ii = ctrlpt_incidence[nen*(element_pointer-1)+ii]-1;
                                for(int jj=0;jj<nen;jj++){
                                    ctrl_pt_jj=ctrlpt_incidence[nen*(element_pointer-1)+jj]-1;

                                    build_k_local(bulk_mod,shear_mod,dR_dx,J*gauss_pts_wts[gp_x]*gauss_pts_wts[gp_y]*gauss_pts_wts[gp_z],
                                                          DOF,K_local_mod,ii,jj);

                                    for(int aa=0;aa<DOF;aa++){
                                    for(int bb=0;bb<DOF;bb++){
                                        K_global_t[DOF*noctrlpt*(DOF*ctrl_pt_ii+aa)+DOF*ctrl_pt_jj+bb]+=K_local_mod[DOF*aa+bb];
                                }}}}

                    //  ********************************************************
                    //  2.4.8 LOOP THROUGH GAUSS POINTS - (END)
                    //  ********************************************************
                        }}}

                    //  ********************************************************
                    //  2.4.8 LOOP THROUGH ELEMENTS - (END)
                    //  ********************************************************
                        //}}}
                        }
                    //  ********************************************************
                    //  2.8 LOOP THROUGH PATCHES - (END)
                    //  ********************************************************
                        }

            //  ********************************************************
            //  2.5 DETERMINE REDUCED GLOBAL STIFFNESS MATRIX
            //  ********************************************************
                K_gg_t.clear();
                delta_F_g_t1.clear();
                    extr_reduced_vector_Kgg(K_gg_t,bound_condition,K_global_t,DOF,noctrlpt);

            //  ********************************************************
            //  2.6 COMPUTE CURRENT RESIDUAL
            //  ********************************************************

                //cout<<"\t f1: "<<norm_previous<<endl;
                calc_residual(K_gg_t,r_current,Fext_g_t,u_g_t1,norm_current,no_unrestrained_dof);
                //cout<<"\t f2: "<<norm_current<<endl;

            //  ********************************************************
            //  2.6 COMPUTE ALPHA
            //  ********************************************************

                if(norm_current>norm_previous)
                {
                    alpha=pow(alpha,2)*norm_previous/(2*(norm_current+alpha*norm_previous-norm_previous));
                }
                else{alpha=1;}

                }
            while(norm_current>norm_previous);

        //  ********************************************************
        //  2.6 UPDATE CONVERGENCE PARAMETERS
        //  ********************************************************
            iterate++;

            for(int i=0;i<no_unrestrained_dof;i++){u_g_t[i]=u_g_t1[i];}
            for(int i=0;i<4*noctrlpt;i++){ctrlpt_coordinates_t[i]=ctrlpt_coordinates_t1[i];}
            //for(int i=0;i<no_unrestrained_dof;i++){initial_value[i]=alpha*delta_u_g_t1[i];}

            TOL=norm_current/norm_initial;
            //cout<<"\t N-R TOL: "<<fixed<<setprecision(18)<<TOL<<endl;

        //  ********************************************************
        //  2.7 EXIT IF CONVERGENCE IS MET
        //  ********************************************************

            }
            //cout<<"TOL: "<<fixed<<setprecision(18)<<TOL<<endl;

            for(int i=0;i<no_unrestrained_dof;i++){
                    u_global_t1[DOF*(unrestrained_dof_data[2*i]-1)+unrestrained_dof_data[2*i+1]]=u_g_t1[i];
                    }
        printf("\t Convergence criterion met with %d iterations , || r_current || / || r_initial || = %.16f \n",iterate,TOL);
        printf("Finished calculating incremental displacements \n");

//----------------------------------------------------------------------------------------------------------------------------------
//                             COMPUTE PARAMETERS FOR TIME t+1 USING CONVERGED VALUES OF INCREMENTAL DISPLACEMENTS
//----------------------------------------------------------------------------------------------------------------------------------

printf("Calculate Parameters for time t+1 using converged values of incremental displacements... \n");
    //  *************************************************
    //  1. UPDATE GLOBAL DISPLACEMENTS FOR TIME t+1
    //  *************************************************
    printf("\t Update Global displacements for time t+1\n");
    printf("\t Update Control Points coordinates for time t+1\n");
    //for(int j=0;j<no_unrestrained_dof;j++){
        //u_global_t[DOF*(unrestrained_dof_data[2*j]-1)+unrestrained_dof_data[2*j+1]]=u_g_t1[j];}

            //for(int i=0;i<noctrlpt;i++){
            //cout<<i+1<<" "<<fixed<<setprecision(8)<<u_global_t1[DOF*i]<<" "<<u_global_t1[DOF*i+1]<<" "<<u_global_t1[DOF*i+2]<<endl;}
            //double sum1=0;
            //for(int i=0;i<noctrlpt*DOF;i++){
                //sum1=0;
                //for(int j=0;j<noctrlpt*DOF;j++){
                  //  sum1+=K_global_t[noctrlpt*DOF*i+j]*u_global_t[j];
                //}
                //cout<<sum1<<endl;
            //}
//----------------------------------------------------------------------------------------------------------------------------------
//                                    COMPUTE FOR STRESSES AT TIME t+1 USING CONVERGED VALUES OF u_t+1
//----------------------------------------------------------------------------------------------------------------------------------
int counter =0;
int phys_node_num=0;
int start_num=0;
int start_elem_eta=0;
vector<int> stress_strain_divisor;
stress_strain_divisor.resize(108,0);

vector<double> stress_t1;
stress_t1.resize(6*108,0);
vector<double> strain_t1;
strain_t1.resize(6*108,0);

vector<double> disp_t1;
disp_t1.resize(3*108,0);

for(int patch=0;patch<nopatch;patch++){

            for(int k=0;k<nozeta;k++){
                counter=0;
                element=0;
                start_num=patch*(noxi*noeta-noxi);

            for(int j=0;j<noeta;j++){
            for(int i=0;i<noxi;i++){

                counter++;
                phys_node_num=noxi*(nopatch*noeta-1)*k+counter+start_num;
                //cout<<"phy:"<<phys_node_num<<endl;
                //cout<<counter<<endl;

                xi=Xi[i+p];
                eta=Eta[j+q];
                zeta=Zeta[k+r];

                if(zeta==1){
                    element=(noxi-1)*(noeta-1)*(nozeta-1)-(noxi-1)*(noeta-1)+noeta*j+i;}

                if(eta==1){
                    element=(noxi-1)*(noeta-1)-(noxi-1)+i;}

                if((xi==1)){
                    element=element;}

                else{element++;}

                int element_pointer=0;
                element_pointer=patch_incidence[nel*patch+element-1];


                //printf("xi,eta,zeta,element: %f %f %f %d\n",xi,eta,zeta,element);
                //printf("node number, element, element pointer: %d %d %d \n",phys_node_num,element,element_pointer);


                stress_strain_divisor[phys_node_num-1]+=1;


                xi_index=findspan(noxi,p,xi,Xi);
                eta_index=findspan(noeta,q,eta,Eta);
                zeta_index=findspan(nozeta,r,zeta,Zeta);

                form_B_splines(N,N_der,xi_index,p,xi,Xi,n);
                form_B_splines(M,M_der,eta_index,q,eta,Eta,m);
                form_B_splines(L,L_der,zeta_index,r,zeta,Zeta,l);

                Shape_func_1(ctrlpt_incidence,ctrlpt_coordinates_t1,R,dR_dxi,N,M,L,N_der,M_der,L_der,p,q,r,element_pointer,nen,DOF,dR_dx,dx_dxi);

                calc_displacement(ctrlpt_incidence,u_global_t1,R,element_pointer,nen,disp_point);
                calc_strain(du_dx_point,dR_dx,u_global_t1,DOF,nen,ctrlpt_incidence,element_pointer,strain_point);
                calc_stress(strain_point,stress_point,D_e);

                for(int aa=0;aa<6;aa++){
                    stress_t1[6*(phys_node_num-1)+aa]+=stress_point[aa];
                    strain_t1[6*(phys_node_num-1)+aa]+=strain_point[aa];
                }

                for(int aa=0;aa<3;aa++){
                    disp_t1[3*(phys_node_num-1)+aa]=disp_point[aa];
                }

                double stress_equiv=0;
                double strain_equiv=0;

                if(phys_node_num==108)
                {
                    stress_equiv=calc_equiv(stress_point);
                    strain_equiv=calc_equiv(strain_point);
                    create_report_stress_strain_pt(step_load,1,stress_equiv,strain_equiv);
                }
                //stress_equiv=0;
                //strain_equiv=0;

                //if(phys_node_num==30)
                //{
                  //  stress_equiv=calc_equiv(stress_point);
                   // strain_equiv=calc_equiv(strain_point);
                    //create_report_stress_strain_pt(step_load,30,stress_equiv,strain_equiv);
                //}

            }}}

    }

 create_report_stress(step_load,108,stress_t1,stress_strain_divisor);
 create_report_deform(step_load,108,disp_t1);



//----------------------------------------------------------------------------------------------------------------------------------
//                                                          END STEP LOAD
//----------------------------------------------------------------------------------------------------------------------------------




}
printf("--------------------------------------- END OF STEP LOAD %d ---------------------------------------\n\n",step_load);


//----------------------------------------------------------------------------------------------------------------------------------
//                                                      POST PROCESSING
//----------------------------------------------------------------------------------------------------------------------------------





    return 0;
}
