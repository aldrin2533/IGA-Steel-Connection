#pragma once
#include <vector>
#include<string>


using namespace std;

//----------------------------------------------------------------------------------------------------------------------------------
//                                              CONSTANT PARAMETERS FOR ANALYSIS
//----------------------------------------------------------------------------------------------------------------------------------
    //  *************************************************
    //  GAUSSIAN POINTS PROPERTIES
    //  DOF                 Degrees of freedom
    //  NGP                 Number of Gauss Points
    //  gauss_pts[]         Gauss Points.
    //  gauss_pts_wts[]     Corresponding weights of gauss
    //                      points.
    //  *************************************************
    int const DOF = 3;
    int const NGP = 3;
    double gauss_pts[NGP]={-0.7745966692,0,0.7745966692};
    double gauss_pts_wts[NGP]={0.555555556,0.8888888889,0.555555556};

//----------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------
//                                                   NURBS MODEL PROPERTIES
//----------------------------------------------------------------------------------------------------------------------------------
    //  *************************************************
    //  MODEL PROPERTIES
    //  p,q,r                   Polynomial orders of Basis functions for xi,eta and zeta direction.
    //  n,m,l                   Number of Basis functions for xi,eta and zeta direction.
    //  nel                     Total number of elements.
    //  nnp                     Total number global basis functions.
    //  nen                     Number of local basis functions per element.
    //  noctrlpt                Total number of control points.
    //  noctrlpt_patch          Number of Control Points per Patch
    //  noxi,noeta,nozeta       Number of Knot points excluding repeated values.
    //  nophyspts               Total number of physical points.
    //  nopatch                 Total number of patches
    //  ctrlpt_coordinates_t    Control Point Coordinates at time t+1
    //  ctrlpt_coordinates      Control Point Coordinates at time t
    //  ctrlpt_incidence        Control point to element connectivity array
    //  bound_condition         Boundary Condition  array (1 restrained 0 unrestrained)
    //  *************************************************
    int p,q,r;
    int n,m,l;
    int nel;
    int nnp;
    int nen;
    int noctrlpt;
    int noctrlpt_patch;
    int noxi;
    int noeta;
    int nozeta;
    int nophyspts;
    int nopatch=2;

    vector<double> ctrlpt_coordinates_t1;
    vector<double> ctrlpt_coordinates_t;
    vector<double> ctrlpt_coordinates_t0;
    vector<int> ctrlpt_incidence;
    vector<int> bound_condition;

    //  *************************************************
    //  DEGREES OF FREEDOM, SUPPORTS
    //  no_dof_release              Total number of unrestrained DOF.
    //  unrestrained_dof_data       Vector to store reference control point and unrestrained DOF
    //  *************************************************
    int no_unrestrained_dof;
    int no_restrained_dof;
    vector<int> unrestrained_dof_data;
    vector<int> IEN;
    vector<int> NURBS_coords;

//----------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------
//                                                      MATERIAL PROPERTIES
//----------------------------------------------------------------------------------------------------------------------------------
    //  *************************************************
    //  MATERIAL PROPERTIES
    //  E           Modulus of Elasticity.
    //  v           Poisson's Ratio.
    //  bulk_mod    Bulk modulus
    //  shear_mod   Shear Modulus
    //  fy          Yield Strength.
    //  *************************************************
    double E = 200000;
    double v = 0.3;
    double bulk_mod = E/(3*(1-2*v));
    double shear_mod=E/(2*(1+v));
    vector<double> D_e(36);
    double fy=275;

//----------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------
//                                                       OPERATOR TENSORS
//----------------------------------------------------------------------------------------------------------------------------------

    vector<int> rank_2_unit{1,1,1,0,0,0};

    vector<int> rank_4_unit{1,1,1,0,0,0,
                            1,1,1,0,0,0,
                            1,1,1,0,0,0,
                            0,0,0,0,0,0,
                            0,0,0,0,0,0,
                            0,0,0,0,0,0,};

    vector<double> dev_tens{0.66666667,-0.33333333,-0.33333333,0,0,0,
                        -0.33333333,0.66666667,-0.33333333,0,0,0,
                        -0.33333333,-0.33333333,0.66666667,0,0,0,
                        0,0,0,0.5,0,0,
                        0,0,0,0,0.5,0,
                        0,0,0,0,0,0.5};

//----------------------------------------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------
//                                                  PARAMETERS FOR ANALYSIS
//----------------------------------------------------------------------------------------------------------------------------------

    //  *************************************************
    //  GLOBAL PARAMETERS
    //  u_global_t              Vector to store all displacements including restrained DOF at time t.
    //  u_g_t                   Vector to store displacements with unrestrained DOF.
    //  Fext_global_t           Vector to store all external load.
    //  Fext_g_t                Reduced vector to store external loads corresponding to u_g_t.
    //  K_global_t              Vector to store global stiffness matrix.
    //  K_gg_t                  Reduced global stiffness matrix corresponding to u_g_t.
    //  delta_F_g_t1            Vector containing incremental external load.
    //  u_global_t1             Vector to store all displacements including restrained DOF at time t+1
    //  *************************************************
    vector<double> u_global_t;
    vector<double> u_g_t1;
    vector<double> u_g_t;
    vector<double> Fext_global_t;
    vector<double> Fext_g_t;
    vector<double> K_global_t;
    vector<double> K_gg_t;
    vector<double> delta_F_g_t1;
    vector<double> u_global_t1;

    //  *************************************************
    //  LOCAL PARAMETERS
    //  element                             Element pointer.
    //  xi,eta,zeta                         Knot points.
    //  xi_index,eta_index,zeta_index       Knot Index pointers.
    //  N[],M[],L[]                         Univariate B-spline Basis functions.
    //  N_der[],M_der[],L_der[]             Derivative of B-spline Basis functions.
    //  R[]                                 NURBS basis functions.
    //  dR_dxi[][]                          Derivative of NURBS basis functions w.r.t. parametric coordinates.
    //  dR_dx[][]                           Derivative of NURBS basis functions w.r.t. physical coordinates.
    //  *************************************************
    int element;
    int xi_index,eta_index,zeta_index;
    double xi,eta,zeta;
    vector<double> N,M,L;
    vector<double> N_der,M_der,L_der;
    vector<double> R;
    vector<double> dR_dxi;
    vector<double> dR_dx;
    vector<double> dx_dxi;
    vector<double> K_local_mod;
    double J;

    //  *************************************************
    //  NEWTON RAPHSON PARAMETERS
    //  r_initial                           Initial residual.
    //  norm_initial                        Norm of initial residual.
    //  *************************************************
    vector<double> r_initial;
    vector<double> r_current;
    double norm_initial;
    double norm_current;
    double norm_previous;

    int iterate;
    double TOL;
    double alpha;

    vector<double> delta_u_g_t1;

    //  *************************************************
    //  LOADING PARAMETERS
    //  *************************************************

    int step_load;

//----------------------------------------------------------------------------------------------------------------------------------
//                                                  PARAMETERS FOR POST - ANALYSIS
//----------------------------------------------------------------------------------------------------------------------------------

    vector<double> u_point{0,0,0};
    vector<double> du_dx_point{0,0,0,0,0,0,0,0,0};
    vector<double> du_dx{0,0,0,0,0,0,0,0,0};
    vector<double> du_dxi{0,0,0,0,0,0,0,0,0};
    vector<double> strain_point{0,0,0,0,0,0};
    vector<double> stress_point{0,0,0,0,0,0};
    vector<double> strain_t;
    vector<double> dev_stress_point{0,0,0,0,0,0};
    vector<double> disp_point{0,0,0};
    vector<double> physical_point{0,0,0};




