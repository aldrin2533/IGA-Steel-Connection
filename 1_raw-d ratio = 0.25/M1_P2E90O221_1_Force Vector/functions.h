#pragma once
#include <vector>
#include<string>

using namespace std;

//----------------------------------------------------------------------------------------------------------------------------------
//                                              CLASS TEMPLATE FOR READING TEXT FILES
//----------------------------------------------------------------------------------------------------------------------------------

class text_file{

public:

    void input_string_file_name();                                // function to write filenames of below strings
        string KV_xi_string,
               KV_eta_string,
               KV_zeta_string,
               ctrl_pt_coords_string,
               ctrl_pt_inc_string,
               bound_cond_string,
               phys_pt_coords_string;

    int determine_size(string);                                   // function to determine number of elements of an array from a text file


    double read_textfile_double(string,int);             // function to read text file with values of type double
    int read_textfile_int(string,int);                    // function to read text file with values of type double

};
//----------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------
//                                        MODEL'S B-SPLINE BASIS FUNCTIONS AND PROPERTIES
//
//  calc_bspline_prop                -function to determine order of the Knot vector.
//  calc_nurbs_mod_prop
//  det_dof_type
//  comp_parametric_coord           -function to determine the equivalent parametric coordinate given the gaussian point.
//  calc_basis_func                 -function to calculate the value of a basis function at a given knot point.
//  calc_Bspline_deriv              -function to calculate the value of the derivative of the basis function at a given knot point.
//  Shape_func_1                    -function to calculate the NURBS basis functions and their derivatives with respect to physical
//                                   coordinates.
//  build_k_local                   -
//
//----------------------------------------------------------------------------------------------------------------------------------

    int findspan(int,int,double,vector<double>);

    void build_IEN(int ,int ,int ,int ,int ,int,int,int,vector<int>& );

    void build_NURBS_coords(vector<double>,vector<double>,vector<double>,vector<int>&);

    void build_NURBS_coords(vector<double>,vector<double>,vector<int>&);

    void calc_bspline_prop(vector<double>,int&,int&);

    void calc_nurbs_mod_prop(int,int,int,int,int,int,int&,int&,int&,int&, int&,int&,int&,int&,vector<double>,vector<double>,vector<double>);

    void calc_nurbs_mod_prop(int,int,int,int,int&,int&,int&,int&, int&,int&,
                         int&,vector<double> ,vector<double>);

    void det_dof_type(int,int, vector<int>&,int&, int&,vector<int>);

    double calc_parametric_coord(vector<double>,int,double);

    double calc_basis_func(int,double,vector<double>,int,int);

    double calc_Bspline_deriv(int,double,std::vector<double>,int,int,int);

    void form_B_splines(vector<double>&,vector<double>&,int,int,double,vector<double>,int);

    void Shape_func_1(vector<int> ,vector<double> , vector<double>& ,vector<double>& ,vector<double> ,vector<double> ,vector<double> ,
                      vector<double> ,int ,int ,int ,int ,vector<double>& ,vector<double>& );

    double Shape_function_2(int,int,vector<double>,
                      vector<double>,vector<double>);

    void build_k_local(double,double,vector<double>,double,int,vector<double>&,int,int);

//----------------------------------------------------------------------------------------------------------------------------------
//                                          PRECONDITIONED CONJUGATE GRADIENT METHOD
//----------------------------------------------------------------------------------------------------------------------------------

    void perform_PCGM(vector<double> , int ,vector<double> ,double ,int ,vector<double>& );

//----------------------------------------------------------------------------------------------------------------------------------
//                                                  GLOBAL OPERATORS
//----------------------------------------------------------------------------------------------------------------------------------

    void load_global_param(vector<double>&,vector<double>&,vector<double>&,vector<double>&,vector<double>&,int,int,int);

    void load_D_e(vector<double>&,double,vector<int>,double,vector<double>);

    void initialize_knot_point_par(double&,double&, double&,int&,int&,int&);

    void extr_reduced_vector(vector<double>&,vector<double>&,vector<double>,vector<int>,vector<double>,int, int);

    void calc_residual(vector<double>,vector<double>&,vector<double>,vector<double>,double&,int);

    void set_curr_residual(vector<double>&,int,vector<double>,double&);

    void update_ctrl_pt(vector<double>&,vector<double>,int,vector<int>);

    void extr_reduced_vector_Kgg(vector<double>&,vector<int>,vector<double>,int, int);

//----------------------------------------------------------------------------------------------------------------------------------
//                                                  POST PROCESSING
//----------------------------------------------------------------------------------------------------------------------------------

    void Shape_func_3(vector<int> ,vector<double> , vector<double>& ,vector<double>& ,
                  vector<double> ,vector<double> ,vector<double> ,vector<double> ,
                  vector<double> ,vector<double> ,int ,int ,int ,int ,int ,
                  int ,vector<double>& ,vector<double>& ,vector<double> ,vector<double>& );

    void calc_strain(vector<double>&, vector<double>,vector<double>,int, int, vector<int>,int, vector<double>&);

    void load_postproc_param(vector<double>&,int);

    void calc_displacement(vector<int>,vector<double>,vector<double> ,int ,int ,vector<double>& );

    void calc_stress(vector<double>,vector<double>&,vector<double>);

    void calc_dev_stress(vector<double>, vector<double>&);

    void evaluate_yield_func(vector<double>,double);




