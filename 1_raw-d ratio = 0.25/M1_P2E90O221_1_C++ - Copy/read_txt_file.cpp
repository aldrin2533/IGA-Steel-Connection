#include<string>
#include<iostream>
#include<vector>
#include<fstream>

#include "functions.h"

using namespace std;

//----------------------------------------------------------------------------------------------------------------------------------
//                                                         READ TEXT FILE
//----------------------------------------------------------------------------------------------------------------------------------

void text_file::input_string_file_name(){
    KV_xi_string="Knot Vector Xi.txt";
    KV_eta_string="Knot Vector Eta.txt";
    KV_zeta_string="Knot Vector Zeta.txt";
    ctrl_pt_coords_string="Control Points Coordinates.txt";
    ctrl_pt_inc_string="Control Point Incidences.txt";
    bound_cond_string="Boundary Conditions.txt";
    phys_pt_coords_string="Nodes Coordinates.txt";
    patch_inc_string="Patch Element Member Incidence.txt";
}

int text_file::determine_size(string file_name){
    int counter=0;
    std::ifstream ifile(file_name, std::ios::in);
        if (!ifile.is_open()){
            std::cerr << "There was a problem opening the input file!\n";
            exit(1);
            }
        double num=0;
        while (ifile >> num){
            counter+=1;
            }
    return counter;
}

double text_file::read_textfile_double(string file_name,int ii){
    double value_i;
    vector<double> vector_i;
    std::ifstream ifile(file_name, std::ios::in);
        if (!ifile.is_open()){
            std::cerr << "There was a problem opening the input file!\n";
            exit(1);
            }
        double num=0;
        while (ifile >> num){
            vector_i.push_back(num);
            }
        value_i=vector_i[ii];
    return value_i;
}

int text_file::read_textfile_int(string file_name,int ii){
    int value_i;
    vector<int> vector_i;
    std::ifstream ifile(file_name, std::ios::in);
        if (!ifile.is_open()){
            std::cerr << "There was a problem opening the input file!\n";
            exit(1);
            }
        int num=0;
        while (ifile >> num){
            vector_i.push_back(num);
            }
        value_i=vector_i[ii];
    return value_i;
}




//----------------------------------------------------------------------------------------------------------------------------------








