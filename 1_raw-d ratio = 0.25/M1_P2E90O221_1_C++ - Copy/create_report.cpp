#include<iostream>
#include<cmath>
#include<iomanip>
#include <fstream>
#include<string>
#include<vector>
#include<sstream>
#include<iterator>
#include <cstdio>


#include "functions.h"

using namespace std;


//----------------------------------------------------------------------------------------------------------------------------------
//                                                    CALCULATE REPORT - STRESSES
//----------------------------------------------------------------------------------------------------------------------------------

void create_report_stress(int step_load,int tot_no_physpts,vector<double> stress_t1,vector<int> stress_strain_divisor)
{
    stringstream ss;
    //string step_load_filename;
    //string string_1="Step Load ";
    //string string_3=".txt";
    string string_2;
    ss<<step_load;
    ss>>string_2;
    //step_load_filename=string_1+string_2+string_3;

    std::ofstream step_load_file("Stress Results.txt",std::ios::app);
    step_load_file<<"result \"Stresses\" \"IGA\" "<<string_2<<" Matrix OnNodes"<<"\n";
    step_load_file<<"ComponentNames \"Sxx\" \"Syy\" \"Szz\" \"Sxy\" \"Syz\" \"Sxz\" "<<"\n";
    step_load_file<<"ResultRangesTable \"My table\" "<<"\n";
    step_load_file<<"values"<<"\n";
    step_load_file<<"#Node S_xx S_yy S_xz S_xy S_yz S_xz "<<"\n";

    for (short int i=0; i<tot_no_physpts; i++)
    {
        step_load_file<<i+1<<"\t"<<setprecision(8)<<stress_t1[6*i]/stress_strain_divisor[i]<<"\t"<<setprecision(8)
        <<stress_t1[6*i+1]/stress_strain_divisor[i]<<"\t"<<setprecision(8)<<stress_t1[6*i+2]/stress_strain_divisor[i]
        <<"\t"<<setprecision(8)<<stress_t1[6*i+5]/stress_strain_divisor[i]<<"\t"<<setprecision(8)<<stress_t1[6*i+4]/stress_strain_divisor[i]<<"\t"
        <<setprecision(8)<<stress_t1[6*i+3]/stress_strain_divisor[i]<<"\n";

    }
    step_load_file<<"end values"<<"\n";
    step_load_file.close();
}

//----------------------------------------------------------------------------------------------------------------------------------
//                                                    CREATE REPORT - DEFORMATIONS
//----------------------------------------------------------------------------------------------------------------------------------

void create_report_deform(int step_load,int tot_no_physpts,vector<double> disp_t1)
{
    stringstream ss;
    //string step_load_filename;
    //string string_1="Step Load ";
    //string string_3=".txt";
    string string_2;
    ss<<step_load;
    ss>>string_2;
    //step_load_filename=string_1+string_2+string_3;

    std::ofstream step_load_file("Deformation Results.txt",std::ios::app);
    step_load_file<<"result \"Deformations\" \"IGA\" "<<string_2<<" Vector OnNodes"<<"\n";
    step_load_file<<"ComponentNames \"U_xx\" \"U_yy\" \"Uzz\" "<<"\n";
    step_load_file<<"values"<<"\n";

    for (short int i=0; i<tot_no_physpts; i++)
    {
        step_load_file<<i+1<<"\t"<<setprecision(8)<<disp_t1[3*i]<<"\t"<<setprecision(8)<<disp_t1[3*i+1]
        <<"\t"<<setprecision(8)<<disp_t1[3*i+2]<<"\n";

    }
    step_load_file<<"end values"<<"\n";
    step_load_file.close();
}


//----------------------------------------------------------------------------------------------------------------------------------
//                                                    CREATE REPORT - STRESS-STRAIN AT A SPECIFIC PT
//----------------------------------------------------------------------------------------------------------------------------------

void create_report_stress_strain_pt(int step_load,int point,double equivalent_stress,double equivalent_strain)
{
    stringstream ss;
    string string_1,string_2,string_3;
    string filename_1;
    string_1="Stress-Strain at Point ";
    ss<<point;
    ss>>string_2;
    string_3=".txt";
    filename_1=string_1+string_2+string_3;

    std::ofstream stress_strain_pt(filename_1,std::ios::app);
    //stress_strain_pt<<"result \"Deformations\" \"IGA\" "<<string_2<<" Vector OnNodes"<<"\n";
    //stress_strain_pt<<"ComponentNames \"U_xx\" \"U_yy\" \"Uzz\" "<<"\n";
    //stress_strain_pt<<"values"<<"\n";

        stress_strain_pt<<equivalent_strain<<"\t"<<equivalent_stress<<"\n";


    stress_strain_pt.close();
}
