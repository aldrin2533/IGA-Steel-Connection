#include<string>
#include<iostream>
#include<vector>
#include<fstream>
#include<cmath>

#include "functions.h"


using namespace std;

//----------------------------------------------------------------------------------------------------------------------------------
//                                         MODEL'S NURBS BASIS FUNCTIONS AND PROPERTIES
//----------------------------------------------------------------------------------------------------------------------------------

//  ********************************************************
//  FIND THE SPAN OF A KNOT POINT
//
//  This function returns the knot index, based on
//  C++ indexing, where the value of xi is contained
//  ********************************************************

int findspan(int noofspandindex,int order,double xi,vector<double> Xi)
{

if(xi==Xi[noofspandindex+1]){return (noofspandindex+1);}
int low=order;
int high=noofspandindex+1;
int mid=0.5*(low+high);
while(xi<Xi[mid]||xi>=Xi[mid+1])
{
    if(xi<Xi[mid]){high=mid;}
    else{low=mid;}
    mid=0.5*(low+high);
}
return mid;
}

//  ********************************************************
//  COMPUTE ORDER OF B-SPLINE BASIS FUNCTION AND NUMBER OF
//  BASIS FUNCTIONS FOR A GIVEN KNOT VECTOR
//  ********************************************************

void calc_bspline_prop(vector<double> vector_i,int& order,int& no_basis_func){
int counter=0;
    for(int ii=0;ii<vector_i.size();ii++){
        if(vector_i[ii]==0){
            counter++;}
        else{break;}
    }
order=counter-1;
no_basis_func=vector_i.size()-order-1;
}

//  ********************************************************
//  COMPUTE NURBS MODEL PROPERTIES
//  ********************************************************
void calc_nurbs_mod_prop(int p,int q,int n,int m,int& nel,int& nnp,int& nen,int& noctrlpt, int& noxi,int& noeta,
                         int& nophyspts,vector<double> Xi,vector<double> Eta)
{
    double counter_1=0;
        for(int i=0;i<Xi.size()-1;i++)
            {
            if(Xi[i]==Xi[i+1]){counter_1=counter_1;}
            else{counter_1++;}
            }

        double counter_2=0;
        for(int i=0;i<Eta.size()-1;i++)
            {
            if(Eta[i]==Eta[i+1]){counter_2=counter_2;}
            else{counter_2++;}
            }

    nel=counter_1*counter_2;
    //nel=(n-p)*(m-q)*(l-r);
    nnp=n*m;
    nen=(p+1)*(q+1);
    noctrlpt=n*m;
    noxi=n-p+1;
    noeta=m-q+1;
    nophyspts=(n-p+1)*(m-q+1);
}

//----------------------------------------------------------------------------------------------------------------------------------
//                                                      BUILD NURBS COORDINATE ARRAY
//----------------------------------------------------------------------------------------------------------------------------------

void build_NURBS_coords(vector<double> Xi,vector<double> Eta,vector<int>& NURBS_coords)
{
    int element = 0;

    vector<double> vector_xi;
    vector<double> vector_eta;

    for(int i=0;i<Xi.size()-1;i++)
    {
        if(Xi[i]!=Xi[i+1]){vector_xi.push_back(i);}
    }


    //cout<<vector_xi.size()<<endl;

    for(int i=0;i<Eta.size()-1;i++)
    {
        if(Eta[i]!=Eta[i+1]){vector_eta.push_back(i);}
    }

    NURBS_coords.resize(vector_xi.size()*vector_eta.size()*2);

        for(int j=0;j<vector_eta.size();j++)
        {
            for(int i=0;i<vector_xi.size();i++)
            {
                element++;
                NURBS_coords[2*(element-1)]=vector_xi[i];
                NURBS_coords[2*(element-1)+1]=vector_eta[j];
                //NURBS_coords[3*(element-1)+2]=vector_zeta[k];
                //cout<<NURBS_coords[3*(element-1)]<<" "<<NURBS_coords[3*(element-1)+1]<<" "<<NURBS_coords[3*(element-1)+2]<<endl;
            }
        }
    vector_xi.clear();
    vector_eta.clear();

}

//  ********************************************************
//  COMPUTE NUMBER OF RESTRAINED AND UNRESTRAINED DOF
//
//      Returns number of restrained and unrestrained DOF.
//      Stores unrestrained DOF to a vector.
//  ********************************************************

void det_dof_type(int DOF,int noctrlpt, vector<int>& unrestrained_dof_data,int& no_unrestrained_dof,int& no_restrained_dof,
                              vector<int> bound_condition)
{
    no_unrestrained_dof=0;
    unrestrained_dof_data.clear();

        for(int j=0;j<noctrlpt;j++){
        for(int i=0;i<DOF;i++){
          if(bound_condition[DOF*j+i]==0){
            no_unrestrained_dof++;
            unrestrained_dof_data.push_back(j+1);
            unrestrained_dof_data.push_back(i);
        }}}
        no_restrained_dof=DOF*noctrlpt-no_unrestrained_dof;
}


double calc_parametric_coord(vector<double> vector_knot,int knot_index,double gauss_pt)
{
    double f;
    f=vector_knot[knot_index]+0.5*(gauss_pt+1)*(vector_knot[knot_index+1]-vector_knot[knot_index]);
    return f;
}

//  ********************************************************
//  COMPUTE B-SPLINE BASIS FUNCTION
//  ********************************************************

double calc_basis_func(int index,double knotpointvalue,vector<double> knotvector,int noofbasisfunctions,int orderofbasisfunctions)
{
double f,saved,Xleft,Xright,temp;
double g[orderofbasisfunctions+1];

if((index==0&&knotpointvalue==knotvector[0])||(index==noofbasisfunctions-orderofbasisfunctions-1&&knotpointvalue==knotvector[noofbasisfunctions])) {f=1;}

if(knotpointvalue==1){f=1;}

for(int j=0; j<=orderofbasisfunctions;j++)
    {
    if(knotpointvalue>=knotvector[index+j]&&knotpointvalue<knotvector[index+j+1]){g[j]=1;}
    else{g[j]=0;}
    }

for(int k=1;k<=orderofbasisfunctions;k++)
    {
    if(g[0]==0){saved=0;}
    else{saved=((knotpointvalue-knotvector[index])*g[0])/(knotvector[index+k]-knotvector[index]);}

    for(int j=0;j<orderofbasisfunctions-k+1;j++)
    {Xleft=knotvector[index+j+1]; Xright=knotvector[index+j+k+1];

        if(g[j+1]==0){g[j]=saved;saved=0;}
        else
        {
            temp=g[j+1]/(Xright-Xleft);
            g[j]=saved+(Xright-knotpointvalue)*temp;
            saved=(knotpointvalue-Xleft)*temp;
        }}}
f=g[0];
return f;
}

// ********************************************************
//  COMPUTE B-SPLINE FUNCTION DERIVATIVE
// ********************************************************

double calc_Bspline_deriv(int index,double knotpointvalue,std::vector<double> knotvector,int noofbasisfunctions,int orderofbasisfunctions,
                          int i)
{
    double temp1=0,temp2=0;
            if(knotvector[index+i]==knotvector[index-orderofbasisfunctions+i]){temp1=0;}
            else{temp1=1/(knotvector[index+i]-knotvector[index-orderofbasisfunctions+i]);}
            if(knotvector[index+1+i]==knotvector[index-orderofbasisfunctions+1+i]){temp2=0;}
            else{temp2=1/(knotvector[index+1+i]-knotvector[index-orderofbasisfunctions+1+i]);}
    double N_1=0,N_2=0;
    N_1=calc_basis_func(index-orderofbasisfunctions+i,knotpointvalue,knotvector,noofbasisfunctions,orderofbasisfunctions-1);
    N_2=calc_basis_func(index-orderofbasisfunctions+i+1,knotpointvalue,knotvector,noofbasisfunctions,orderofbasisfunctions-1);

    double N_deriv=0;
    N_deriv =(orderofbasisfunctions)*N_1*temp1-orderofbasisfunctions*N_2*temp2;

    return N_deriv;

}

// ********************************************************
//  STORE B-SPLINE BASIS FUNCTIONS AND ITS' DERIVATIVE
// ********************************************************

void form_B_splines(vector<double>& B_spline,vector<double>& B_spline_deriv,int knot_pt_index,int order,double knot_pt_val,vector<double> knot_vect,
     int no_basis_funcs)
{
    B_spline.clear();
    B_spline_deriv.clear();
        for(int i=0;i<order+1;i++){

            if(knot_pt_val==1&&(knot_pt_index-order+i)==knot_pt_index){B_spline.push_back(1);}
            if((knot_pt_val==1)&&((knot_pt_index-order+i)<knot_pt_index)){B_spline.push_back(0);}
            else{
            B_spline.push_back(calc_basis_func(knot_pt_index-order+i,knot_pt_val,knot_vect,no_basis_funcs,order));
            B_spline_deriv.push_back(calc_Bspline_deriv(knot_pt_index,knot_pt_val,knot_vect,no_basis_funcs,order,i));}
        }

}

// ********************************************************
//  COMPUTE NURBS BASIS FUNCTION AND ITS DERIVATIVE W.R.T.
//  PHYSICAL COORDINATES
//                                  ** verified 04052021 **
// ********************************************************

void Shape_func_1(vector<int> ctrlpt_incidence,vector<double> ctrlpt_coordinates, vector<double>& NURBS_bas_func,vector<double>& NURBS_deriv,
                  vector<double> basis_func_N,vector<double> basis_func_M,vector<double> basis_deriv_N,
                  vector<double> basis_deriv_M,int order_p,int order_q,int element,int nen,vector<double>& dR_dx,vector<double>& dx_dxi)
{
    int loc_num=0;
    double sumtot=0;
    double sumxi=0;
    double sumeta=0;
    double sumzeta=0;
    int e;
    int cc;
    int deg_of_freedom=2;

    NURBS_bas_func.clear();
    NURBS_deriv.clear();

            // Determine NURBS basis function and its derivative w.r.t parametric coordinate
            for(int j=0;j<order_q+1;j++){
            for(int i=0;i<order_p+1;i++){

                loc_num+=1;
                e = element-1;

                cc=ctrlpt_incidence[nen*e+(loc_num-1)];
                //cout<<cc<<endl;

                NURBS_bas_func.push_back(basis_func_N[i]*basis_func_M[j]*ctrlpt_coordinates[4*(cc-1)+3]);
                sumtot+=NURBS_bas_func[loc_num-1];

                NURBS_deriv.push_back(basis_deriv_N[i]*basis_func_M[j]*ctrlpt_coordinates[4*(cc-1)+3]);
                sumxi+=NURBS_deriv[(loc_num-1)*deg_of_freedom];

                NURBS_deriv.push_back(basis_func_N[i]*basis_deriv_M[j]*ctrlpt_coordinates[4*(cc-1)+3]);
                sumeta+=NURBS_deriv[(loc_num-1)*deg_of_freedom+1];

                }}

            for(int a=0;a<loc_num;a++){
                NURBS_bas_func[a]=NURBS_bas_func[a]/sumtot;
                NURBS_deriv[deg_of_freedom*a]=(NURBS_deriv[deg_of_freedom*a]*sumtot-NURBS_bas_func[a]*sumxi)/pow(sumtot,2);
                NURBS_deriv[deg_of_freedom*a+1]=(NURBS_deriv[deg_of_freedom*a+1]*sumtot-NURBS_bas_func[a]*sumeta)/pow(sumtot,2);
            }
        //cout<<"ELEMENT:"<<e+1<<endl;
        //for(int i=0;i<nen;i++)
        //{
          //  cout<<NURBS_bas_func[i]<<" "<<NURBS_deriv[3*i]<<" "<<NURBS_deriv[3*i+1]<<" "<<NURBS_deriv[3*i+2]<<endl;
        //}

        loc_num=0;

            // Determine mapping from parameter space to physical space
            dx_dxi.clear();
            dx_dxi.resize(deg_of_freedom*deg_of_freedom);

            for(int j=0;j<order_q+1;j++){
            for(int i=0;i<order_p+1;i++){

                loc_num+=1;
                for(int aa=0;aa<deg_of_freedom;aa++){
                for(int bb=0;bb<deg_of_freedom;bb++){
                        int cc=ctrlpt_incidence[nen*e+(loc_num-1)];
                        //cout<<cc<<endl;
                        //cout<<ctrlpt_coordinates[4*(cc-1)+aa]<<endl;
                        dx_dxi[deg_of_freedom*aa+bb]+=(ctrlpt_coordinates[4*(cc-1)+aa]*NURBS_deriv[deg_of_freedom*(loc_num-1)+bb]);
                        //cout<<bb<<endl;
                        //printf("cc, aa, bb : %d %d %d %f \n",cc,aa,bb,ctrlpt_coordinates[4*(cc-1)+aa]);
                        //printf("dx_dxi , aa , bb : %f %d %d \n",dx_dxi,aa,bb);


                }}
            }}

            //cout<<dx_dxi[0]<<" "<<dx_dxi[1]<<endl;
            //cout<<dx_dxi[2]<<" "<<dx_dxi[3]<<endl;

            //for(int i=0;i<9;i++){cout<<dx_dxi[i]<<endl;}

            // Determine inverse of mapping from parameter space to physical space
            vector<double> dxi_dx;
            dxi_dx.clear();
            dxi_dx.resize(deg_of_freedom*deg_of_freedom);

            double determinant=0;
            determinant =dx_dxi[0]*dx_dxi[3]-dx_dxi[1]*dx_dxi[2];

                dxi_dx[0]=dx_dxi[3]/determinant;
                dxi_dx[1]=-dx_dxi[1]/determinant;
                dxi_dx[2]=-dx_dxi[2]/determinant;
                dxi_dx[3]=dx_dxi[0]/determinant;

            //  Determine derivative of NURBS basis functions w.r.t physical coordinate
            dR_dx.clear();
            dR_dx.resize(nen*deg_of_freedom);

                for(int loc_num=0;loc_num<nen;loc_num++){
                for(int aa=0;aa<deg_of_freedom;aa++){
                for(int bb=0;bb<deg_of_freedom;bb++){
                    dR_dx[deg_of_freedom*loc_num+aa]+=(NURBS_deriv[deg_of_freedom*loc_num+bb]*dxi_dx[deg_of_freedom*bb+aa]);
                }}}


}

// ********************************************************
//  COMPUTE FOR THE DETERMINANT OF THE JACOBIAN MATRIX
//  OR MAPPING FROM PARENT TO PHYSICAL COORDINATES
// ********************************************************
double Shape_function_2(int xi_index,int eta_index,vector<double> Knot_Vector_xi,
                      vector<double> Knot_Vector_eta,vector<double> dx_dxi)
{
    //  Calculate mapping from Parent Element to Parameter Space.

    int deg_of_freedom=2;
    vector<double> dxi_dtildexi;
    dxi_dtildexi.clear();

        dxi_dtildexi.push_back(0.5*(Knot_Vector_xi[xi_index+1]-Knot_Vector_xi[xi_index]));
        dxi_dtildexi.push_back(0.5*(Knot_Vector_eta[eta_index+1]-Knot_Vector_eta[eta_index]));

    // Calculate Jacobian Matrix
    vector<double> dx_dtildexi;
    dx_dtildexi.clear();
    dx_dtildexi.resize(deg_of_freedom*deg_of_freedom);

    //for (int aa=0;aa<deg_of_freedom;aa++){
    //for(int bb=0;bb<deg_of_freedom;bb++){
            //dx_dtildexi[deg_of_freedom*aa+bb]+=dx_dxi[deg_of_freedom*aa+bb]*dxi_dtildexi[bb];
    //}}

            dx_dtildexi[0]=dx_dxi[0]*dxi_dtildexi[0];
            dx_dtildexi[1]=dx_dxi[1]*dxi_dtildexi[1];
            dx_dtildexi[2]=dx_dxi[2]*dxi_dtildexi[0];
            dx_dtildexi[3]=dx_dxi[3]*dxi_dtildexi[1];
            //cout<<dx_dtildexi[0]<<" "<<dx_dtildexi[1]<<endl;
            //cout<<dx_dtildexi[2]<<" "<<dx_dtildexi[3]<<endl;

    //  Compute Determinant of Jacobian Matrix
    double determinant=0;
    determinant=dx_dtildexi[0]*dx_dtildexi[3]-dx_dtildexi[1]*dx_dtildexi[2];
return determinant;
}



//----------------------------------------------------------------------------------------------------------------------------------







