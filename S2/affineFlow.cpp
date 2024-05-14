#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <array>
#include <random>
#include <fstream>
#include <climits>
#include <Eigen/Core> // Core Eigen functionality
#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>
#include <Eigen/Dense>
#include "triangle.h"


using namespace std;


#define L 32
#define TWOPI  6.283185307179586
#define Root2 1.4142135623730951
#define Two 2
#define Three 3
#define Debug 0
#define Zero 0.0
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;


int main(int argc, char *argv[])
{
  double xvec[L+1][L+1][Three]; //  vectors on plane at fixed z
  double rvec[L+1][L+1][Three];  // vector on sphere
  double xivec[L+1][L+1][Three]; // Barycentric coordinates
  int    N_iter=200000; 
  //int N_iter=50;
    double z = sqrt((7 + 3 * sqrt(5))/8); 
    //double z = 0.0;
    double r1_target[3] =  { sqrt(3)/2.0 , -1/2.0,  z };
    double r2_target[3] =  {           0,   1.0 ,   z };
    double r3_target[3] =  { -sqrt(3)/2.0, -1/2.0,  z };


  for(int mu =0; mu < 3; mu++)  // setting target r vectors to unit norm 
    {
      r1_target[mu] = L*r1_target[mu]/sqrt(1 + z*z);
      r2_target[mu] = L*r2_target[mu]/sqrt(1 + z*z);
      r3_target[mu] = L*r3_target[mu]/sqrt(1 + z*z);
    }

  

  for(int ny = 0; ny<= L; ny++) //initialize Barycentric Coordinates
    {
      for(int nx = 0; nx <= L; nx++)
	{   
    int nz = L - ny-nx; 
	  xivec[ny][nx][0] = double(nx)/double(L); 
    xivec[ny][nx][1] = double(ny)/double(L); 
    xivec[ny][nx][2] = double(nz)/double(L); 
	}
    }

    SetPosition<L>(xivec, rvec, r1_target, r2_target,r3_target, xvec);
  
    int map[L*L][3] = {}; 
    int vecInf[L+1][L+1][3];
    int orbitclass[L+1][L+1][2];

    for (int i=0; i<=L; i++){
        for (int j=0; j<=L; j++){
            vecInf[i][j][0]=-1; 
            vecInf[i][j][1]=-1; 
            vecInf[i][j][2]=-1; 
            orbitclass[i][j][0] = -1; 
            orbitclass[i][j][0] = -1; 
        }
    }
    TriangleVertexMapping<L>(map);
    int dof =BuildBasisVectorInfo<L>(vecInf);  
    int eqclass = BuildEquivalenceClass<L>(orbitclass);
    //printf("dof: %d\n", dof); 
    std::unordered_map<std::string, int> orbit_map;
    EquivalenceClassID<L>(orbit_map, orbitclass,  dof);

#if 1 //Optimization Routine
double r1[3] =  { r1_target[0] , r1_target[1],  r1_target[2] };
double r2[3] =  { r2_target[0] , r2_target[1],  r2_target[2] };
double r3[3] =  { r3_target[0] , r3_target[1],  r3_target[2] };

SetPosition<L>(xivec, rvec, r1, r2,r3, xvec); //Set x and r coordinates based on barycentric coordinates

double r12[3] = {}; 
double r13[3] = {}; 
double r23[3] = {}; 

for (int mu = 0; mu<3; mu++){
        r12[mu] = r1[mu]-r2[mu];
        r23[mu] = r2[mu]-r3[mu];
        r13[mu] = r1[mu]-r3[mu];
        }
////////////////////Start of Modification


double res = 1; 
int counter = 0; 

printf("0 %.12f %.16f", Zero, Zero); 
PrintGeometry<L>(map, rvec); 
PrintRenormalizedRMS<L>(map, rvec); 
PrintDualInfo<L>(map, rvec); 
printf("\n"); 

for (int iter = 0; iter<N_iter; iter++){

    std::vector<T> tripletList;
    BuildAreaOperator<L>(tripletList, map, vecInf, rvec, xvec, r12, r13, r23);
    Eigen::VectorXd avec = CurrentAreaList<L>(map, rvec);


    Eigen::VectorXd sol(dof);
    setUpGMRESsolver<L>(tripletList, avec, sol, dof, 1);
    Eigen::VectorXd Gradient = returnGradient<L>(tripletList, avec, dof);


    double epsilon_max= 1; //Maximum candidate step size
    double AG_condition= -sol.dot(Gradient)/2.0;
    double ratio = 0.8;
    double action = avec.dot(avec); 
    int condition = 0;   
    int counter = 0; 
    while (condition!=1){
        double xivecCopy[L+1][L+1][3] = {};
        for (int nx = 0; nx<=L; nx++){
          for (int ny = 0; ny<=L-nx; ny++){
            for (int mu = 0; mu<3; mu++){
              xivecCopy[ny][nx][mu]=xivec[ny][nx][mu]; 
            }
          }
        }
        double epsilon = epsilon_max*pow(ratio, counter); 
        Eigen::VectorXd testupdate = sol*epsilon; 

        updateBarycentricwSymmetry<L>(xivecCopy, vecInf,orbitclass, testupdate, dof);
        SetPosition<L>(xivecCopy, rvec, r1,r2,r3,xvec);
        Eigen::VectorXd temp= CurrentAreaList<L>(map, rvec); 
        if (action-temp.dot(temp)>=AG_condition*epsilon){           
        epsilon_max = epsilon;
        condition=1;}
        else{
          counter+=1; 
        }
    }    
    res = sol.dot(sol); 
    sol*=epsilon_max; 
    updateBarycentricwSymmetry<L>(xivec, vecInf,orbitclass, sol, dof);
    SetPosition<L>(xivec, rvec, r1,r2,r3,xvec);
    counter = 0; 
    printf("%d %.12f %.16f ", iter+1, epsilon_max, res); 
    PrintGeometry<L>(map, rvec); 
    PrintRenormalizedRMS<L>(map, rvec); 
    PrintDualInfo<L>(map, rvec); 
    printf("\n");  
    if(res<=1e-16 || sqrt(res*epsilon_max)<=1e-16){//printf("ending at iter = %d\n", iter);
    break;}
}


  FILE* testFull = NULL;  // C style
  char test[64];
  sprintf(test,"data/xivec_%d.dat",L); // filename
  testFull = fopen(test,"w");
  
  SaveBarycentric<L>(testFull, xivec, orbitclass, eqclass);
  fclose(testFull);

#if 0
  FILE* finptrFull = NULL;  // C style
  char final[64];
  sprintf(final,"data/CurvedTriangleupdate_%d_final.dat",L); // filename
  SavePosition<L>(finptrFull, rvec); 
  fclose(finptrFull);
#endif 
#endif

////////////////////End of Modification
return 0;
}

