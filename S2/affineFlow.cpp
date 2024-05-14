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


#define L 96
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




    FILE* fptrFull = NULL;  // C style
  #if 0
  char in_nameFull[64];
  sprintf(in_nameFull,"Jin_data/CurvedTriangleupdate_%d_initial.dat",L); // filename
  fptrFull = fopen(in_nameFull, "w"); 
  SavePosition<L>(fptrFull, rvec); 
  fclose(fptrFull);
  #endif

    FILE* pos = NULL;
    char file[64];
    sprintf(file,"Area(new)/data/xivec_%d.dat",L); 
    //sprintf(file,"resultsS2new/orbit/q5k%d_orbit.dat",L); 
    //pos = fopen(file,"w");
    //SavePosition<L>(pos, rvec); 
    //printf("%d ", L); 
    //fclose(pos);

  
    //ReadRData<L>(file, rvec);

    //ReadBarycentric<L>(file, xivec, orbitclass, eqclass); 
    //SetPosition<L>(xivec, rvec, r1_target, r2_target,r3_target, xvec);


    //fclose(pos);
    #if 0
    char ct[64]; 
    sprintf(ct, "ct/ct_5_%dnoref.dat", L); 
    double counterterm[L+1][L+1] = {}; 
    ReadOrbitData<L>(ct, orbitclass, eqclass, counterterm);

    pos = fopen("counterterm.dat","w");
    SaveVertexFunction<L>(pos, counterterm); 
    fclose(pos);
    #endif


#if 0
  char out_nameFull[64];
  sprintf(out_nameFull,"P(new)/data/CurvedTriangleupdate_%d_final.dat",L); // filename
  fptrFull = fopen(out_nameFull, "w"); 
  SavePosition<L>(fptrFull, rvec); 
  fclose(fptrFull);
#endif 

    double dualAr[L+1][L+1] = {}; 
    double defAng[L+1][L+1] = {}; 
    CurrentDualArea<L>(map, rvec, dualAr); 
    CurrentDeficitAngle<L>(map, rvec, defAng); 
    Eigen::VectorXd avec = CurrentAreaList<L>(map, rvec); 
    Eigen::VectorXd pvec = CurrentPerimeterList<L>(map, rvec); 
    Eigen::VectorXd crvec = CurrentCircumradiusList<L>(map, rvec); 

    // for (int ny = 0; ny<=L; ny++){
    //     for (int nx = 0; nx<=L-ny; nx++){
    //       printf("%.4f %.4f %.12f %.12f %.12f\n", nx + ny/2.0,   sqrt(3.0)*ny/2.0, rvec[ny][nx][0],rvec[ny][nx][1], rvec[ny][nx][2]); 
    // }
    // }



    #if 0
    pos = fopen("dualArr.dat","w");
    SaveVertexFunction<L>(pos, dualAr); 
    fclose(pos);

    pos = fopen("deficitAngle.dat","w");
    SaveVertexFunction<L>(pos, defAng); 
    fclose(pos);
    
    pos = fopen("Area.dat","w");
    SaveFaceFunction<L>(avec, pos, map); 
    fclose(pos);
    
    pos = fopen("Perimeter.dat","w");
    SaveFaceFunction<L>(pvec, pos, map); 
    fclose(pos);

    pos = fopen("Circumradius.dat","w");
    SaveFaceFunction<L>(crvec, pos, map); 
    //SaveCurvedFaceFunction<L>(crvec, pos, map, rvec); 
    fclose(pos);
    #endif

    #if 1
    printf("%d ", L); 
    PrintGeometry<L>(map, rvec); 
    //std::cout<<avec<<endl; 
    PrintDualInfo<L>(map, rvec);
    PrintDeficitInfo<L>(map, rvec); 
    printf("\n"); 
    #endif

#if 1
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
    BuildCircumradiusOperator<L>(tripletList, map, vecInf, rvec, xvec, r12, r13, r23);
    SubtractAverage<L>(tripletList, dof); 
    Eigen::VectorXd avec = CurrentCircumradiusList<L>(map, rvec);
    avec = avec.array()-avec.mean(); 

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
        Eigen::VectorXd temp= CurrentCircumradiusList<L>(map, rvec); 
        temp = temp.array()-temp.mean(); 
        if (action-temp.dot(temp)>=AG_condition*epsilon){           
        epsilon_max = epsilon;
        condition=1;}
        else{
          counter+=1; 
        }
    }    

    #if 0
    char file[64];
    sprintf(file,"temparea/L_%diter_%d.dat",L, iter); 
    pos = fopen(file,"w");
    avec=CurrentAreaList<L>(map, rvec); 
    SaveFaceFunction<L>(avec, pos, map); 
    fclose(pos);
    #endif

     #if 0
    std::vector<T> tripletList;
    std::vector<T> tripletListB;
    BuildPerimeterOperator<L>(tripletList, map, vecInf, rvec, xvec, r12, r13, r23);
    BuildCircumradiusOperator<L>(tripletListB, map, vecInf, rvec, xvec, r12, r13, r23);


    Eigen::VectorXd avec = CurrentPerimeterList<L>(map, rvec);
    Eigen::VectorXd bvec = CurrentCircumradiusList<L>(map, rvec);

    Eigen::VectorXd target(2*L*L); 
    target<<avec, bvec; 

    Eigen::VectorXd sol(dof);
    setUpStackedGMRESsolver<L>(tripletList,tripletListB, target, sol, dof, 1);
    Eigen::VectorXd Gradient = returnGradient<L>(tripletList, avec, dof)+returnGradient<L>(tripletListB,bvec, dof);

    //printf("Gradient^2: %.12f\n", Gradient.dot(Gradient));

    double epsilon_max= 1e-1; //Maximum candidate step size
    double AG_condition= -sol.dot(Gradient)/2.0;
    double ratio = 0.8;
    double action = target.dot(target); 
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
        //cout<< testupdate<<endl; 
        updateBarycentricwSymmetry<L>(xivecCopy, vecInf,orbitclass, testupdate, dof);
        SetPosition<L>(xivecCopy, rvec, r1,r2,r3,xvec);
        Eigen::VectorXd temp(2*L*L); 
        temp<< CurrentPerimeterList<L>(map, rvec), CurrentCircumradiusList<L>(map, rvec); 
        if (action-temp.dot(temp)>=AG_condition*epsilon){           
        epsilon_max = epsilon;
      condition=1;}
        else{
          counter+=1; 
        }
    }    
     #endif
    
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

