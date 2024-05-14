/**********************************
THESE ARE THE  BASIC TRIANGLE IDENTITIES:

Triangle(a,b,c) edge lenths squared  x = a^2, y = b^2 , z = c^2
I THINK WE SHOULD GO BACK TO a,b,c NOT the Squares! Particularly since
lstar  internal can be negative.

4 A R = \sqrt{x y z}   4 A R = a b c
4 A[x,y,z] = Sqrt[(2 x y + 2 y z + 2 z x - x x - y y - z z)]
R[x,y,z] = Sqrt[x y z]/Sqrt[(2 x y + 2 y z + 2 z x - x x - y y - z z)]
D[R[x,y,z],x] = R^3[x,y,z] ( x^2 - (y-z)^2)/(2 x^2 y z)


***************************************/

#pragma once
#include <cmath>
#include <Eigen/Core> // Core Eigen functionality
#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>
#include <Eigen/Dense>
#include <unordered_map>
#include <string>
#include "util.h"
#define TWOPI  6.283185307179586

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

inline double A(double x, double y, double z)
{
  return sqrt(2*x*y + 2*y*z + 2*z*x - x*x - y*y - z*z)/4;
}

inline double Cr(double x, double y, double z)
{
  //printf("Printing Area in the Cr function: %.12f\n", A(x,y,z));
  return sqrt(x*y*z)/(4*A(x,y,z));
}


inline double lstar(double x, double y, double z, double targetedge)
{
  return 0;
}

inline  double R(double x, double y, double z)
{
  return  sqrt(x*y*z)/sqrt(2*x*y + 2*y*z + 2*z*x - x*x - y*y - z*z);
}

inline double DRx(double x, double y, double z)
{
  return  R(x,y,z) * R(x,y,z) * R(x,y,z) * ( x*x - (y-z)*(y-z))/(2*x*x*y*z);
}

inline double Compute_ls(double a, double b, double c, double edge_sq){
  double l_sq_sum = a+b+c;
  return sqrt(edge_sq)*(l_sq_sum-2.0*edge_sq)/(16*A(a,b,c));
} //Remember to double check this

inline double dotproduct3D(double r1[3], double r2[3])
{
  return r1[0]*r2[0]+r1[1]*r2[1]+r1[2]*r2[2];
}

int returnSerial(int nx, int ny, int N){
    return nx+(N+1)*ny;
}
/*********************************************************
Angle at AB (or xy) vertex:

angAB(x,y,z) = acos( (x + y - z) / (2.0*sqrt(x*y)) )
D[theta,x] =  (-x + y - z)/(2 x Sqrt[4 x y - (x + y - z)^2]) =  (-x + y - z)/(8 x A(x,y,z)) 
D[theta,z] = 1/(4 A[x,y,z])

La* = a(b^2 + c ^2 - a^2)/(8 A[a,b,c])

dualAB[x,y,z]  = la^* la/4 + lb* lb/4 
dualAB[x_, y_, z_] :=  ( -x^2 + 2 x y - y^2 + x z + y z ) /( 32 A[x, y, z])

*************************************************************************/
inline double angAB(double x, double y, double z)
{
  return acos( (x + y - z) / (2.0*sqrt(x*y))) ;
}

inline double DABx(double x, double y, double z)
{
  return  (-x + y - z)/(8* x* A(x,y,z));
}

inline double DABy(double x, double y, double z)
{
  return  (x - y - z)/(8 * y * A(x,y,z));
}

inline double DABz(double x, double y, double z)
{
  return 1.0/( 4.0 * A(x,y,z));
}

inline double DABdualz(double x, double y, double z)
{
  return  x*y*z/( 128 * A(x,y,z) * A(x,y,z)* A(x,y,z));
}

inline double DABdualx(double x, double y, double z)
{
  return  ((x - y)*  (x - y)* (x - y)- 3 * (x - y) * (x + y) *z + (3*x - y)* z * z - z*z*z)/( 512 * A(x,y,z) * A(x,y,z)* A(x,y,z));
}

bool inRange(double low, double high, double x)
{
    return ((x-high)*(x-low) <= 0);
}


double returnRMS(Eigen::VectorXd& vec){
  double mean = vec.array().mean();
  double sq_mean = vec.array().square().mean();
  return sq_mean - mean*mean;
}

template <int L>
void TriangleVertexMapping(int (&arr)[L*L][3]);

template <int L>
int BuildBasisVectorInfo(int (&arr)[L+1][L+1][3]);

template <int L>
void BuildAreaOperator(std::vector<T>& (tripleList),int (&map)[L*L][3],int (&vecInf)[L+1][L+1][3],double (&rvec)[L+1][L+1][3],double(&xvec)[L+1][L+1][3],double r12[3], double r13[3], double r23[3]);

template <int L>
Eigen::VectorXd CurrentAreaList(int (&arr)[L*L][3], double (&rvec)[L+1][L+1][3]);

template <int L>
Eigen::VectorXd CurrentPerimeterList(int (&arr)[L*L][3], double (&rvec)[L+1][L+1][3]);


template <int L>
Eigen::VectorXd  CurrentCircumradiusList(int (&arr)[L*L][3], double (&rvec)[L+1][L+1][3]);

template <int L>
void SetPosition(double (&xivec)[L+1][L+1][3],double (&rvec)[L+1][L+1][3], double r1[3], double r2[3],double r3[3],double (&xvec)[L+1][L+1][3]);

template <int L>
void BuildPerimeterOperator(std::vector<T>& (tripletList),int (&map)[L*L][3],int (&vecInf)[L+1][L+1][3],double(&rvec)[L+1][L+1][3], double (&xvec)[L+1][L+1][3],double r12[3], double r13[3], double r23[3]);

template <int L>
void BuildCircumradiusOperator(std::vector<T>& (tripletList),int (&map)[L*L][3],int (&vecInf)[L+1][L+1][3],double(&rvec)[L+1][L+1][3],double (&xvec)[L+1][L+1][3],double r12[3], double r13[3], double r23[3]);


template <int L>
int BuildEquivalenceClass(int (&orbitclass)[L+1][L+1][2]);

template <int L>
void updateBarycentricwSymmetry(double (&xivec)[L+1][L+1][3], int (&arr)[L+1][L+1][3] , int (&classes)[L+1][L+1][2], Eigen::VectorXd& sol, int dof);

template<int L> //The orbit class object is initialized to -1 
void UpdateOrbit(int (&orbitclass)[L+1][L+1][2],double  (&xivec)[L+1][L+1][3], int nx, int ny);

template <int L>
bool CheckSymmetry(double (&xivec)[L+1][L+1][3],int (&vecInf)[L+1][L+1][3], int (&orbitclass)[L+1][L+1][2], int dof);

template <int L>
void setUpGMRESsolver(std::vector<T>& (tripletList),  Eigen::VectorXd& b, Eigen::VectorXd& sol,int dof, int iter_num);

template <int L>
void setUpStackedGMRESsolver(std::vector<T>& (tripletList),std::vector<T>& (tripletList2), Eigen::VectorXd& b, Eigen::VectorXd& sol,int dof, int iter_num);

template <int L>
void PrintGeometry(int (&arr)[L*L][3], double (&rvec)[L+1][L+1][3]);

template <int L>
void SubtractAverage(std::vector<T>& tripletList, int dof);
template <int L>
void PrintRenormalizedRMS(int (&arr)[L*L][3], double (&rvec)[L+1][L+1][3]);

template<int L>
Eigen::VectorXd returnGradient(std::vector<T>& tripletList, Eigen::VectorXd& target, int dof); 

template <int L>
void EquivalenceClassID(std::unordered_map<std::string, int> &orbit_map, int (&orbitclass)[L+1][L+1][2], int dof); 


template <int L>
void SaveBarycentric(FILE* file, double (&xivec)[L+1][L+1][3], int (&orbitclass)[L+1][L+1][2], int classnumber);

template <int L>
void SaveMetric(FILE* file, double (&rvec)[L+1][L+1][3], int (&map)[L*L][3]);

template <int L>
void CurrentDualArea(int (&arr)[L*L][3], double (&rvec)[L+1][L+1][3], double (&dualArea)[L+1][L+1]); 

template<int L>
void PrintDualInfo(int (&arr)[L*L][3], double (&rvec)[L+1][L+1][3]);

template <int L>
void ReadRData(const std::string& filePath, double (&rvec)[L+1][L+1][3]);

template <int L>
void ReadBarycentric(const std::string& filePath, double (&xivec)[L+1][L+1][3], int (&orbitclass)[L+1][L+1][2], int classnumber);


template <int L>
void CurrentDeficitAngle(int (&arr)[L*L][3], double (&rvec)[L+1][L+1][3], double (&defangle)[L+1][L+1]);

template<int L>
void PrintDeficitInfo(int (&arr)[L*L][3], double (&rvec)[L+1][L+1][3]);

template <int L>
void SavePosition(FILE* file, double (&rvec)[L+1][L+1][3]);

template <int L>
void SaveVertexFunction(FILE* file, double (&func)[L+1][L+1]);

template <int L>
void SaveFaceFunction(Eigen::VectorXd& quantity, FILE* file, int (&map)[L*L][3]);

template <int L>
void SaveCurvedFaceFunction(Eigen::VectorXd& quantity, FILE* file, int (&map)[L*L][3], double (&rvec)[L+1][L+1][3]);

template<int L>
void ReadOrbitData(const std::string& filePath, int (&orbitclass)[L+1][L+1][2], int classnumber, double (&target)[L+1][L+1]);
#include "traingle.tpp"
