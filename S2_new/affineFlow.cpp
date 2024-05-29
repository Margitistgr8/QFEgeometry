#include "triangle.h"


using namespace std;


#define TWOPI  6.283185307179586
#define Root2 1.4142135623730951
#define Two 2
#define Three 3
#define Debug 0
#define Zero 0.0
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

//Implement size as a command line parameter
int main(int argc, char *argv[])
{
    int N_iter=5000; 
    int L = 4; 
    if(argc ==1)
    printf("In program  %s no input default %d \n",argv[0],L);
    else
    L= atoi(argv[1]);
    printf("Simulation for L=%d\n", L); 

    double z = sqrt((7 + 3 * sqrt(5))/8); 
    //double z = 0; 
    Eigen::Vector3d r1 =  { sqrt(3)/2.0 , -1/2.0,  z };
    Eigen::Vector3d r2 =  {           0,   1.0 ,   z };
    Eigen::Vector3d r3 =  { -sqrt(3)/2.0, -1/2.0,  z };

    r1*=(L/r1.norm());
    r2*=(L/r2.norm());
    r3*=(L/r3.norm());

    double R= r1.norm(); 
    std::vector<Eigen::Vector3d> r={r1,r2,r3}; 

    Lattice lattice = GenerateLattice(L, r);
    std::vector<Triangle> TList = FindTriangleRepresentative(lattice, L); 
    Tensor2D<Eigen::Vector3d> D = FindDerivative(r, lattice); 
    Tensor2D<int> TTList = TriangleVertexMapping(lattice.Vertices, L);

    for (int n=0; n< N_iter; n++)
    {
    double epsilon =1; 
    
    SpMat Op = AreaOperator(lattice, D, TTList, R); 
    Eigen::VectorXd Area = returnCurrentArea(lattice, TTList); 
    double rms = returnRMS(Area); 
    Eigen::VectorXd sol = setUpGMRESsolver(Op, Area, 10);
    printf("iter = %d, rms = %.12f,grad  = %.12f, action = %.12f\n", n, rms, sol.norm(), Area.norm()); 
    //printf("iter = %d\n", n);
    sol*=epsilon;  
    UpdateLattice(lattice, r, sol); 
    if (sol.norm()<=1e-12){break;}
    }
   
  
return 0;
}

