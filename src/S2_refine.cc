#include "S2.h"


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
    int stepsize= 1; 
    int q = 5; 
    std::string data_dir = "";

    if(argc ==1)
    printf("In program  %s no input default %d \n", argv[0], L);
    else if (argc<5)
    {
    L= atoi(argv[1]);
    N_iter = atoi(argv[2]);
    stepsize = atoi(argv[3]);
    data_dir = argv[4];
    }
    else
    {
    L= atoi(argv[1]);
    N_iter = atoi(argv[2]);
    stepsize = atoi(argv[3]);
    data_dir = argv[4];   
    q = atoi(argv[5]);
    }
    printf("Simulation for L=%d with q=%d \n", L, q); 
    double z = z = sqrt((7 + 3 * sqrt(5))/8);
    if (q==4){z = 1/sqrt(2); }
    if (q==3){z = 1/(2*sqrt(2)); }
    
    printf("z=%.12f\n", z);
    Eigen::Vector3d r1 =  { sqrt(3)/2.0 , -1/2.0,  z };
    Eigen::Vector3d r2 =  {           0,   1.0 ,   z };
    Eigen::Vector3d r3 =  { -sqrt(3)/2.0, -1/2.0,  z };

    r1*=(L/r1.norm());
    r2*=(L/r2.norm());
    r3*=(L/r3.norm());

    double R= r1.norm(); 
    std::vector<Eigen::Vector3d> r={r1,r2,r3}; 

    Lattice lattice = GenerateLattice(L, r);
    Tensor2D<Eigen::Vector3d> D = FindDerivative(r, lattice); 
    Tensor2D<int> TTList = TriangleVertexMapping(lattice.Vertices, L);

    for (int n=0; n< N_iter; n++)
    {
    double epsilon =1; 
    
    SpMat Op = AreaOperator(lattice, D, TTList, R); 
    Eigen::VectorXd Area = returnCurrentArea(lattice, TTList); 
    double rms = returnRMS(Area); 
    Eigen::VectorXd sol = setUpGMRESsolver(Op, Area, stepsize);
    printf("iter = %d ", n); 
    PrintGeometry(lattice, TTList, L);
    
    //printf("iter = %d\n", n);
    sol*=epsilon;  
    UpdateLattice(lattice, r, sol); 
    if (sol.norm()<=1e-12){break;}
    }
   
    FILE* testFull = NULL;  
    char test[64];
    sprintf(test, "%s/xivec_%d.dat", data_dir.c_str(), L); // filename
    testFull = fopen(test,"w");

    SaveBarycentric(lattice, testFull);
    fclose(testFull);
return 0;
}

