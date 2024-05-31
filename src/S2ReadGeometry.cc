#include "S2.h"
using namespace std;

int main(int argc, char *argv[])
{
    int L = 4; 
    int q = 5; 
    std::string data_dir = "";
    if (argc==1)
    { std::cout << "Using default values: L = " << L << ", q = " << q  << std::endl;}
    else{
    if (argc<=3&&argc>1)
    {
    L= atoi(argv[1]);
    q= atoi(argv[2]);
    }
    else
    {
    L= atoi(argv[1]);
    q= atoi(argv[2]);
    data_dir = argv[3];   
    }}

    double z = z = sqrt((7 + 3 * sqrt(5))/8);
    if (q==4){z = 1/sqrt(2); }
    if (q==3){z = 1/(2*sqrt(2)); }
    Eigen::Vector3d r1 =  { sqrt(3)/2.0 , -1/2.0,  z };
    Eigen::Vector3d r2 =  {           0,   1.0 ,   z };
    Eigen::Vector3d r3 =  { -sqrt(3)/2.0, -1/2.0,  z };

    r1*=(L/r1.norm());
    r2*=(L/r2.norm());
    r3*=(L/r3.norm());

    double R= r1.norm(); 
    std::vector<Eigen::Vector3d> r={r1,r2,r3}; 

    Lattice lattice = GenerateLattice(L, r);
    Tensor2D<int> TTList = TriangleVertexMapping(lattice.Vertices, L);
    
    if (data_dir.empty())
    {
    //printf("Initialize lattice-----\n"); 
    printf("%d ", L); 
    PrintGeometry(lattice, TTList, L);
    }

    else
    {
    //printf("Downloading Orbit-----\n"); 
    ReadBarycentric(lattice, data_dir, double(L), r); 
    printf("%d ", L); 
    PrintGeometry(lattice, TTList, L);
    }

    // FILE* testFull = NULL;  
    // char test[128];
    // sprintf(test, "/Users/jinyunlin/QFE-Research/QFEgeometry/S2_new/q4/test_xivec_%d.dat", L); // filename
    // testFull = fopen(test,"w");
    // SaveBarycentric(lattice, testFull);
    // fclose(testFull);
return 0;
}

