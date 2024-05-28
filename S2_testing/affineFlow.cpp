#include "triangle.h"


using namespace std;


#define L 3
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
    int N_iter=10; 
  
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
    
    std::cout<< lattice.rvec.size()<<' '<<lattice.xvec.size()<<' '<<lattice.xivec.size()<< "\n";
    std::cout<< lattice.VertexLabels.size()<<' '<<lattice.OrbitPos.size()<<' '<<lattice.Basis.size()<<"\n";
    std::cout<<lattice.ndof<<'\n'; 

    for (int i=0; i<lattice.Vertices.size(); i++)
    {
      printf("Vertex Position label %d with nx, ny, nz = %d %d %d\n", i, lattice.Vertices[i][0], lattice.Vertices[i][1], lattice.Vertices[i][2]);
    }

    for (int i=0; i<lattice.Vertices.size(); i++)
    {
      printf("Vertex Position label %d with nx, ny, nz = %d %d %d\n", i, lattice.Vertices[i][0], lattice.Vertices[i][1], lattice.Vertices[i][2]);
    }

    // for (Triangle TT: TList)
    // {
    //   std::cout<< TT.Vert[0]<<' '<< TT.Vert[1]<<' '<< TT.Vert[2] <<'\n'; 
    //    std::cout<< TT.ClassLabels[0]<<' '<< TT.ClassLabels[1]<<' '<< TT.ClassLabels[2] <<'\n'; 
    // }

    for (std::vector<int> TT: TTList)
    {
      std::cout<< TT[0]<<' '<< TT[1]<<' '<< TT[2] <<'\n'; 
     }


    Eigen::VectorXd totalArea = TotalAreaList(lattice, L);
    // std::cout<<totalArea<<'\n'; 
    for (int TT: lattice.VertexLabels)
    {
      std::cout<< TT <<'\n'; 
     }
    std::cout<< "Orbit Pos"<<'\n'; 
    for (int TT: lattice.OrbitPos)
    {
      std::cout<< TT <<'\n'; 
     }


    for (int n=0; n< N_iter; n++)
    {
    Eigen::VectorXd totalArea = TotalAreaList(lattice, L);
    double rms = returnRMS(totalArea); 
    printf("iter = %d, rms = %.12f\n", n, rms); 
    SpMat Op = AreaOperator(lattice, D, TList, R); 
    Eigen::VectorXd Area = returnCurrentArea(lattice, TList); 
    Eigen::VectorXd sol = setUpGMRESsolver(Op, Area, 1); 
    UpdateLattice(lattice, r, sol); 

    }
   
   totalArea = TotalAreaList(lattice, L);
    std::cout<<totalArea<<'\n'; 
    for (Eigen::Vector3d TT: lattice.xivec)
    {
      std::cout<< TT[0]<<' '<< TT[1]<<' '<< TT[2] <<'\n'; 
     }

     std::cout<<"Vectors\n";
    for (std::vector<int> TT: lattice.Vertices)
    {
      std::cout<< TT[0]<<' '<< TT[1]<<' '<< TT[2] <<'\n'; 
     }

return 0;
}

