#include "CtIsing.h"
#include <numeric>
using namespace std;


//To Do: Checking that the statistics is respected when cutting ladder
//       Add wrappers to prevent contamination from nodes in different lists. 
int main(int argc, char *argv[])
{
    unsigned int seed = 1234u;
    QfeLattice lattice;
    int Nx = 2;
    int Ny = 2; 
    double k1 = 1.0; 
    double k2 = 1.0; 
    double invT = 5.0; 
    double t_coup = 2.0;
    lattice.SeedRng(seed);
    lattice.InitRect(Nx, Ny, k1, k2);
    ContinuousTimeLattice field(&lattice, invT, t_coup);
    CtIsing Ising(&field); 
    Ising.HotStart(); 
    for (int i=0; i<Ising.n_sites; i++)
    {
        std::cout<< i <<std::endl;
        field.TimeLadders[i].display();
    }
    printf("Cutting Ladder------\n");
    Ising.CutLadder();
    // field.TimeLadders[0].insertAtEnd(8.0);
    for (int i=0; i<Ising.n_sites; i++)
    {
        std::cout<< i <<std::endl;
        field.TimeLadders[i].display();
    }

    return 0; 
}
