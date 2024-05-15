/**********

r1 = (0,1,1), r2 = (1,0,1), r3 = (1,1,0)

Xindex(nx,ny,nz) = n1 * r1 + n2*r2 + n3 *r3
n1 + n2 + n3 <= L

integer n's gerated even lattice
n's all half integer are odd lattice


Nlayer[L_] := Sum[l + 1, {l, 0, L}] = 1/2 (1 + L) (2 + L)
NumVert[L_] := Sum[1/2  (1 + l)  (2 + l), {l, 0, L}]  = 1/6 (1 + L) (2 + L) (3 + L)

**********/

#include <iostream>
#include <math.h>
#include <vector>
#include <stack>
#include <set>
using namespace std;

//#define L 4

struct param{
  int L;
  int N;
  int EvenNum;
};

void printEven(int tetraVertex[][3],  param p);
std::vector<int> FindOrbitLabels(int tetraVertex[][3], param p);
std::vector<std::vector<int>>FindOrbitRepresentative(std::vector<int>& OrbitLabels, int tetraVertex[][3], param p);
int FindRep(std::vector<std::vector<int>>& orbit_reps, int n1, int n2, int n3, param p);
std::vector<std::vector<int>> FindRedTetrahedrons(std::vector<std::vector<int>>& orbit_reps, int tetraVertex[][3] ,param p);
//Need to fix FindRedTetrahedron Function , it gives label for tetrahedrons but also need a list with embedded vertices. 

int main( int argc, char *argv[])
{
  param p;
  p.L = 4;
  
  if(argc ==1)
    printf("In program  %s no input default %d \n",argv[0],p.L);
  else
    p.L = atoi(argv[1]);
  p.N = p.L*p.L*p.L;
   
  printf("L = %d  N = %d  \n", p.L,p.N);

  p.EvenNum = ((p.L)*(p.L+1)*(p.L +2))/6;

  cout << "  p.EvenNum "<<  p.EvenNum << endl;
  
  int tetraVertex[p.EvenNum][3];

  // int nx,ny, nz;

  int n = 0;
  
  for(int lev =0; lev < p.L; lev++)
      {
	for(int n1 = 0; n1<= lev; n1++)
	  for(int n2 = 0; n2<= lev; n2++)
	    for(int n3 = 0; n3<= lev; n3++)
	      {	
		if (n1 + n2 + n3 == lev)
		  { // nx = n2 + n3; ny = n1 + n3;  nz = n1 + n2;
		    printf(" (%d,%d,%d) ",  n2 + n3, n1 + n3, n1 + n2);
		    // This convention for printing out the lattice coordinates: 
        //tetraVertex[n][0] =  n2 + n3; tetraVertex[n][1] =  n1 + n3; tetraVertex[n][2] =  n1 + n2;
		    tetraVertex[n][0] =  n1; tetraVertex[n][1] =  n2; tetraVertex[n][2] =  n3;
        n++; 
		  }
	      }
	printf("\n \n");
      }

  printEven(tetraVertex, p);
  std::vector<int> orbitlabels = FindOrbitLabels(tetraVertex,p);
  for (int i = 0; i<p.EvenNum; i++)
  {std::cout<<orbitlabels[i]<<" "; }
  printf("\n"); 
  std::vector<std::vector<int>> orbit_reps = FindOrbitRepresentative(orbitlabels, tetraVertex, p); 
  for (size_t i = 0; i < orbit_reps.size(); ++i) {  // Loop over rows
        std::cout << i << " ";
    for (size_t j = 0; j < orbit_reps[i].size(); ++j) {  // Loop over columns
        std::cout << orbit_reps[i][j] << " ";
    }
    std::cout << std::endl;  // New line after each row
}
  std::cout << "Printing Tetrahedrons...\n"; 
  std::vector<std::vector<int>> tetra= FindRedTetrahedrons(orbit_reps, tetraVertex,p); 
  for (size_t i = 0; i < tetra.size(); ++i) {  // Loop over rows
        std::cout << i << " ";
    for (size_t j = 0; j < tetra[i].size(); ++j) {  // Loop over columns
        std::cout << tetra[i][j] << " ";
    }
    std::cout << std::endl;  // New line after each row
}  
  return 0;
}




std::vector<std::vector<int>>FindOrbitRepresentative(std::vector<int>& OrbitLabels, int tetraVertex[][3], param p){
  std::set<int> unique_labels(OrbitLabels.begin(), OrbitLabels.end());
  int n_orbit = unique_labels.size(); 
  std::vector<std::vector<int>> result; 
  for (int o=0; o<n_orbit; o++){
    auto it = std::find(OrbitLabels.begin(), OrbitLabels.end(), o);
    int index = std::distance(OrbitLabels.begin(), it);
    std::vector<int> pos = {tetraVertex[index][0],tetraVertex[index][1],tetraVertex[index][2],
                     p.L-1-tetraVertex[index][0]-tetraVertex[index][1]-tetraVertex[index][2]};
    result.push_back(pos); 
  }
  return result;
}



int FindRep(std::vector<std::vector<int>>& orbit_reps, int n1, int n2, int n3, param p){
  std::vector<int> v = {n1, n2, n3, p.L-1-n1-n2-n3}; 
  std::sort(v.begin(), v.end()); 
  auto it = std::find(orbit_reps.begin(), orbit_reps.end(), v);
  int index = std::distance(orbit_reps.begin(), it);
  return index;
}

//Check if we get correct total number of tetrahedrons!!
std::vector<std::vector<int>> FindRedTetrahedrons(std::vector<std::vector<int>>& orbit_reps, int tetraVertex[][3] ,param p){
  std::vector<std::vector<int>> result; 
  int vertexcounter = 0; 
  for(int lev =0; lev < p.L-1; lev++)
    {
	for(int n1 = 0; n1<= lev; n1++)
	  for(int n2 = 0; n2<= lev; n2++)
	    for(int n3 = 0; n3<= lev; n3++)
	      {	
		if (n1 + n2 + n3 == lev)
		  { 
        int v1 = FindRep(orbit_reps, n1, n2, n3, p); 
        int v2 = FindRep(orbit_reps, n1+1, n2, n3, p);
        int v3 = FindRep(orbit_reps, n1, n2+1, n3, p);
        int v4 = FindRep(orbit_reps, n1, n2, n3+1, p);
        std::vector<int> tetrahedron = {v1,v2,v3,v4}; 
        std::sort(tetrahedron.begin(), tetrahedron.end()); 
        auto it = std::find(result.begin(), result.end(), tetrahedron);
        if (it==result.end()){result.push_back(tetrahedron); }
		  }
	      }

      }
  return result; 
}

std::vector<int> FindOrbitLabels(int tetraVertex[][3], param p){

  std::vector<int> counter(p.EvenNum); //to keep track of points already labeled
  std::vector<int> classlabels(p.EvenNum); //a list that returns the label for all points 
  int dof_counter = 0; 
  for(int n = 0; n<p.EvenNum; n++){
      if (counter[n]!=0){continue;}
      else{
        classlabels[n]=dof_counter;
        counter[n] = 1; 
        int ncounted = 1; 
        std::vector<int> xi= {tetraVertex[n][0],tetraVertex[n][1],tetraVertex[n][2],
                    p.L-1-tetraVertex[n][0]-tetraVertex[n][1]-tetraVertex[n][2]};
        std::sort(xi.begin(), xi.end()); 
      // std::cout << xi[0] << ' ' << xi[1] << ' ' << xi[2]<< ' ' << xi[3] << '\n';
       // std::cout << "begin Check....\n";
        for (int n1=0; n1<p.EvenNum; n1++){
            if (ncounted==24){break;}
            if (counter[n1]!=0){continue;}
            else{
            std::vector<int>  xi2 = {tetraVertex[n1][0],tetraVertex[n1][1],tetraVertex[n1][2],
                    p.L-1-tetraVertex[n1][0]-tetraVertex[n1][1]-tetraVertex[n1][2]}; 
            std::sort(xi2.begin(),xi2.end());

            //std::cout << xi2[0] << ' ' << xi2[1] << ' ' << xi2[2]<< ' ' << xi2[3] << '\n';

            if (xi==xi2){
              classlabels[n1]=dof_counter;
              counter[n1] = 1; 
              ncounted += 1; 
              //std::cout << "Success....\n";
            }
            }
        }
        //std::cout << "end Check....\n";
        dof_counter+=1; 
      }
  }

  return classlabels;}
  
void printEven(int tetraVertex[][3], param p)
{
  FILE* fptr = NULL;  // C style
  char out_name[64];
  //  sprintf(out_name,"data/CGstate_%d_%d.dat",N,frame); // filename
  sprintf(out_name,"data/EvenVert_%d.dat",p.L); // filename
  fptr = fopen(out_name,"w");
  if(fptr == NULL)
    {
      printf("Error!");   
      exit(1);             
    }
  
  for(int n = 0; n<p.EvenNum; n++)
	{
	  fprintf(fptr, " %d ", n);
	  for(int k = 0; k<3; k++)
	    {
	      fprintf(fptr,"  %d ", tetraVertex[n][k]);
	    }
	  fprintf(fptr,"\n");
	}
  fclose(fptr);
}
