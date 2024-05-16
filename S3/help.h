#include <iostream>
#include <math.h>
#include <vector>
#include <stack>
#include <set>



struct param{
  int L;
  int N;
  int EvenNum;
};

struct T{
  std::vector<int> Vert; 
  std::vector<int> ClassLabels; 
  int wt; 
};




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

int FindInd(int tetraVertex[][3], int n1, int n2, int n3, param p){
  int indx = 0; 
  for (int i =0; i<p.EvenNum; i++){
      if (n1==tetraVertex[i][0]&&n2==tetraVertex[i][1]&&n3==tetraVertex[i][2]){
        indx = i; 
        break;
      }
  }
  return indx;
}

//Check if we get correct total number of tetrahedrons!!
//std::vector<std::vector<int>> 
std::vector<T> FindTetrahedronRepresentative(std::vector<std::vector<int>>& orbit_reps, int tetraVertex[][3] ,param p){
  std::vector<std::vector<int>> tmp; 
  std::vector<T> result; 
  for(int lev =0; lev < p.L-1; lev++)
    {
	for(int n1 = 0; n1<= lev; n1++)
	  for(int n2 = 0; n2<= lev; n2++)
	    for(int n3 = 0; n3<= lev; n3++)
	      {	
		if (n1 + n2 + n3 == lev)
		  { 
        int v1 = FindRep(orbit_reps, n1, n2, n3, p); 
        int v2 = FindRep(orbit_reps, n1, n2, n3+1, p);
        int v3 = FindRep(orbit_reps, n1, n2+1, n3, p);
        int v4 = FindRep(orbit_reps, n1+1, n2, n3, p);
        std::vector<int> tetrahedron = {v1,v2,v3,v4}; 
        std::sort(tetrahedron.begin(), tetrahedron.end()); 
        auto it = std::find(tmp.begin(), tmp.end(), tetrahedron);
        if (it==tmp.end()){
          int ind1 = FindInd(tetraVertex, n1, n2, n3, p);    
          int ind2 = FindInd(tetraVertex, n1, n2, n3+1,p); 
          int ind3 = FindInd(tetraVertex, n1, n2+1, n3, p); 
          int ind4 = FindInd(tetraVertex, n1+1, n2, n3,p); 
          // std::vector<int> indList = {ind1, ind2, ind3, ind4}; 
          // std::sort(indList.begin()+1, indList.end());
          T tetra = {{ind1, ind2, ind3, ind4}, {v1,v2,v3, v4}, 1};
          tmp.push_back(tetrahedron);
          result.push_back(tetra); }
        else{
          int index = std::distance(tmp.begin(), it);
          result[index].wt+=1;
        }
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
