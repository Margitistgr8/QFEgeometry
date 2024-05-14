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
#include <stack>
using namespace std;

//#define L 4

struct param{
  int L;
  int N;
  int EvenNum;
};

void printEven(int tetraVertex[][3],  param p);

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
		    tetraVertex[n][0] =  n2 + n3; tetraVertex[n][1] =  n1 + n3; tetraVertex[n][2] =  n1 + n2;
		    n++; 
		  }
	      }
	printf("\n \n");
      }

  printEven(tetraVertex, p);
     
  return 0;
}

      
  
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
