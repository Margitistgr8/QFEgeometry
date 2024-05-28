template <int L>
void SubtractAverage(std::vector<T>& tripletList, int dof){
  double av[dof] ; 
for (int i=0; i<dof; i++){
  av[i]=0; 
}
for (size_t i = 0; i < tripletList.size(); ++i) {
    const auto& triplet = tripletList[i];
    av[triplet.row()]+=triplet.value(); 
}
for (int i=0; i<dof; i++){
  av[i]/=L*L; 
}
for (size_t i = 0; i < tripletList.size(); ++i) {
    const auto& triplet = tripletList[i];
    tripletList[i] = Eigen::Triplet<double>(triplet.row(), triplet.col(), triplet.value()-av[triplet.row()]);
}
}

template <int L>
void CurrentDualArea(int (&arr)[L*L][3], double (&rvec)[L+1][L+1][3], double (&dualArea)[L+1][L+1]){
for (int Tri=0; Tri<L*L; Tri++){
    double r[3][3] =  {};
      for (int vert = 0; vert<3; vert++){
          int s = arr[Tri][vert];
          int nx = s % (L+1); 
          int ny = s / (L+1);
          for (int mu=0; mu <3; mu++){
          r[vert][mu] = rvec[ny][nx][mu];
          //printf("%.12f\n", r[vert][mu]);
      }
      }
          double mside1[3] =  {0,0,0};
          double mside2[3] =  {0,0,0};
          double mside3[3] =  {0,0,0};

          for (int mu=0; mu<3; mu++){
            mside1[mu] = r[0][mu]-r[1][mu];
            mside2[mu] = r[0][mu]-r[2][mu];
            mside3[mu] = mside1[mu]-mside2[mu];
          } 

          double a = dotproduct3D(mside1, mside1);
          double b = dotproduct3D(mside2, mside2);
          double c = dotproduct3D(mside3, mside3);
          double astar = Compute_ls(a,b,c,a); 
          double bstar = Compute_ls(a,b,c,b); 
          double cstar = Compute_ls(a,b,c,c); 
          
          int s = arr[Tri][0];
          int nx = s % (L+1); 
          int ny = s / (L+1);
          dualArea[ny][nx] += (sqrt(a)*astar+sqrt(b)*bstar)/2.0; 

          int s1 = arr[Tri][1];
          int nx1 = s1 % (L+1); 
          int ny1 = s1 / (L+1);
          dualArea[ny1][nx1] += (sqrt(a)*astar+sqrt(c)*cstar)/2.0; 

          int s2 = arr[Tri][2];
          int nx2 = s2 % (L+1); 
          int ny2 = s2 / (L+1);
          dualArea[ny2][nx2] += (sqrt(b)*bstar+sqrt(c)*cstar)/2.0; 
    }
  for (int nx=0; nx<=L; nx++){
      for (int ny=0; ny<=L-nx; ny++){
        int nz = L-ny-nx; 
        if (nx==L||ny==L||nz==L){
          dualArea[ny][nx]*=5.0; 
        }
        else if (nx==0||ny==0||nz==0){
          dualArea[ny][nx]*=2.0; 
        }
      }
  }
#if 0
    for (int nx=0; nx<=L; nx++){
      for (int ny=0; ny<=L-nx; ny++){
        printf("dual area at nx,ny (%d %d) is %.12f\n", nx, ny, dualArea[ny][nx]); 
      }}  
#endif
  
}

template <int L>
void CurrentDeficitAngle(int (&arr)[L*L][3], double (&rvec)[L+1][L+1][3], double (&defangle)[L+1][L+1]){
for (int Tri=0; Tri<L*L; Tri++){
    double r[3][3] =  {};
      for (int vert = 0; vert<3; vert++){
          int s = arr[Tri][vert];
          int nx = s % (L+1); 
          int ny = s / (L+1);
          for (int mu=0; mu <3; mu++){
          r[vert][mu] = rvec[ny][nx][mu];
          //printf("%.12f\n", r[vert][mu]);
      }
      }
          double mside1[3] =  {0,0,0};
          double mside2[3] =  {0,0,0};
          double mside3[3] =  {0,0,0};

          for (int mu=0; mu<3; mu++){
            mside1[mu] = r[0][mu]-r[1][mu];
            mside2[mu] = r[0][mu]-r[2][mu];
            mside3[mu] = mside1[mu]-mside2[mu];
          } 

          double a = sqrt(dotproduct3D(mside1, mside1));
          double b = sqrt(dotproduct3D(mside2, mside2));
          double c = sqrt(dotproduct3D(mside3, mside3));
          double perimeter = a+b+c; 
          double ratio; 

          int s = arr[Tri][0];
          int nx = s % (L+1); 
          int ny = s / (L+1); 
          ratio = (perimeter-2*b)*(perimeter-2*a)/(perimeter*(perimeter-2*c)); 
          defangle[ny][nx] -= 2.0*atan(sqrt(ratio)); 

          int s1 = arr[Tri][1];
          int nx1 = s1 % (L+1); 
          int ny1 = s1 / (L+1);
          ratio = (perimeter-2*c)*(perimeter-2*a)/(perimeter*(perimeter-2*b)); 
          defangle[ny1][nx1] -= 2.0*atan(sqrt(ratio)); 

          int s2 = arr[Tri][2];
          int nx2 = s2 % (L+1); 
          int ny2 = s2 / (L+1);
          ratio = (perimeter-2*c)*(perimeter-2*b)/(perimeter*(perimeter-2*a)); 
          defangle[ny2][nx2] -= 2.0*atan(sqrt(ratio)); 
    }
  for (int nx=0; nx<=L; nx++){
      for (int ny=0; ny<=L-nx; ny++){
        int nz = L-ny-nx; 
        if (nx==L||ny==L||nz==L){
          defangle[ny][nx]*=5.0; 
        }
        else if (nx==0||ny==0||nz==0){
          defangle[ny][nx]*=2.0; 
        }
      }
  }  
    for (int nx=0; nx<=L; nx++){
      for (int ny=0; ny<=L-nx; ny++){
       defangle[ny][nx]+=TWOPI; 
      }
  }  
}


template <int L>
void SaveFaceFunction(Eigen::VectorXd& quantity, FILE* file, int (&map)[L*L][3]){
  //default is to save the face quantities on the center of each face
  double r1[2] = {1.0,0}; 
  double r2[2] = {1/2.0, sqrt(3)/2}; 
  double r3[2] = {0,0}; 
  for (int Tri=0; Tri<L*L; Tri++){
    int xi[3] = {map[Tri][0], map[Tri][1], map[Tri][2]}; 
    for (int mu = 0; mu<3; mu++){
    int x = xi[mu]% (L+1); 
    int y = xi[mu]/ (L+1);
    fprintf(file, "%.4f %.4f ", r1[0]*x+r2[0]*y,r1[1]*x+r2[1]*y);
    }
    fprintf(file, "%.16f\n",quantity(Tri)); 
  }
}

template <int L>
void SaveCurvedFaceFunction(Eigen::VectorXd& quantity, FILE* file, int (&map)[L*L][3], double (&rvec)[L+1][L+1][3]){
  //default is to save the face quantities on the center of each face
  double r1[2] = {1.0,0}; 
  double r2[2] = {1/2.0, sqrt(3)/2}; 
  double r3[2] = {0,0}; 
  for (int Tri=0; Tri<L*L; Tri++){
    int xi[3] = {map[Tri][0], map[Tri][1], map[Tri][2]}; 
    for (int mu = 0; mu<3; mu++){
    int x = xi[mu]% (L+1); 
    int y = xi[mu]/ (L+1);
    fprintf(file, "%.12f %.12f %.12f ", rvec[y][x][0],rvec[y][x][1],rvec[y][x][2]);
    }
    fprintf(file, "%.16f\n",quantity(Tri)); 
  }
}
template<int L>
void PrintDualInfo(int (&arr)[L*L][3], double (&rvec)[L+1][L+1][3]){
  double dualArea[L+1][L+1] = {};
  CurrentDualArea<L>(arr, rvec, dualArea); 
  double mean = 0.0; 
  double mean_sq = 0.0; 
  for (int nx=0; nx<=L; nx++){
      for (int ny=0; ny<=L-nx; ny++){
        int nz = L-ny-nx; 
        if (nx==L||ny==L||nz==L){
          mean += 4.0*dualArea[ny][nx]; 
          mean_sq+=4.0*dualArea[ny][nx]*dualArea[ny][nx]; 
        }
        else if (nx==0||ny==0||nz==0){
          mean +=  10.0*dualArea[ny][nx];
          mean_sq += 10.0*dualArea[ny][nx]*dualArea[ny][nx];
        }
        else{
          mean += 20.0*dualArea[ny][nx];
          mean_sq += 20.0*dualArea[ny][nx]*dualArea[ny][nx];
        }     }
  } 
  mean/=double(2+10*L*L); 
  mean_sq/=double(2+10*L*L); 
  double rms = mean_sq-mean*mean; 
  printf(" %.16f %.16f %.16f ", mean, mean_sq, rms); 
}

template<int L>
void PrintDeficitInfo(int (&arr)[L*L][3], double (&rvec)[L+1][L+1][3]){
  double defangle[L+1][L+1] = {};
  CurrentDeficitAngle<L>(arr, rvec, defangle); 
  double mean = 0.0; 
  double mean_sq = 0.0; 
  for (int nx=0; nx<=L; nx++){
      for (int ny=0; ny<=L-nx; ny++){
        int nz = L-ny-nx; 
        if (nx==L||ny==L||nz==L){
          mean += 4.0*defangle[ny][nx]; 
          mean_sq+=4.0*defangle[ny][nx]*defangle[ny][nx]; 
        }
        else if (nx==0||ny==0||nz==0){
          mean +=  10.0*defangle[ny][nx];
          mean_sq += 10.0*defangle[ny][nx]*defangle[ny][nx];
        }
        else{
          mean += 20.0*defangle[ny][nx];
          mean_sq += 20.0*defangle[ny][nx]*defangle[ny][nx];
        }
      }
  } 
  mean/=double(2+10*L*L); 
  mean_sq/=double(2+10*L*L); 
  double rms = mean_sq-mean*mean; 
  printf(" %.16f %.16f %.16f ", mean, mean_sq, rms); 
}

template <int L>
int BuildBasisVectorInfo(int (&arr)[L+1][L+1][3]){
     int Poscounter = 0; 
    for (int ny =0; ny<=L; ny++){
        for (int nx =0; nx<=L-ny; nx++){
            int nz = L-ny-nx; 
            if (nz==L||nx==L||ny==L){continue;} //No update if the vertex is one of the icosahedral vertices
            //No update if the vertex is one of the midpoints
            if (L%2==0&&nx== L/2 && ny==L/2){continue;} 
            if (L%2==0&& nx==0 && ny==L/2){continue;} 
            if (L%2==0&&nx==L/2 && ny==0){continue;} 
            if (nx==ny&& ny==nz){continue;}//Skip midpoint
            else{
                arr[nx][ny][0]=Poscounter; 
                if (nx==0){
                    arr[nx][ny][1]=1; 
                    arr[nx][ny][2]=1; //Derivative w.r.t xi 2
                    Poscounter+=1; 
                }
                else if (ny==0){
                    arr[nx][ny][1]=1; 
                    arr[nx][ny][2]=0;  //Derivative w.r.t xi 1
                    Poscounter+=1; 
                }
                else if (nz==0){            
                    arr[nx][ny][1]=1; 
                    arr[nx][ny][2]=0; 
                    Poscounter+=1; }
                else {               
                    if (nx==ny){//Derivative w.r.t. xi1
                    arr[nx][ny][1]=1; 
                    arr[nx][ny][2]=0; 
                    Poscounter+=1;
                    } 
                    else if (nx==nz){//Derivative w.r.t. xi1 
                    arr[nx][ny][1]=1; 
                    arr[nx][ny][2]=0; 
                    Poscounter+=1;
                    } 
                    else if (ny==nz){//Derivative w.r.t. xi2 
                    arr[nx][ny][1]=1; 
                    arr[nx][ny][2]=1; 
                    Poscounter+=1;
                    }  
                    else{
                    arr[nx][ny][1]=2; 
                    arr[nx][ny][2]=2; 
                    Poscounter+=2;}
                    }

            }
        }
    }
    return Poscounter;
}

template <int L>
Eigen::VectorXd CurrentAreaList(int (&arr)[L*L][3], double (&rvec)[L+1][L+1][3]){
    Eigen::VectorXd vec(L*L); 
    for (int Tri=0; Tri<L*L; Tri++){
    double r[3][3] =  {};

      for (int vert = 0; vert<3; vert++){
          int s = arr[Tri][vert];
          int nx = s % (L+1); 
          int ny = s / (L+1);
          for (int mu=0; mu <3; mu++){
          r[vert][mu] = rvec[ny][nx][mu];
          //printf("%.12f\n", r[vert][mu]);
      }
      }
          double mside1[3] =  {0,0,0};
          double mside2[3] =  {0,0,0};
          double mside3[3] =  {0,0,0};

          for (int mu=0; mu<3; mu++){
            mside1[mu] = r[0][mu]-r[1][mu];
            mside2[mu] = r[0][mu]-r[2][mu];
            mside3[mu] = mside1[mu]-mside2[mu];
          } 

          double a = dotproduct3D(mside1, mside1);
          double b = dotproduct3D(mside2, mside2);
          double c = dotproduct3D(mside3, mside3);
          //printf("Sides are given %.4f %.4f %.4f\n", a,b,c);
     vec(Tri) = A(a,b,c);   
    }
    return vec; 
}

template <int L>
Eigen::VectorXd CurrentPerimeterList(int (&arr)[L*L][3], double (&rvec)[L+1][L+1][3]){
  Eigen::VectorXd vec(L*L); 
  for (int Tri=0; Tri<L*L; Tri++){
    double r[3][3] =  {};

      for (int vert = 0; vert<3; vert++){
          int s = arr[Tri][vert];
          int nx = s % (L+1); 
          int ny = s / (L+1);
          for (int mu=0; mu <3; mu++){
          r[vert][mu] = rvec[ny][nx][mu];
      }
      }
          double mside1[3] =  {0,0,0};
          double mside2[3] =  {0,0,0};
          double mside3[3] =  {0,0,0};

          for (int mu=0; mu<3; mu++){
            mside1[mu] = r[0][mu]-r[1][mu];
            mside2[mu] = r[0][mu]-r[2][mu];
            mside3[mu] = mside1[mu]-mside2[mu];
          } 

          double a = dotproduct3D(mside1, mside1);
          double b = dotproduct3D(mside2, mside2);
          double c = dotproduct3D(mside3, mside3);
          //printf("Sides are given %.4f %.4f %.4f\n", a,b,c);
          vec(Tri) = (sqrt(a)+sqrt(b)+sqrt(c));   
    }
    return vec; 
}

template <int L>
void SetPosition(double (&xivec)[L+1][L+1][3],double (&rvec)[L+1][L+1][3], double r1[3], double r2[3],double r3[3],double (&xvec)[L+1][L+1][3]){
double Rsqr = r1[0]*r1[0]+r1[1]*r1[1]+r1[2]*r1[2]; 
//printf("Radius of point in update given: %.12f\n", sqrt(Rsqr)); 
for(int ny = 0; ny<= L; ny++)
    {
      for(int nx = 0; nx <= L; nx++)
	{   
    double x_sqr = 0.0; 
        for (int mu = 0; mu<3; mu++){
          xvec[ny][nx][mu] =( xivec[ny][nx][0]*r1[mu] +   xivec[ny][nx][1]*r2[mu] +  xivec[ny][nx][2]*r3[mu] ); 
          x_sqr+= xvec[ny][nx][mu]*xvec[ny][nx][mu];
        }
	      
	  for(int mu =0; mu < 3; mu++)
	    {	    
        rvec[ny][nx][mu]=sqrt(Rsqr)*xvec[ny][nx][mu]/sqrt(x_sqr);
	    }
	}
  } 

}

template <int L>
void TriangleVertexMapping(int (&arr)[L*L][3]){
    int rowcounter = 0; 
    int arrPosCounter = 0; 
    for (int layer = L; layer>0; layer--){
        
        int numOftriangle = 1+(layer-1)*2;
        // printf("Layer: %d, Num of Triangle: %d\n", layer, numOftriangle);
        int nx = 0; 
        int ny = rowcounter; 
        int poscounter = 0; 
        // int oddcounter = 0; 
        // int evencounter = 0; 
    for (int point = 0; point<numOftriangle; point++){
        if (point%2==0){
            // int x1 = nx+evencounter; 
            // int y1 = ny+evencounter; 
            arr[arrPosCounter][0] = returnSerial(nx, ny,L);
            arr[arrPosCounter][1] = returnSerial(nx+1, ny, L);
            arr[arrPosCounter][2] = returnSerial(nx, ny+1, L);
            ny+=1; 
        }
        else{
            arr[arrPosCounter][0] = returnSerial(nx, ny, L);
            arr[arrPosCounter][1] = returnSerial(nx+1, ny, L);
            arr[arrPosCounter][2] = returnSerial(nx+1, ny-1, L);
            ny-=1; 
            nx+=1;
        }
        arrPosCounter +=1; 
    }
        rowcounter+=1; 
    }
}


template <int L>
void BuildAreaOperator(std::vector<T>& (tripletList),int (&map)[L*L][3],int (&vecInf)[L+1][L+1][3],double(&rvec)[L+1][L+1][3],double(&xvec)[L+1][L+1][3], double r12[3], double r13[3], double r23[3])
{
  double b[3] = {rvec[0][0][0],rvec[0][0][1],rvec[0][0][2]};
  double R = sqrt(dotproduct3D(b,b)); 
  
  for (int Tri = 0; Tri<L*L; Tri++){//Curved Space Formula
      double r[3][3] =  {};
      for (int vert = 0; vert<3; vert++){
          int s = map[Tri][vert];
          int nx = s % (L+1); 
          int ny = s / (L+1);
          for (int mu=0; mu <3; mu++){
          r[vert][mu] = rvec[ny][nx][mu];}
      }
      for (int vert = 0; vert<3; vert++){
  
          int s = map[Tri][vert];
          int nx = s % (L+1); 
          int ny = s / (L+1);
          int nz = L-nx-ny;
          int pos = vecInf[nx][ny][0];
          //printf("Looking at vertex: %d %d at %.8f %.8f %.8f\n", nx, ny, rvec[ny][nx][0],rvec[ny][nx][1],rvec[ny][nx][2]);
          if (pos==-1){continue;}    //Skip if vertex is not a d.o.f

          double side1[3] =  {0,0,0};
          double side2[3] =  {0,0,0};
          double side3[3] =  {0,0,0};

          int ind1 = (vert+1)%3;
          int ind2 = (vert+2)%3;

          for (int mu=0; mu<3; mu++){
            side1[mu] = r[vert][mu]-r[ind1][mu];
            side2[mu] = r[vert][mu]-r[ind2][mu];
            side3[mu] = side1[mu]-side2[mu];
          } 

          double l1 = dotproduct3D(side1, side1);
          double l2 = dotproduct3D(side2, side2);
          double l3 = dotproduct3D(side3, side3);
          // printf("Subtraingle given: %.4f %.4f %.4f\n", l1, l2, l3);
          double l1star = Compute_ls(l1, l2, l3, l1);
          double l2star = Compute_ls(l1, l2, l3, l2);

          double dRnx[3]={0,0,0};
          double dRny[3]={0,0,0};
          double dRnz[3]={0,0,0};

          double xpos[3]={xvec[ny][nx][0], xvec[ny][nx][1], xvec[ny][nx][2]};
          const double xnorm = sqrt(dotproduct3D(xpos, xpos));

          for (int mu=0; mu<3; mu++){
            dRnx[mu]+=(r13[mu])/xnorm - (xpos[mu]* dotproduct3D(r13, xpos))/(pow(xnorm, 3.0));
            dRny[mu]+=(r23[mu])/xnorm - (xpos[mu]* dotproduct3D(r23, xpos))/(pow(xnorm, 3.0));
            dRnz[mu]+=(r12[mu])/xnorm - (xpos[mu]* dotproduct3D(r12, xpos))/(pow(xnorm, 3.0));
            dRnx[mu]*=R; 
            dRny[mu]*=R; 
            dRnz[mu]*=R;
          }
    
          if (vecInf[nx][ny][1]==1){
            if (vecInf[nx][ny][2]==0){//Take derivative w.r.t xi 1           
            if (ny==0){
              double val = 2*dotproduct3D(dRnx, side1)*l1star/sqrt(l1)+2*dotproduct3D(dRnx, side2)*l2star/sqrt(l2);
              tripletList.push_back(Eigen::Triplet<double>(pos, Tri, val));
              }
            else if (nz==0){
              double val = 2*dotproduct3D(dRnz, side1)*l1star/sqrt(l1)+2*dotproduct3D(dRnz, side2)*l2star/sqrt(l2);
              tripletList.push_back(Eigen::Triplet<double>(pos, Tri, val)); 
              }
            else if (nx==ny){
            double val = 2*dotproduct3D(dRnx, side1)*l1star/sqrt(l1)+2*dotproduct3D(dRnx, side2)*l2star/sqrt(l2);
            val += 2*dotproduct3D(dRny, side1)*l1star/sqrt(l1)+2*dotproduct3D(dRny, side2)*l2star/sqrt(l2);
            tripletList.push_back(Eigen::Triplet<double>(pos, Tri, val));              
            }
            else{//nx==nz
            double val = 2*dotproduct3D(dRnz, side1)*l1star/sqrt(l1)+2*dotproduct3D(dRnz, side2)*l2star/sqrt(l2);
            val -= 2*dotproduct3D(dRny, side1)*l1star/sqrt(l1)+2*dotproduct3D(dRny, side2)*l2star/sqrt(l2);
            tripletList.push_back(Eigen::Triplet<double>(pos, Tri, val));  
            }
            }
            else{//Take derivative w.r.t xi 2
            if (nx==0){
              double val = 2*dotproduct3D(dRny, side1)*l1star/sqrt(l1)+2*dotproduct3D(dRny, side2)*l2star/sqrt(l2);
              tripletList.push_back(Eigen::Triplet<double>(pos, Tri, val));}
              else{//ny==nz
              double val = -2*dotproduct3D(dRnz, side1)*l1star/sqrt(l1)-2*dotproduct3D(dRnz, side2)*l2star/sqrt(l2);
              val -= 2*dotproduct3D(dRnx, side1)*l1star/sqrt(l1)+2*dotproduct3D(dRnx, side2)*l2star/sqrt(l2);
              tripletList.push_back(Eigen::Triplet<double>(pos, Tri, val));
              }         
            }
          }
          else{
          double val1 = 2*dotproduct3D(dRnx, side1)*l1star/sqrt(l1)+2*dotproduct3D(dRnx, side2)*l2star/sqrt(l2);
          tripletList.push_back(Eigen::Triplet<double>(pos, Tri, val1)); 
          //Take derivative w.r.t xi 1 and xi 2
          double val2 = 2*dotproduct3D(dRny, side1)*l1star/sqrt(l1)+2*dotproduct3D(dRny, side2)*l2star/sqrt(l2);
          tripletList.push_back(Eigen::Triplet<double>(pos+1, Tri, val2));              
          }
      }
    }
}

template <int L>
void BuildPerimeterOperator(std::vector<T>& (tripletList),int (&map)[L*L][3],int (&vecInf)[L+1][L+1][3],double(&rvec)[L+1][L+1][3],double(&xvec)[L+1][L+1][3],double r12[3], double r13[3], double r23[3])
{
  double b[3] = {rvec[0][0][0],rvec[0][0][1],rvec[0][0][2]};
  double R = sqrt(dotproduct3D(b,b)); 
 
  for (int Tri = 0; Tri<L*L; Tri++){//Curved Space Formula
      double r[3][3] =  {};
      for (int vert = 0; vert<3; vert++){
          int s = map[Tri][vert];
          int nx = s % (L+1); 
          int ny = s / (L+1);
          for (int mu=0; mu <3; mu++){
          r[vert][mu] = rvec[ny][nx][mu];}
      }
      for (int vert = 0; vert<3; vert++){
  
          int s = map[Tri][vert];
          int nx = s % (L+1); 
          int ny = s / (L+1);
          int nz = L-nx-ny;
          int pos = vecInf[nx][ny][0];
          //printf("Looking at vertex: %d %d at %.8f %.8f %.8f\n", nx, ny, rvec[ny][nx][0],rvec[ny][nx][1],rvec[ny][nx][2]);
          if (pos==-1){continue;}    //Skip if vertex is not a d.o.f
          double side1[3] =  {0,0,0};
          double side2[3] =  {0,0,0};
          double side3[3] =  {0,0,0};

          int ind1 = (vert+1)%3;
          int ind2 = (vert+2)%3;

          for (int mu=0; mu<3; mu++){
            side1[mu] = r[vert][mu]-r[ind1][mu];
            side2[mu] = r[vert][mu]-r[ind2][mu];
            side3[mu] = side1[mu]-side2[mu];
          } 

          double l1 = dotproduct3D(side1, side1);
          double l2 = dotproduct3D(side2, side2);
          double l3 = dotproduct3D(side3, side3);
          double perimeter = sqrt(l1)+sqrt(l2)+sqrt(l3);

          double dRnx[3]={0,0,0};
          double dRny[3]={0,0,0};
          double dRnz[3]={0,0,0};

          double xpos[3]={xvec[ny][nx][0], xvec[ny][nx][1], xvec[ny][nx][2]};
          const double xnorm = sqrt(dotproduct3D(xpos, xpos));

          for (int mu=0; mu<3; mu++){
            dRnx[mu]+=(r13[mu])/xnorm - (xpos[mu]* dotproduct3D(r13, xpos))/(pow(xnorm, 3.0));
            dRny[mu]+=(r23[mu])/xnorm - (xpos[mu]* dotproduct3D(r23, xpos))/(pow(xnorm, 3.0));
            dRnz[mu]+=(r12[mu])/xnorm - (xpos[mu]* dotproduct3D(r12, xpos))/(pow(xnorm, 3.0));
            dRnx[mu]*=R; 
            dRny[mu]*=R; 
            dRnz[mu]*=R;          
          }
    
          if (vecInf[nx][ny][1]==1){
            if (vecInf[nx][ny][2]==0){//Take derivative w.r.t xi 1 
            if (ny==0){
              double val = dotproduct3D(dRnx, side1)/sqrt(l1)+dotproduct3D(dRnx, side2)/sqrt(l2);
       
              tripletList.push_back(Eigen::Triplet<double>(pos, Tri, val));
              }
            
            else if (nz==0){//nz==0
              double val = dotproduct3D(dRnz, side1)/sqrt(l1)+dotproduct3D(dRnz, side2)/sqrt(l2);
        
              tripletList.push_back(Eigen::Triplet<double>(pos, Tri, val)); 
              }
            else if (nx==ny){
              double val = dotproduct3D(dRnx, side1)/sqrt(l1)+dotproduct3D(dRnx, side2)/sqrt(l2);
              val += dotproduct3D(dRny, side1)/sqrt(l1)+dotproduct3D(dRny, side2)/sqrt(l2);
        
              tripletList.push_back(Eigen::Triplet<double>(pos, Tri, val)); 

            }
            else{
              double val = dotproduct3D(dRnz, side1)/sqrt(l1)+dotproduct3D(dRnz, side2)/sqrt(l2);
              val -= dotproduct3D(dRny, side1)/sqrt(l1)+dotproduct3D(dRny, side2)/sqrt(l2);
    
              tripletList.push_back(Eigen::Triplet<double>(pos, Tri, val)); 

            }
            }
            else{//Take derivative w.r.t xi 2
            if (nx==0){
              double val = dotproduct3D(dRny, side1)/sqrt(l1)+dotproduct3D(dRny, side2)/sqrt(l2);
           
              tripletList.push_back(Eigen::Triplet<double>(pos, Tri, val));}
            else{
              double val = -dotproduct3D(dRnz, side1)/sqrt(l1)-dotproduct3D(dRnz, side2)/sqrt(l2);
              val -= dotproduct3D(dRnx, side1)/sqrt(l1)+dotproduct3D(dRnx, side2)/sqrt(l2);
             
              tripletList.push_back(Eigen::Triplet<double>(pos, Tri, val)); 
            }          
            }
          }
          else{
          double val1 = dotproduct3D(dRnx, side1)/sqrt(l1)+dotproduct3D(dRnx, side2)/sqrt(l2);
     
          tripletList.push_back(Eigen::Triplet<double>(pos, Tri, val1)); 
          //Take derivative w.r.t xi 1 and xi 2
          double val2 = dotproduct3D(dRny, side1)/sqrt(l1)+dotproduct3D(dRny, side2)/sqrt(l2);
         
          tripletList.push_back(Eigen::Triplet<double>(pos+1, Tri, val2));              
          }
      }
    }
}


template <int L>
Eigen::VectorXd CurrentCircumradiusList(int (&arr)[L*L][3], double (&rvec)[L+1][L+1][3]){
    Eigen::VectorXd vec(L*L); 
    for (int Tri=0; Tri<L*L; Tri++){
    double r[3][3] =  {};

      for (int vert = 0; vert<3; vert++){
          int s = arr[Tri][vert];
          int nx = s % (L+1); 
          int ny = s / (L+1);
          for (int mu=0; mu <3; mu++){
          r[vert][mu] = rvec[ny][nx][mu];
          //printf("%.12f\n", r[vert][mu]);
      }
      }
          double mside1[3] =  {0,0,0};
          double mside2[3] =  {0,0,0};
          double mside3[3] =  {0,0,0};

          for (int mu=0; mu<3; mu++){
            mside1[mu] = r[0][mu]-r[1][mu];
            mside2[mu] = r[0][mu]-r[2][mu];
            mside3[mu] = mside1[mu]-mside2[mu];
          } 

          double a = dotproduct3D(mside1, mside1);
          double b = dotproduct3D(mside2, mside2);
          double c = dotproduct3D(mside3, mside3);
          //printf("Sides are given %.4f %.4f %.4f\n", a,b,c);
     vec(Tri) = Cr(a,b,c);
    }
    return vec; 
}



template <int L>
void BuildCircumradiusOperator(std::vector<T>& (tripletList),int (&map)[L*L][3],int (&vecInf)[L+1][L+1][3],double(&rvec)[L+1][L+1][3], double(&xvec)[L+1][L+1][3],double r12[3], double r13[3], double r23[3])
{
  double b[3] = {rvec[0][0][0],rvec[0][0][1],rvec[0][0][2]};
  double R = sqrt(dotproduct3D(b,b)); 
  
  for (int Tri = 0; Tri<L*L; Tri++){//Curved Space Formula
      double r[3][3] =  {};
      for (int vert = 0; vert<3; vert++){
          int s = map[Tri][vert];
          int nx = s % (L+1); 
          int ny = s / (L+1);
          for (int mu=0; mu <3; mu++){
          r[vert][mu] = rvec[ny][nx][mu];}
      }
      for (int vert = 0; vert<3; vert++){
  
          int s = map[Tri][vert];
          int nx = s % (L+1); 
          int ny = s / (L+1);
          int nz = L-nx-ny;
          int pos = vecInf[nx][ny][0];
          //printf("Looking at vertex: %d %d at %.8f %.8f %.8f\n", nx, ny, rvec[ny][nx][0],rvec[ny][nx][1],rvec[ny][nx][2]);
          if (pos==-1){continue;}    //Skip if vertex is not a d.o.f
          


          double side1[3] =  {0,0,0};
          double side2[3] =  {0,0,0};
          double side3[3] =  {0,0,0};

          int ind1 = (vert+1)%3;
          int ind2 = (vert+2)%3;

          for (int mu=0; mu<3; mu++){
            side1[mu] = r[vert][mu]-r[ind1][mu];
            side2[mu] = r[vert][mu]-r[ind2][mu];
            side3[mu] = side1[mu]-side2[mu];
          } 

          double l1 = dotproduct3D(side1, side1);
          double l2 = dotproduct3D(side2, side2);
          double l3 = dotproduct3D(side3, side3);
          double Area = A(l1, l2, l3);
          double l1star = Compute_ls(l1, l2, l3, l1);
          double l2star = Compute_ls(l1, l2, l3, l2);
          double cr = Cr(l1,l2,l3);

          double edge1factor = cr*(1/l1 - (2*l1star)/(sqrt(l1)*Area));
          double edge2factor = cr*(1/l2 - (2*l2star)/(sqrt(l2)*Area));


          double dRnx[3]={0,0,0};
          double dRny[3]={0,0,0};
          double dRnz[3]={0,0,0};

          double xpos[3]={xvec[ny][nx][0], xvec[ny][nx][1], xvec[ny][nx][2]};
          const double xnorm = sqrt(dotproduct3D(xpos, xpos));

          for (int mu=0; mu<3; mu++){
            dRnx[mu]+=(r13[mu])/xnorm - (xpos[mu]* dotproduct3D(r13, xpos))/(pow(xnorm, 3.0));
            dRny[mu]+=(r23[mu])/xnorm - (xpos[mu]* dotproduct3D(r23, xpos))/(pow(xnorm, 3.0));
            dRnz[mu]+=(r12[mu])/xnorm - (xpos[mu]* dotproduct3D(r12, xpos))/(pow(xnorm, 3.0));
            dRnx[mu]*=R; 
            dRny[mu]*=R; 
            dRnz[mu]*=R; 
          }
   
          if (vecInf[nx][ny][1]==1){
            if (vecInf[nx][ny][2]==0){//Take derivative w.r.t xi 1 
            if (ny==0){
              double val = dotproduct3D(dRnx, side1)*edge1factor+dotproduct3D(dRnx, side2)*edge2factor;
              tripletList.push_back(Eigen::Triplet<double>(pos, Tri, val));
              }
            else if (nz==0){//nz==0
              double val = dotproduct3D(dRnz, side1)*edge1factor+dotproduct3D(dRnz, side2)*edge2factor;
              tripletList.push_back(Eigen::Triplet<double>(pos, Tri, val)); 
              }
            else if (nx==ny){//Take derivative w.r.t xi 2
              double val = dotproduct3D(dRnx, side1)*edge1factor+dotproduct3D(dRnx, side2)*edge2factor;
              val += dotproduct3D(dRny, side1)*edge1factor+dotproduct3D(dRny, side2)*edge2factor;
              tripletList.push_back(Eigen::Triplet<double>(pos, Tri, val));          
            }
            else{
              double val = dotproduct3D(dRnz, side1)*edge1factor+dotproduct3D(dRnz, side2)*edge2factor;
              val-= dotproduct3D(dRny, side1)*edge1factor+dotproduct3D(dRny, side2)*edge2factor;
              tripletList.push_back(Eigen::Triplet<double>(pos, Tri, val)); 
            }
            }
            else{
              if (nx==0){
              double val = dotproduct3D(dRny, side1)*edge1factor+dotproduct3D(dRny, side2)*edge2factor;
              tripletList.push_back(Eigen::Triplet<double>(pos, Tri, val));
     
              }
              else{
              double val = -dotproduct3D(dRnz, side1)*edge1factor-dotproduct3D(dRnz, side2)*edge2factor;
              val-= dotproduct3D(dRnx, side1)*edge1factor+dotproduct3D(dRnx, side2)*edge2factor;
              tripletList.push_back(Eigen::Triplet<double>(pos, Tri, val)); 

              }
            }

          }
          else{

          double val1 = dotproduct3D(dRnx, side1)*edge1factor+dotproduct3D(dRnx, side2)*edge2factor;
          tripletList.push_back(Eigen::Triplet<double>(pos, Tri, val1)); 
          //Take derivative w.r.t xi 1 and xi 2
          double val2 = dotproduct3D(dRny, side1)*edge1factor+dotproduct3D(dRny, side2)*edge2factor;
          tripletList.push_back(Eigen::Triplet<double>(pos+1, Tri, val2));              
          }
      }
    }
}


template <int L>
void SavePosition(FILE* file, double (&rvec)[L+1][L+1][3]){
  for (int nx=0; nx<=L; nx++){
    for (int ny=0; ny<=L-nx; ny++){
      fprintf(file, "%d %d %.12f %.12f %.12f\n", nx, ny, rvec[ny][nx][0],rvec[ny][nx][1],rvec[ny][nx][2]); 
      }
    }
}


template <int L>
void SaveVertexFunction(FILE* file, double (&func)[L+1][L+1]){
  for (int nx=0; nx<=L; nx++){
    for (int ny=0; ny<=L-nx; ny++){
      fprintf(file, "%d %d %.12f \n", nx, ny, func[ny][nx]); 
      }
    }
}

template <int L>
void SaveMetric(FILE* file, double (&rvec)[L+1][L+1][3], int (&arr)[L*L][3]){
    for (int Tri=0; Tri<L*L; Tri++){
    double r[3][3] =  {};

      for (int vert = 0; vert<3; vert++){
          int s = arr[Tri][vert];
          int nx = s % (L+1); 
          int ny = s / (L+1);
          for (int mu=0; mu <3; mu++){
          r[vert][mu] = rvec[ny][nx][mu];
      }
      }
          double mside1[3] =  {0,0,0};
          double mside2[3] =  {0,0,0};
          double mside3[3] =  {0,0,0};

          for (int mu=0; mu<3; mu++){
            mside1[mu] = r[0][mu]-r[1][mu];
            mside2[mu] = r[0][mu]-r[2][mu];
            mside3[mu] = mside1[mu]-mside2[mu];
          } 

          double a = dotproduct3D(mside1, mside1);
          double b = dotproduct3D(mside2, mside2);
          double c = dotproduct3D(mside3, mside3);
          fprintf(file, "%d %d %d %.12f %.12f %.12f\n", arr[Tri][0],arr[Tri][1],arr[Tri][2], sqrt(a),  sqrt(c), sqrt(b));
          //Prints the traingle edge length in cyclic permutation to the order of the traingles. 
    }
}
template<int L>
void ReadOrbitData(const std::string& filePath, int (&orbitclass)[L+1][L+1][2], int classnumber, double (&target)[L+1][L+1]){
   std::ifstream file(filePath);
    std::string line;
  for (int classn=0; classn<classnumber; classn++){
    if (!getline(file, line)) {
          std::cerr << "Error reading line for class number" << classn << std::endl;
          return;
      }
    std::istringstream iss(line);
    double data; 
    iss>>data; 
  for (int nx = 0; nx<=L; nx++){
      for (int ny = 0; ny<=L-nx; ny++){
        int nz = L-nx-ny; 
       if (orbitclass[ny][nx][0]==classn){
        target[ny][nx] = data;  
       }
        }
      }
    }
}

template <int L>
void ReadBarycentric(const std::string& filePath, double (&xivec)[L+1][L+1][3], int (&orbitclass)[L+1][L+1][2], int classnumber){
    std::ifstream file(filePath);
    std::string line;
    int row = 0;
  for (int classn=0; classn<classnumber; classn++){
    if (!getline(file, line)) {
          std::cerr << "Error reading line for class number" << classn << std::endl;
          return;
      }
    std::istringstream iss(line);
    int orbitID; 
    iss>>orbitID; 
  for (int nx = 0; nx<=L; nx++){
      for (int ny = 0; ny<=L-nx; ny++){
        int nz = L-nx-ny; 
       if (orbitclass[ny][nx][0]==orbitID){
        int xi[3];
        xi[0] = nx;
        xi[1] = ny;
        xi[2] = L - nx - ny;
        std::sort(xi, xi + 3, std::greater<int>());
        double xi1, xi2; 
        iss>>xi1; 
        iss>>xi2; 
        xivec[xi[1]][xi[0]][0] = xi1; 
        xivec[xi[1]][xi[0]][1] = xi2; 
        xivec[xi[1]][xi[0]][2] = 1-xi1-xi2; 
        UpdateOrbit<L>(orbitclass, xivec, xi[0], xi[1]); 
        break; 
       }
        }
      }
    }
}

template <int L>
void SaveBarycentric(FILE* file, double (&xivec)[L+1][L+1][3], int (&orbitclass)[L+1][L+1][2], int classnumber){
  int classcounter[classnumber];
  for (int i = 0; i<classnumber; i++){
    classcounter[i]=0; 
  } 
  for (int nx=0; nx<=L; nx++){
    for (int ny=0; ny<=L-nx; ny++){
      int orbit_num = orbitclass[ny][nx][0]; 
      if (classcounter[orbit_num]==0){
      int xi[3]; 
      xi[0] = nx;
      xi[1] = ny;
      xi[2] = L - nx - ny;
      std::sort(xi, xi + 3, std::greater<int>());
        fprintf(file, "%d %.12f %.12f\n", orbit_num, xivec[xi[1]][xi[0]][0], xivec[xi[1]][xi[0]][1]); 
        classcounter[orbit_num] = 1; }
      }
    }
}
template <int L>
void EquivalenceClassID(std::unordered_map<std::string, int> &orbit_map, int (&orbitclass)[L+1][L+1][2], int dof){
  int classcounter[dof];
  for (int i =0; i<dof; i++){
    classcounter[i] = -1;
  } 
  for (int nx = 0; nx<=L; nx++){
      for (int ny = 0; ny<=L-nx; ny++){
        int nz = L-nx-ny; 
        int orbitID=orbitclass[ny][nx][0];
        if (classcounter[orbitID]==-1){
        int xi[3];
        xi[0] = nx;
        xi[1] = ny;
        xi[2] = L - nx - ny;
        std::sort(xi, xi + 3, std::greater<int>());
        std::string orbit_name = string_format("%d_%d_%d", xi[0], xi[1], xi[2]);
        orbit_map[orbit_name] = orbitID; 
        classcounter[orbitID] = 0; 
        }
      }
    }

}
template<int L> //The orbit class object is initialized to -1 
int BuildEquivalenceClass(int (&orbitclass)[L+1][L+1][2]){
  //Choose permutation of basis to be
  //  {(1,2,3), (1, 3, 2),  (2,1,3), (2, 3, 1), (3, 1, 2), (3, 2, 1)}
  int dof = 0;
  for (int nx = 0; nx<=L; nx++){
    for (int ny = 0; ny <=L-nx; ny++){
      int nz = L-nx-ny; 
      if (orbitclass[ny][nx][0] == -1){
        orbitclass[ny][nx][0] = dof; 
        orbitclass[ny][nx][1] = 0; //label permutation 

        orbitclass[nz][nx][0] = dof; 
        orbitclass[nz][nx][1] = 1; //label permutation 

        orbitclass[nx][ny][0] = dof; 
        orbitclass[nx][ny][1] = 2; //label permutation 

        orbitclass[nz][ny][0] = dof; 
        orbitclass[nz][ny][1] = 3; //label permutation 

        orbitclass[nx][nz][0] = dof; 
        orbitclass[nx][nz][1] = 4; //label permutation 

        orbitclass[ny][nz][0] = dof; 
        orbitclass[ny][nz][1] = 5; //label permutation 
        dof+=1; 
      }
    }
  }
  return dof; 
}

template<int L> //The orbit class object is initialized to -1 
void UpdateOrbit(int (&orbitclass)[L+1][L+1][2], double (&xivec)[L+1][L+1][3], int nx, int ny){
  double values[3] = {xivec[ny][nx][0],xivec[ny][nx][1],xivec[ny][nx][2]};
  int permutation[6][3] = {{1,2,3}, {1, 3, 2},  {2,1,3}, {2, 3, 1}, {3, 1, 2},{3, 2, 1}};
  int nz = L-nx-ny; 
  int position[3] = {nx, ny, nz}; 
  //printf("Main Starting Point %d %d %d with values %.4f %.4f %.4f------------\n", nx,ny,nz, values[0], values[1], values[2]);
  for (int permID = 0; permID<6; permID++){ 
    int new_nx = position[permutation[permID][0]-1]; 
    int new_ny = position[permutation[permID][1]-1];
    int new_nz = position[permutation[permID][2]-1];
    //printf("Target Permutation %d %d %d with ID %d \n", new_nx,new_ny,new_nz, xiID);
    //printf("Working on Permutation %d %d %d with ID %d \n", new_nx,new_ny,new_nz, xiID);
    xivec[new_ny][new_nx][0]=  values[permutation[permID][0]-1]; 
    xivec[new_ny][new_nx][1]=  values[permutation[permID][1]-1];
    xivec[new_ny][new_nx][2]=  values[permutation[permID][2]-1];
  }
}
template <int L>
void updateBarycentricwSymmetry(double (&xivec)[L+1][L+1][3], int (&vecInf)[L+1][L+1][3] , int (&classes)[L+1][L+1][2], Eigen::VectorXd& sol, int dof){
int classcounter[dof]; 
for (int i =0; i<dof; i++){
  classcounter[i] = -1; 
}
for (int ny =0; ny<=L; ny++){
        for (int nx =0; nx<=L-ny; nx++){
          int nz = L-nx-ny;
          int pos = vecInf[nx][ny][0];
          int classID = classes[ny][nx][0];
          if (pos==-1){continue;}
          else if (classcounter[classID]==1){continue;}//Orbit has already been updated
          else{
          //printf("Updating Site %d %d with original values %.4f %.4f %.4f\n", nx, ny, xivec[ny][nx][0],xivec[ny][nx][1],xivec[ny][nx][2]);
          //Update the primary site
          if (nx==0){//xi 2 is updated
              //printf("case 1\n");
                    xivec[ny][nx][1]+=sol[pos];
                    xivec[ny][nx][2]-=sol[pos];
              }
              else if (ny==0){//xi 1 is updated
              //printf("case 2\n");
                    xivec[ny][nx][0]+=sol[pos];
                    xivec[ny][nx][2]-=sol[pos];
              }
              else if (nz==0){//xi 1 is updated
              //printf("case 3\n");
                    xivec[ny][nx][0]+=sol[pos];
                    xivec[ny][nx][1]-=sol[pos];
              }

              else{
                if (nx==ny){
                    xivec[ny][nx][0]+=sol[pos];
                    xivec[ny][nx][1]+=sol[pos];
                    xivec[ny][nx][2]-=2*sol[pos];
                }
                else if (ny==nz){
                    xivec[ny][nx][0]-=2*sol[pos];
                    xivec[ny][nx][1]+=sol[pos];
                    xivec[ny][nx][2]+=sol[pos];
                }
                else if (nx==nz){
                    xivec[ny][nx][0]+=sol[pos];
                    xivec[ny][nx][1]-=2*sol[pos];
                    xivec[ny][nx][2]+=sol[pos];
                }
                //printf("case 4\n");
                else{
                xivec[ny][nx][0]+=sol[pos];
                xivec[ny][nx][1]+=sol[pos+1];
                xivec[ny][nx][2]-=(sol[pos]+sol[pos+1]);}
              }  

            //printf("Updated Site %d %d with new values %.4f %.4f %.4f\n", nx, ny, xivec[ny][nx][0],xivec[ny][nx][1],xivec[ny][nx][2]);
            UpdateOrbit<L>(classes, xivec, nx, ny);     
            classcounter[classID] +=1; 
          }
        }}
}




template <int L>
void setUpGMRESsolver(std::vector<T>& (tripletList),  Eigen::VectorXd& b, Eigen::VectorXd& sol,int dof, int iter_num){
    SpMat sparseMat(dof, L*L);
    sparseMat.setFromTriplets(tripletList.begin(), tripletList.end());
    sparseMat = -1*sparseMat; 
    SpMat AtA = sparseMat*(sparseMat.transpose());
    Eigen::VectorXd target = sparseMat*b;
    Eigen::GMRES<SpMat> solver; 
    solver.compute(AtA);
    solver.setMaxIterations(iter_num);
    sol = solver.solve(target);
}

template <int L>
void setUpStackedGMRESsolver(std::vector<T>& (tripletList),std::vector<T>& (tripletList2),  Eigen::VectorXd& b, Eigen::VectorXd& sol,int dof, int iter_num){
    SpMat A(dof, L*L);
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    //std::cout << "Original 1:\n" << Eigen::MatrixXd(A) << std::endl;

    SpMat B(dof, L*L);
    B.setFromTriplets(tripletList2.begin(), tripletList2.end());
    //std::cout << "Original 2:\n" << Eigen::MatrixXd(B) << std::endl;
    SpMat C(dof, 2*L*L); 
    C.reserve(A.nonZeros()+B.nonZeros());
    for (int k=0; k<A.outerSize(); ++k)
    for (SpMat::InnerIterator it(A,k); it; ++it)
        C.insert(it.row(), it.col()) = it.value();
    // Copy elements from B, offsetting the column index by y
    for (int k=0; k<B.outerSize(); ++k)
        for (SpMat::InnerIterator it(B,k); it; ++it)
            C.insert(it.row(), it.col() + L*L) = it.value();
    //Eigen::MatrixXd denseMat = Eigen::MatrixXd(C);
    //std::cout << "Combined:\n" << denseMat << std::endl;

    C = -1*C; 
    SpMat AtA = C*(C.transpose());
    Eigen::VectorXd target = C*b;
    Eigen::GMRES<SpMat> solver; 
    solver.compute(AtA);
    solver.setMaxIterations(iter_num);
    sol = solver.solve(target);
}


template <int L>
void PrintGeometry(int (&arr)[L*L][3], double (&rvec)[L+1][L+1][3]){
  Eigen::VectorXd a = CurrentAreaList<L>(arr, rvec);
  Eigen::VectorXd p = CurrentPerimeterList<L>(arr, rvec);
  Eigen::VectorXd cr = CurrentCircumradiusList<L>(arr, rvec); 
  double aRMS = returnRMS(a); 
  double pRMS = returnRMS(p); 
  double crRMS = returnRMS(cr); 
  printf(" %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f", 
  a.array().mean(), p.array().mean(), cr.array().mean(), 
  a.array().square().mean(), p.array().square().mean(), cr.array().square().mean(),
  aRMS, pRMS, crRMS  
  );
}

template <int L>
void PrintRenormalizedRMS(int (&arr)[L*L][3], double (&rvec)[L+1][L+1][3]){
  Eigen::VectorXd a = CurrentAreaList<L>(arr, rvec);
  Eigen::VectorXd p = CurrentPerimeterList<L>(arr, rvec);
  Eigen::VectorXd cr = CurrentCircumradiusList<L>(arr, rvec); 
  double aRMS = returnRMS(a)/(a.array().mean()*a.array().mean()); 
  double pRMS = returnRMS(p)/(p.array().mean()*p.array().mean());
  double crRMS = returnRMS(cr)/(cr.array().mean()*cr.array().mean());
  printf(" %.16f %.16f %.16f", 
  aRMS, pRMS, crRMS  
  );
}

template <int L>
bool CheckSymmetry(double (&xivec)[L+1][L+1][3],int (&vecInf)[L+1][L+1][3], int (&orbitclass)[L+1][L+1][2], int dof){
int classcounter[dof]; 
for (int i =0; i<dof; i++){
  classcounter[i] = -1; 
}
//printf("Printing Class Counter = (%d %d %d)\n", classcounter[0][0],classcounter[1][0],classcounter[2][0]);
int permutation[6][3] = {{1,2,3}, {1, 3, 2},  {2,1,3}, {2, 3, 1}, {3, 1, 2},{3, 2, 1}};
int countViolation = 0;  
  for (int nx = 0; nx<=L; nx++){
    for (int ny = 0; ny<=L-nx; ny++){
        int nz = L-nx-ny;
        int position[3] = {nx, ny, nz}; 
        int pos = vecInf[nx][ny][0];
        int classID = orbitclass[ny][nx][0];
        //printf("Refernce Point %d %d %d with classID: %d\n", nx, ny, nz, classID);
        if (pos==-1){continue;}
        else if (classcounter[classID]!=-1){continue;}
        else{
          double refValue[3] ={xivec[ny][nx][0],xivec[ny][nx][1],xivec[ny][nx][2]}; 
          classcounter[classID] = 1; 
          for (int permID = 0; permID<6; permID++){
           // printf("position Permutation ID: %d\n", permID); 
            int new_nx = position[permutation[permID][0]-1]; 
            int new_ny = position[permutation[permID][1]-1];
            int new_nz = position[permutation[permID][2]-1];
            //printf("Permutation ID: %d\n", xiID); 
            //printf("Reference Position v.s. new Position: (%d %d %d) and (%d %d %d)\n)", nx, ny, nz, new_nx, new_ny, new_nz);
            if (xivec[new_ny][new_nx][0]==refValue[permutation[permID][0]-1]&&
                xivec[new_ny][new_nx][1]==refValue[permutation[permID][1]-1]&&
                xivec[new_ny][new_nx][2]==refValue[permutation[permID][2]-1]){
                continue;
              }
            else{
              printf("Fails with target class %d %d %d %.4f %.4f %.4f\n", nx, ny, nz, refValue[0],refValue[1],refValue[2]);
              printf("Transformation %d %d %d\n", permutation[permID][0],permutation[permID][1],permutation[permID][2]);
              printf("permutation class %d %d %d %.4f %.4f %.4f\n", new_nx, new_ny, new_nz, xivec[new_ny][new_nx][0],xivec[new_ny][new_nx][1],xivec[new_ny][new_nx][2]);
              printf("Transformed Ref Values %.4f %.4f %.4f\n", refValue[permutation[permID][0]-1],refValue[permutation[permID][1]-1],refValue[permutation[permID][2]-1]);
              countViolation+=1; 
            }
          
          }
        }
    }
  }
if (countViolation != 0){return false;}
else{return true;}
}


template <int L>
void ReadRData(const std::string& filePath, double (&rvec)[L+1][L+1][3]){
    std::ifstream file(filePath);
    std::string line;
    int row = 0;
    for (int ny = 0; ny<=L; ny++){
      for (int nx = 0; nx<=L-ny; nx++){
          if (!getline(file, line)) {
                std::cerr << "Error reading line for ny=" << ny << ", nx=" << nx << std::endl;
                return;
            }
        std::istringstream iss(line);
        double value;
        for (int i = 0; i < 2; ++i) {
            iss >> value; 
        }
        int col = 0;
        while (iss >> value && col < 3) {
            rvec[ny][nx][col++] = value;
        }
      }
    }
}

template<int L>
Eigen::VectorXd returnGradient(std::vector<T>& tripletList, Eigen::VectorXd& target, int dof){
    SpMat sparseMat(dof, L*L);
    sparseMat.setFromTriplets(tripletList.begin(), tripletList.end()); 
    Eigen::VectorXd Gradient = 2*(sparseMat*target).rowwise().sum(); 
    return Gradient; 
}
