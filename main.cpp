#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<sstream>
#include<cmath>
#include "Eigen/Core"
#include "Eigen/LU"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"
double calcDeterminant_2x2(const double (&a)[2][2]);
void calcInverseMatrix_2x2(double (&inv_a)[2][2],const double (&a)[2][2]);

using namespace Eigen;
using namespace std;
typedef Triplet<double> T;
void export_vtu(const std::string &file, vector<vector<double>> node, vector<vector<int>> element, vector<double> C)
{
    FILE *fp;
  if ((fp = fopen(file.c_str(), "w")) == NULL)
  {
    cout << file << " open error" << endl;
    exit(1);
  }
  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", node.size(), element.size());
  fprintf(fp, "<Points>\n");
  int offset = 0;
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * node.size() * 3;
  fprintf(fp, "</Points>\n");
  fprintf(fp, "<Cells>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
  for (int i = 0; i < element.size(); i++){
    for (int j = 0; j < element[i].size(); j++) fprintf(fp, "%d ", element[i][j]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
  int num = 0;
  for (int i = 0; i < element.size(); i++)
  {
    num += element[i].size();
    fprintf(fp, "%d\n", num);
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (int i = 0; i < element.size(); i++) fprintf(fp, "%d\n", 5);
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Cells>\n");
  fprintf(fp, "<PointData>\n");
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"pressure[Pa]\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * node.size();
  fprintf(fp, "</PointData>\n");
  fprintf(fp, "<CellData>\n");
  fprintf(fp, "</CellData>\n");
  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</UnstructuredGrid>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fclose(fp);
  fstream ofs;
  ofs.open(file.c_str(), ios::out | ios::app | ios_base::binary);
  double *data_d = new double[node.size()*3];
  num = 0;
  int size=0;
  for (int ic = 0; ic < node.size(); ic++){
    data_d[num] = node[ic][0];
    num++;
    data_d[num] = node[ic][1];
    num++;
    data_d[num] = 0.0;
    num++;
  }
  size=sizeof(double)*node.size()*3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_d, size);
  num=0;
  for (int ic = 0; ic < node.size(); ic++){
      data_d[num]   = C[ic];
      num++;
  }
  size=sizeof(double)*node.size();
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_d, size);
  delete data_d;
  ofs.close();
  if ((fp = fopen(file.c_str(), "a")) == NULL)
  {
    cout << file << " open error" << endl;
    exit(1);
  }
  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

   
int main()
{
    string str;
    ifstream ifs("node.dat");
    //vector<double> t(2) は double t[2]と同じ
    //vector<vector<double>> t(2, vector<double>(2))は double t[2][2]と同じ
    vector<vector<double>> x;
    while(getline(ifs,str)){
        istringstream ss(str);
        string tmp;
        vector<double> tmp_x;
        for(int j=0; j<3; j++){
            getline(ss, tmp, ' ');
            tmp_x.push_back(stod(tmp));
        }
        x.push_back(tmp_x);
    }
    ifs.close();
    ifs.open("element.dat");
    vector<vector<int>> element;
    while(getline(ifs,str)){
        istringstream ss(str);
        string tmp;
        vector<int> tmp_element;
        for(int j=0; j<4; j++){
            getline(ss, tmp, ' ');
            if(j==0) continue;
            tmp_element.push_back(stoi(tmp));
        }
        element.push_back(tmp_element);
    }
    ifs.close();
    
    
    vector<double> C1(x.size(),0.0);
    double minimum = 100000.0;
    for(int i=0; i<x.size(); i++){
        minimum=min(minimum,x[i][0]);
    }
    ofstream ofs1("boundary_left.dat");
    for(int i=0; i<x.size(); i++){
        if(fabs(x[i][0]-minimum)<0.000001){
            ofs1 << i << endl;
            C1[i] = 1.0;
        }
    }
    ofs1.close();

    vector<double> C2(x.size(),0.0);
    double maximum = 0.001;
    for(int i=0; i<x.size(); i++){
        maximum=max(maximum,x[i][0]);
    }
    ofstream ofs2("boundary_right.dat");
    for(int i=0; i<x.size(); i++){
        if(fabs(x[i][0]-maximum)<0.000001){
            ofs2 << i << endl;
            C2[i] = 1.0;
        }
    }
    ofs2.close();

    vector<double> C3(x.size(),0.0);
    double minimum2 = 100000.0;
    for(int i=0; i<x.size(); i++){
        minimum2=min(minimum2,x[i][1]);
    }
    ofstream ofs3("boundary_bottom.dat");
    for(int i=0; i<x.size(); i++){
        if(fabs(x[i][1]-minimum2)<0.000001){
            ofs3 << i << endl;
            C3[i] = 1.0;
        }
    }
    ofs3.close();

    vector<double> C4(x.size(),0.0);
    double maximum2 = 0.001;
    for(int i=0; i<x.size(); i++){
        maximum2=max(maximum2,x[i][1]);
    }
    ofstream ofs4("boundary_upper.dat");
    for(int i=0; i<x.size(); i++){
        if(fabs(x[i][1]-maximum2)<0.000001){
            ofs4 << i << endl;
            C4[i] = 1.0;
        }
    }
    ofs4.close();



    
   

    ifs.open("boundary_left.dat");
    vector<int> boundary_left;
    while(getline(ifs,str)){
        istringstream ss(str);
        string tmp;
        getline(ss, tmp, ' ');
        boundary_left.push_back(stoi(tmp));
    }
    ifs.close();

    ifs.open("boundary_right.dat");
    vector<int> boundary_right;
    while(getline(ifs,str)){
        istringstream ss(str);
        string tmp;
        getline(ss, tmp, ' ');
        boundary_right.push_back(stoi(tmp));
    }
    ifs.close();

    ifs.open("boundary_bottom.dat");
    vector<int> boundary_bottom;
    while(getline(ifs,str)){
        istringstream ss(str);
        string tmp;
        getline(ss, tmp, ' ');
        boundary_bottom.push_back(stoi(tmp));
    }
    ifs.close();

    ifs.open("boundary_upper.dat");
    vector<int> boundary_upper;
    while(getline(ifs,str)){
        istringstream ss(str);
        string tmp;
        getline(ss, tmp, ' ');
        boundary_upper.push_back(stoi(tmp));
    }
    ifs.close();

    

    MatrixXd K(x.size()*3,x.size()*3);
    MatrixXd Ke(9,9);
    MatrixXd Ke6(3,3);
    VectorXd U(x.size()*3);
    VectorXd R(x.size()*3);
    U= VectorXd::Zero(x.size()*3);
    K = MatrixXd::Zero(x.size()*3,x.size()*3);
    R = VectorXd::Zero(x.size()*3);

    
    




    

    double dNdr[2][3]; 
    dNdr[0][0] = -1.0; dNdr[0][1] = 1.0; dNdr[0][2] = 0.0; 
    dNdr[1][0] = -1.0; dNdr[1][1] = 0.0; dNdr[1][2] = 1.0; 
    

    double dxdr[2][2];

   
    
    for(int ic=0; ic<element.size(); ic++){
        
     
      for(int i=0; i<2; i++){
        for(int j=0; j<2; j++){
        dxdr[i][j]=0.0;
          for(int k=0; k<3; k++){
            dxdr[i][j] += dNdr[i][k]*x[element[ic][k]][j];
          }
        }
      }

      


      double inv_dxdr[2][2];
      double det;


      det = calcDeterminant_2x2(dxdr);

     

     calcInverseMatrix_2x2(inv_dxdr,dxdr);

      
      
        
      
     double dndx[2][3];

     /*for(int i=0; i<2; i++){
        for(int j=0; j<2; j++){
            cout<<dxdr[i][j]<< " " ;
        }
     }
     cout << endl;*/

    
      //cout<<det<<endl;
   

     for(int i=0;i<2;i++){
            for(int j=0;j<3;j++){
                dndx[i][j]=0;
                for(int k=0;k<2;k++){
                    dndx[i][j]+=inv_dxdr[i][k]*dNdr[k][j];
                }
            }
     }
    
     

        double Ke1[3][3];
        double r1[3][2];
        double r2[2][3];

        for(int i=0;i<3;i++){
            r1[i][0]=dndx[0][i];
        }
        for(int i=0;i<3;i++){
            r1[i][1]=dndx[1][i];
        }
        for(int i=0;i<3;i++){
            r2[0][i]=dndx[0][i];
        }
        for(int i=0;i<3;i++){
            r2[1][i]=dndx[1][i];
        }

        for(int i=0;i<3;i++){
                for(int j=0;j<3;j++){
                    Ke1[i][j]=0.0;
                    for(int k=0;k<2;k++){
                        Ke1[i][j]-=r1[i][k]*r2[k][j]*0.0035*det*0.5;
                    }
                }
        }

        double Ke2[3][3];
        double r3[3][1];
        double r4[1][3];
        for(int i=0;i<3;i++){
            r3[i][0]=dndx[0][i];
        }
        for(int i=0;i<3;i++){
            r4[0][i]=0.333;
        }

        for(int i=0;i<3;i++){
                for(int j=0;j<3;j++){
                    Ke2[i][j]=0.0;
                    for(int k=0;k<1;k++){
                        Ke2[i][j]+=r3[i][k]*r4[k][j]*det*0.5;
                    }
                }
        }

        double Ke3[3][3];
        double r5[3][1];
        double r6[1][3];
        for(int i=0;i<3;i++){
            r5[i][0]=dndx[1][i];
        }
        for(int i=0;i<3;i++){
            r6[0][i]=0.333;
        }

        for(int i=0;i<3;i++){
                for(int j=0;j<3;j++){
                    Ke3[i][j]=0;
                    for(int k=0;k<1;k++){
                        Ke3[i][j]+=r5[i][k]*r6[k][j]*det*0.5;
                    }
                }
        }

        double Ke4[3][3];
        double Ke5[3][3];

        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                Ke4[i][j]=Ke2[j][i];
                Ke5[i][j]=Ke3[j][i];
            }
        }


        

        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                Ke(i,j)=0.0;
                Ke(i,j)+=Ke1[i][j];
            }
        }
        for(int i=0;i<3;i++){
            for(int j=6;j<9;j++){
                 Ke(i,j)=0.0;
                Ke(i,j)+=Ke2[i][j-6];
            }
        }
        for(int i=3;i<6;i++){
            for(int j=3;j<6;j++){
                 Ke(i,j)=0.0;
                Ke(i,j)+=Ke1[i-3][j-3];
            }
        }
        for(int i=3;i<6;i++){
            for(int j=6;j<9;j++){
                 Ke(i,j)=0.0;
                Ke(i,j)+=Ke3[i-3][j-6];
            }
        }
        for(int i=6;i<9;i++){
            for(int j=0;j<3;j++){
                Ke(i,j)=0.0;
                Ke(i,j)+=Ke4[i-6][j];
            }
        }
        for(int i=6;i<9;i++){
            for(int j=3;j<6;j++){
                Ke(i,j)=0.0;
                Ke(i,j)+=Ke5[i-6][j-3];
            }
        }

      double dx = x[element[ic][0]][0] - x[element[ic][1]][0];
      double dy = x[element[ic][0]][1] - x[element[ic][1]][1];
      double h = sqrt(dx * dx + dy * dy);
      double tau = h * h / (4e0 * 0.0035) / 3e0;
      
      
      for (int p = 0; p <3; p++)
      {
        for (int q = 0; q < 3; q++)
        {
           Ke6(p,q)=0;
          for (int i = 0; i < 2; i++)
          {
            Ke6(p, q) -= tau * dndx[p][i] * dndx[q][i] * det *0.5;
          }
          
          Ke(p + 6, q + 6) += Ke6(p, q);
        }
      }
       

        

        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
               
                K(element[ic][i],element[ic][j]) += Ke(i,j);
            }

        }
        for(int i=0;i<3;i++){
            for(int j=3;j<6;j++){
                
                K(element[ic][i],element[ic][j-3]+x.size()) += Ke(i,j);
            }
        }
        for(int i=0;i<3;i++){
            for(int j=6;j<9;j++){
             
                K(element[ic][i],element[ic][j-6]+2*x.size()) += Ke(i,j);
            }
        }
        for(int i=3;i<6;i++){
            for(int j=0;j<3;j++){
                 
                K(element[ic][i-3]+x.size(),element[ic][j])+= Ke(i,j);
            }

        }
        for(int i=3;i<6;i++){
            for(int j=3;j<6;j++){
                
                K(element[ic][i-3]+x.size(),element[ic][j-3]+x.size()) += Ke(i,j);
            }
        }
        for(int i=3;i<6;i++){
            for(int j=6;j<9;j++){
                
                K(element[ic][i-3]+x.size(),element[ic][j-6]+2*x.size()) +=Ke(i,j);
            }
        }
        for(int i=6;i<9;i++){
            for(int j=0;j<3;j++){
                
                K(element[ic][i-6]+2*x.size(),element[ic][j]) += Ke(i,j);
            }

        }
        for(int i=6;i<9;i++){
            for(int j=3;j<6;j++){
                 
                K(element[ic][i-6]+2*x.size(),element[ic][j-3]+x.size()) += Ke(i,j);
            }
        }
        for(int i=6;i<9;i++){
            for(int j=6;j<9;j++){
               
                K(element[ic][i-6]+2*x.size(),element[ic][j-6]+2*x.size()) +=Ke(i,j);
            }
        }
       
    
    }

   for(int i=0; i<boundary_left.size(); i++){
      for(int j=0; j<x.size()*3; j++){
        K(boundary_left[i],j) = 0.0;
      }
     }
    /*for(int i=0; i<boundary_right.size(); i++){
      for(int j=0; j<x.size()*3; j++){
        K(boundary_right[i]+2*x.size(),j) = 0;
      }
    }*/
    for(int i=0; i<boundary_bottom.size(); i++){
      for(int j=0; j<x.size()*3; j++){
        K(boundary_bottom[i],j) = 0.0;
        K(boundary_bottom[i]+x.size(),j) = 0.0;
      }
    }
    for(int i=0; i<boundary_upper.size(); i++){
      for(int j=0; j<x.size()*3; j++){
        K(boundary_upper[i],j) = 0.0;
        K(boundary_upper[i]+x.size(),j) = 0.0;
      }
    }
    for(int i=0; i<boundary_left.size(); i++){
        K(boundary_left[i],boundary_left[i]) = 1.0;
      }
     /* for(int i=0; i<boundary_right.size(); i++){
        K(boundary_right[i]+2*x.size(),boundary_right[i]+2*x.size()) = 1.0;
      }*/
      for(int i=0; i<boundary_bottom.size(); i++){
        K(boundary_bottom[i],boundary_bottom[i]) = 1.0;
        K(boundary_bottom[i]+x.size(),boundary_bottom[i]+x.size()) = 1.0;
      }
       for(int i=0; i<boundary_upper.size(); i++){
        K(boundary_upper[i],boundary_upper[i]) = 1.0;
        K(boundary_upper[i]+x.size(),boundary_upper[i]+x.size()) = 1.0;
      }
    for(int i=0; i<boundary_right.size(); i++){
        R(boundary_left[i]) = 1.0;
        
     }

    SparseMatrix<double> K_sparse(x.size()*3, x.size()*3);
      vector<T> tripletList;
      for (int i=0;i<x.size()*3 ;i++){
        for(int j=0;j<x.size()*3 ;j++){
          if(K(i,j)!=0){
            tripletList.push_back(T(i,j,K(i,j)));
          }

        }
      }
      K_sparse.setFromTriplets(tripletList.begin(),tripletList.end());

     SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solver;
    solver.compute(K_sparse);
    if(solver.info()!=Success) {
    // decomposition failed
    cout << "decomposition failed" << endl;
    return -1;
    }
    U = solver.solve(R);
    if(solver.info()!=Success) {
    // solving failed
    std::cout << "solving failed" << endl;
    return -1;
    }
    ofstream outputfile("U.dat");
    outputfile << U;
    outputfile.close();

    vector<double> U1(x.size()*3);


    for(int i=0; i<x.size(); i++){
        U1[i]=U[i];
    }

    export_vtu("result.vtu", x, element, U1);




      

    


}







double calcDeterminant_2x2( const double (&a)[2][2])
{
  double det  = a[0][0] * a[1][1]-a[1][0]*a[0][1];
  return det;
}

void calcInverseMatrix_2x2(double (&inv_a)[2][2],const double (&a)[2][2])
{
  double k;

  k =1.0/ calcDeterminant_2x2(a);

  inv_a[0][0] = a[1][1];
  inv_a[0][1] = -1.0*a[0][1];
  inv_a[1][0] =-1.0*a[1][0];
  inv_a[1][1] = a[0][0];
  

  for(int i=0;i<2;i++){
    for(int j=0;j<2;j++) inv_a[i][j] = inv_a[i][j] * k;
}
}
   // vector<vector<double>> K(





    




