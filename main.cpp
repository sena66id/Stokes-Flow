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
    //export_vtu("test.vtu", x, element, C1);

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


    MatrixXd K(x.size()*3, x.size()*3);
    VectorXd U(x.size()*3);
    VectorXd R(x.size()*3);
    U = VectorXd::Zero(x.size()*3);
    R = VectorXd::Zero(x.size()*3);


    for(int i=0; i<x.size()*3; i++){
      for(int j=0; j<x.size()*3; j++){
        K(i,j) = 0;
      }
    }

    

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

       
     

      calcInverseMatrix_2x2(inv_dxdr,dxdr);

     double dndx[2][3];

   

     for(int i=0;i<2;i++){
            for(int j=0;j<3;j++){
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
        r1[1][i]=dndx[1][i];
    }

    for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                for(int k=0;k<2;k++){
                    Ke1[i][j]+=r1[i][k]*r2[k][j];
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
                for(int k=0;k<1;k++){
                    Ke2[i][j]+=r3[i][k]*r4[k][j];
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
                for(int k=0;k<1;k++){
                    Ke3[i][j]+=r5[i][k]*r6[k][j];
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



    MatrixXd Ke(9, 9);
    Ke = VectorXd::Zero(9,9);

    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            Ke(i,j)+=Ke1[i][j];
        }
    }
     for(int i=0;i<3;i++){
        for(int j=6;j<9;j++){
            Ke(i,j)+=Ke2[i][j-6];
        }
    }
     for(int i=3;i<6;i++){
        for(int j=3;j<6;j++){
            Ke(i,j)+=Ke1[i-3][j-3];
        }
    }
     for(int i=3;i<6;i++){
        for(int j=6;j<9;j++){
            Ke(i,j)+=Ke3[i-3][j-6];
        }
    }
     for(int i=6;i<9;i++){
        for(int j=0;j<3;j++){
            Ke(i,j)+=Ke4[i-6][j];
        }
    }
     for(int i=6;i<9;i++){
        for(int j=3;j<6;j++){
            Ke(i,j)+=Ke5[i-6][j-3];
        }
    }

    cout<<Ke<<endl;


    }

return 0;


}


double calcDeterminant_2x2( const double (&a)[2][2])
{
  double det  = a[0][0] * a[1][1]-a[1][0]*a[0][1];
  return det;
}

void calcInverseMatrix_2x2(double (&inv_a)[2][2],const double (&a)[2][2])
{
  double det;

  det =1.0/ calcDeterminant_2x2(a);

  inv_a[0][0] = a[1][1];
  inv_a[0][1] = -1.0*a[0][1];
  inv_a[1][0] =-1.0*a[1][0];
  inv_a[1][1] = a[0][0];
  

  for(int i=0;i<2;i++){
    for(int j=0;j<2;j++) inv_a[i][j] = inv_a[i][j] / det;
}
}
   // vector<vector<double>> K(





    




