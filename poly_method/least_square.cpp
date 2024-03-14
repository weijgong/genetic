#include<iostream>
#include<math.h>
#include<vector>
#include <fstream>
#include <string>

using namespace std;

void gaussEliminationLS(int m, int n, vector<vector<double>>&a, vector<double>& x);
void printMatrix(int m, int n, vector<vector<double>> matrix);
vector<double> LeastSquarePoly(vector<double>x, vector<double>y, int n);

void gaussEliminationLS(int m, int n, vector<vector<double>>&a, vector<double>& x){
    int i,j,k;
    for(i=0;i<m-1;i++){
        for(k=i+1;k<m;k++){
            if(fabs(a[i][i])<fabs(a[k][i])){
                for(j=0;j<n;j++){
                    double temp;
                    temp=a[i][j];
                    a[i][j]=a[k][j];
                    a[k][j]=temp;
                }
            }
        }
        //Begin Gauss Elimination
        for(k=i+1;k<m;k++){
            double  term=a[k][i]/ a[i][i];
            for(j=0;j<n;j++){
                a[k][j]=a[k][j]-term*a[i][j];
            }
        }
         
    }
    //Begin Back-substitution
    for(i=m-1;i>=0;i--){
        x[i]=a[i][n-1];
        for(j=i+1;j<n-1;j++){
            x[i]=x[i]-a[i][j]*x[j];
        }
        x[i]=x[i]/a[i][i];
    }
             
}

void printMatrix(int m, int n, vector<vector<double>> matrix){
    int i,j;
    for(i=0;i<m;i++){
        for(j=0;j<n;j++){
            printf("%lf\t",matrix[i][j]);
        }
        printf("\n");
    } 
}

// x自变量,y因变量,n拟合多项式的阶数
vector<double> LeastSquarePoly(vector<double>x, vector<double>y, int n){
    int N = x.size();
    vector<double>X(2*n+1,0);
    for(int i=0;i<=2*n;i++)
        for(int j=0;j<N;j++)
            X[i]=X[i]+pow(x[j],i);
    vector<vector<double>> B(n+1,vector<double>(n+2)); 
    vector<double> Y(n+1,0);      
        for(int i=0;i<=n;i++)
            for(int j=0;j<N;j++)
                Y[i]=Y[i]+pow(x[j],i)*y[j];
        
        for(int i=0;i<=n;i++)
            for(int j=0;j<=n;j++)
                B[i][j]=X[i+j];
                
        for(int i=0;i<=n;i++)
            B[i][n+1]=Y[i];
        
        vector<double> A(n+1,0);
        // printMatrix(n+1,n+2,B);
        gaussEliminationLS(n+1,n+2,B,A);
    // A:shape(n+1,)
    return A;
}

vector<double> poly_main(vector<double> x_arr,int poly_n){
    int N =x_arr.size();
    vector<double> x(N,0);
    for(int i=0;i<N;i++)x[i]=i;
    // int poly_n = 0;
    vector<double> a = LeastSquarePoly(x,x_arr,poly_n);
    // for(auto i:a)cout<<i<<"\t";
    // cout<<endl;
    return a;
}

vector<double> slice_vector(vector<double> arr,int slice_index){
    vector<double>new_arr(slice_index);
    for(int i=0;i<slice_index;i++){
        new_arr[i]=arr[i];
    }
    return new_arr;
}

int main(){
    const int POLY_MAX      = 1500;
    int       DATA_LINE_NUM = 0;
    int       POLY_NUM      = 10;

    ifstream data_fsile ("/home/gwj/genetic/sat_data/track_simu_yaogan35.csv");

    vector<string>time_arr;
    vector<double>x_arr;
    vector<double>y_arr;
    vector<double>z_arr;
    vector<double>vx_arr;
    vector<double>vy_arr;
    vector<double>vz_arr;
    if (data_fsile.is_open())
    {
        string line;
        while ( getline (data_fsile,line) )
        {
            if(line.size()<10)break;
            DATA_LINE_NUM++;
            int place = line.find(',');
            time_arr.push_back(line.substr(0,place));
            line = line.substr(place+1,line.size());
            place = line.find(',');
            x_arr.push_back(stoi(line.substr(0,place)));
            line = line.substr(place+1,line.size());
            place = line.find(',');
            y_arr.push_back(stoi(line.substr(0,place)));
            line = line.substr(place+1,line.size());
            place = line.find(',');
            z_arr.push_back(stoi(line.substr(0,place)));
            line = line.substr(place+1,line.size());
            place = line.find(',');
            vx_arr.push_back(stoi(line.substr(0,place)));
            line = line.substr(place+1,line.size());
            place = line.find(',');
            vy_arr.push_back(stod(line.substr(0,place)));
            line = line.substr(place+1,line.size());
            place = line.find(',');
            vz_arr.push_back(stod(line.substr(0,place)));
            // cout<<line<<"\n";
        }
        data_fsile.close();
    }
    else cout << "Unable to open file"; 
    
    vector<double> a_x;
    vector<double> a_y;
    vector<double> a_z;
    vector<double> a_vx;
    vector<double> a_vy;
    vector<double> a_vz;

    vector<double> p_x =slice_vector(x_arr ,POLY_MAX);
    vector<double> p_y =slice_vector(y_arr ,POLY_MAX);
    vector<double> p_z =slice_vector(z_arr ,POLY_MAX);
    vector<double> p_vx=slice_vector(vx_arr,POLY_MAX);
    vector<double> p_vy=slice_vector(vy_arr,POLY_MAX);
    vector<double> p_vz=slice_vector(vz_arr,POLY_MAX);

    a_x  = poly_main(p_x ,POLY_NUM);
    a_y  = poly_main(p_y ,POLY_NUM);
    a_z  = poly_main(p_z ,POLY_NUM);
    a_vx = poly_main(p_vx,POLY_NUM);
    a_vy = poly_main(p_vy,POLY_NUM);
    a_vz = poly_main(p_vz,POLY_NUM);

    vector<double> pre_x(DATA_LINE_NUM-POLY_MAX+1,0);

    for(int x=0;x<DATA_LINE_NUM;x++){
        for(int i=0;i<a_x.size();i++){
            pre_x[x]+=a_x[i]*pow(x,i);            
        }
        cout<<pre_x[x]<<endl;
    }
}