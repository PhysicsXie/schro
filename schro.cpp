#include<iostream>
#include<Eigen/Dense>
#include<cmath>
#include<fstream>
#include<string>

using namespace std;
using namespace Eigen;


int main()
{
    //setting const
    const int N=100;
    const double pi=3.1415926535;
    const double xmin=0.0;
    const double xmax=1.0; // unit \circ{A}
    const double m=511.0; // electron mass, unit keV
    const double dx=(xmax-xmin)/(N-1); //step of x

    //defining spacial x
    VectorXd x(N);

    for(int i=0;i<N;++i){
        x(i)=xmin+i*dx;
    }

    //defining kinetic matrix
    MatrixXd H = MatrixXd::Zero(N,N);

    //filling the kinetic
    for(int i=1;i<N-1;++i){
        H(i,i-1)=1.0;
        H(i,i+1)=1.0;
        H(i,i)=-2.0;
    }

    H=(-1/(2*m*dx*dx))*H;

    //defining potential matrix
    MatrixXd V = VectorXd::Zero(N);
    for(int i=0;i<N;++i){

        V(i)=500*(1/(2*m))*sin(pi*x(i));
        H(i,i)+=V(i);
    }

    //solving the eigen equation
    SelfAdjointEigenSolver<MatrixXd> solver(H);
    VectorXd eigenvalues=solver.eigenvalues();
    MatrixXd eigenvectors = solver.eigenvectors();

    ofstream eigenval("eigenvalues.csv");
    if(!eigenval.is_open()){
        cerr<<"Failed to open eigenvalues.csv!!"<<endl;
        return 1;
    }
    
    for(int i=0;i<eigenvalues.size();++i){

        eigenval<<eigenvalues(i)<<endl;
    }

    eigenval.close();

    ofstream eigenfun("eigenfunction.csv");
        if(!eigenfun.is_open()){
        cerr<<"Failed to open eigenfunction.csv!!"<<endl;
        return 1;
    }

    for(int i=0;i<N;++i){

        eigenfun<<x(i)<<","<<eigenvectors(i,0)<<","<<eigenvectors(i,1)<<","<<eigenvectors(i,2)<<endl;
    }

    eigenfun.close();

    system("pause");
    return 0;
}
