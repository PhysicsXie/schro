#include<iostream>
#include<eigen3/Eigen/Dense>
#include<cmath>
#include<fstream>
#include<string>

using namespace std;
using namespace Eigen;


int main()
{
    //setting const
    const int N=500;
    const double pi=3.1415926535;
    const double xmin=0.0;
    const double xmax=1; // unit A
    const double m=511; // electron mass, unit keV
    const double dx=(xmax-xmin)/(N-1); //step of x
    const double hb=1.9733*1.9733;

    //defining spacial x
    VectorXd x(N);

    for(int i=0;i<N;++i){
        x(i)=xmin+i*dx;
    }

    //defining kinetic matrix
    MatrixXd H = MatrixXd::Zero(N,N);

    //filling the kinetic
    for(int i=0;i<N;++i){

        if(i==0){
            H(i,i)=-2.0;
            H(i,i+1)=1.0;
        }
        else if (i==N-1)
        {
            H(i,i)=-2.0;
            H(i,i-1)=1;
        }
        else{
            H(i,i-1)=1.0;
            H(i,i+1)=1.0;
            H(i,i)=-2.0;
        }
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

    eigenvalues=eigenvalues*hb;
    eigenvectors=eigenvectors*(1.0/sqrt(dx));

//    MatrixXd normalized_eigenvectors=eigenvectors.colwise().normalized();

    ofstream eigenval("eigenvalues.csv");
    if(!eigenval.is_open()){
        cerr<<"Failed to open eigenvalues.csv!!"<<endl;
        return 1;
    }
    
    for(int i=0;i<=10;++i){

        eigenval<<eigenvalues(i)<<endl;
    }

    eigenval.close();

    ofstream eigenfun("eigenfunction.csv");
        if(!eigenfun.is_open()){
        cerr<<"Failed to open eigenfunction.csv!!"<<endl;
        return 1;
    }

    for(int i=0;i<N;++i){

        eigenfun<<x(i)<<","<<eigenvectors(i,0)<<","<<eigenvectors(i,1)<<endl;
    }

    eigenfun.close();

    ofstream potential("potential.csv");
    
        if (!potential.is_open())
        {
            cerr<<"Failed to open potential.csv!!"<<endl;
            return 1;
        }
    
    for(int i=0;i<N;++i){
        potential<<x(i)<<","<<hb*V(i)<<endl;
    }
    potential.close();
    

    system("pause");
    return 0;
}
