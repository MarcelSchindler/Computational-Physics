
#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <random>

using namespace std;


void leapfrog(double array_p_phi[], int N_md, double J, double h, int N)
{
    double epsilon= 1.0/N_md;

    array_p_phi[1] = array_p_phi[1] + epsilon*array_p_phi[0]/2;
    for (int i = 0; i < N_md-1; i++)
    {
    array_p_phi[0]= array_p_phi[0] - epsilon*(array_p_phi[1]/J-N*tanh(h+array_p_phi[1]));
    array_p_phi[1] = array_p_phi[1] + epsilon*array_p_phi[0];
    }
    array_p_phi[0]= array_p_phi[0] - epsilon*(array_p_phi[1]/J-N*tanh(h+array_p_phi[1]));
    array_p_phi[1] = array_p_phi[1] + epsilon*array_p_phi[0]/2;
}

double Hamiltionian(double p,double phi, double J, int N, double h)
{
    return p*p/2+phi*phi/(2*J)-N*log(2*cosh(h+phi));
}

double magnitization(double phi, double h)
{
    return tanh(h+phi);
}

double energy(double phi, double h, int N, double J)
{
    double beta =1;
    return -1/N-phi*phi/(2*J*N)-h*tanh(h+phi);
}

int factorial(int n)
{
    if(n > 1)
        return n * factorial(n - 1);
    else
        return 1;
}

double function(double J, double h, double x)
{
    return exp(0.5*J*x*x+h*x);
}
double Z_analytic(double J, double h, int N)
{
    double sum=0;
    for (int i = 0; i < N; i++)
    {
        sum+= (double)factorial(N)/(factorial(i)*factorial(N-i))*function(J,h,N-2*i);
        
    }
    return sum;
}
double magnitization_analytic(double J, double h, int N)
{
    double sum = 0;
    for (int i = 0; i < N; i++)
    {
        sum+=(double)factorial(N)/(factorial(i)*factorial(N-i))*(N-2*i)*function(J,h,N-2*i);
    }
    sum = sum/(N*Z_analytic(J,h,N));
    return sum;
    
}

double energy_analytic(double J, double h, int N)
{
    double sum = 0;
    for (int i = 0; i < N; i++)
    {
        sum+=(double)factorial(N)/(factorial(i)*factorial(N-i))*(0.5*J*(N-2*i)*(N-2*i)+h*(N-2*i))*function(J,h,N-2*i);
    }
    sum = -sum/(N*Z_analytic(J,h,N));
    return sum;
}

double standard_deviation(int R, int start_value, double mu, double array[])
{
    double deviation=0;
    for (int i = start_value-1; i < R; i++)
    {
        deviation+=(array[i]-mu)*(array[i]-mu);
    }
    deviation= sqrt(deviation/(R-start_value-1));
    return deviation;
}



int N;
int N_md;
double h;
double J;

int main()
{
fstream f;//this let us open/write in external Data
    f.open("convergence.dat", ios::out);
//first we start with the convergence of the leapfrog 
    double phi_0=2;
    double p_0=2;
    N=20;
    h=0.5;
    J=1;
    double array_p_phi[2];
    double difference_of_hamiltonian;
    for (int i = 1; i <= 100; i++)// we use 100 different values for N_md and calculate the difference
    {
        array_p_phi[0]=p_0;
        array_p_phi[1]=phi_0;
        leapfrog(array_p_phi, i,J,  h, N);
        difference_of_hamiltonian= (Hamiltionian(p_0,phi_0, J, N, h)-Hamiltionian(array_p_phi[0],array_p_phi[1], J, N, h))/Hamiltionian(p_0,phi_0, J, N, h);
        f<< i <<' '<< difference_of_hamiltonian <<endl;
    }
    
f.close();
// now we want to calculate the long range ising modell
    f.open("long_range_ising.dat", ios::out);
    N_md= 100;
    h=0.5;
    double array[5000];
    int success_rate;
    double magnitization_per_site[4],energy_per_site[4], error[4], error_m[4],error_e[4];
    double p_start, phi_start, random_number,Hstart,Hend, expactation_value;

    default_random_engine generator;// we need this to get an random number from a normal deistribution
    normal_distribution<double> distribution(0,1);
    
    for (int k = 0; k < 50; k++)// We choose 50 different for J between 0.2-2 
    {
        J= 0.2+1.8*k/49;
    for (int n = 0; n < 4; n++)// we choose 4 different values for N between 5-20
    {
        N=5*(n+1);
            
        success_rate=0;
    array_p_phi[1]= 100;
    
    for (int i = 0; i < 5000; i++)
    {
    random_number=((double)rand() / (double)RAND_MAX);
    array_p_phi[0]= distribution(generator);
    phi_start =array_p_phi[1];//start values stored
    p_start=array_p_phi[0];
    leapfrog(array_p_phi, N_md, J, h, N);//leapfrog algortihm
    Hstart=Hamiltionian(p_start,phi_start, J, N,  h);
    Hend=Hamiltionian(array_p_phi[0],array_p_phi[1], J, N,  h);

    if(random_number<=exp(-(Hstart-Hend)))//calculate the diffrence and choose if we want to accept it
    {
        success_rate+=1;
        array[i]=array_p_phi[1];//we did accept and write the start value in an array
    }
    else
    {
    array_p_phi[0]=p_start;
    array_p_phi[1]=phi_start;
    array[i]=array_p_phi[1];// we reject the values and write the old ones in the array
    }
    }
    magnitization_per_site[n]=0;
    energy_per_site[n]=0;
    expactation_value=0;
    for (int i = 0; i < 4000; i++)// we consider thermalization. so we take only 4000 samples
    {
        magnitization_per_site[n]+=magnitization(array[i+1000],h);//calculate the magnitization
        energy_per_site[n]+=energy(array[i+1000],h,N,J);//calculate the energy
        expactation_value+=array[i+1000];//calculate the error of phi
    }
    expactation_value=expactation_value/4000;
    magnitization_per_site[n]=magnitization_per_site[n]/4000;
    energy_per_site[n]=energy_per_site[n]/4000;
    error[n]=standard_deviation(5000,1000, expactation_value,array);//error analysis
    error_m[n]= 1/(cosh(h+expactation_value)*cosh(h+expactation_value))*error[n];
    error_e[n]= (expactation_value/J*h/(cosh(h+expactation_value)*cosh(h+expactation_value)))*error[n];
    }
        f<< J << ' '<< magnitization_per_site[0]<< ' ' <<error_m[0]<<  ' '<< magnitization_per_site[1]<< ' ' <<error_m[1]<< ' '<< magnitization_per_site[2]<< ' ' <<error_m[2]<< ' '<< magnitization_per_site[3]<< ' ' <<error_m[3]<< ' '<< energy_per_site[0]<< ' ' <<error_e[0]<< ' ' << energy_per_site[1]<< ' ' <<error_e[1]<< ' '<< energy_per_site[2]<< ' ' <<error_e[2]<< ' '<< energy_per_site[3]<< ' ' <<error_e[3]<< ' '<<magnitization_analytic(J,h,5) << ' ' <<magnitization_analytic(J,h,10) << ' '<<magnitization_analytic(J,h,15) << ' '<<magnitization_analytic(J,h,20) << ' '<< energy_analytic(J,h,5)<< ' '<< energy_analytic(J,h,10)<< ' '<< energy_analytic(J,h,15)<< ' '<< energy_analytic(J,h,20)<<endl;
        
    }
    
    

f.close();

}