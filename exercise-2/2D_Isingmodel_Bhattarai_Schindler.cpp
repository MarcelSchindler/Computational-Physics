#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>

using namespace std;


double elliptische_integral(double m)//Elliptische Intrgral numerical. I did vary the number of points and saw that it is stable
{
    double integral_value =0;
    double t;
    for (int i = 0; i < 1000; i++)
    {   
        t= (double)i/1000;
        integral_value += (double)1/1000* 1.0/sqrt((1-t*t)*(1-m *t*t));
    }
    return integral_value;
}
long double H(int s[][20], double h, int N, double J)// Hamiltonian
{
    long double value=0;
    for (int j = 0; j < N-1; j++)
    {
        
       for (int i=0; i<N-1; i++)
    {
         value+=-J*s[i][j]*s[(i+1)%N][j]-h*s[i][j];
         value+=-J*s[i][j]*s[i][(j+1)%N];

    }
    }
    return value;
}


void energiedifference(int s[][20], int N, double T, double J, double h)//calculates energiedifference(formal in the report)
{
    int random_number1 = rand() % N;
    int random_number2 = rand() % N;
    double random_number_comparison=((double)rand() / (double)RAND_MAX);
    double delta_S;
    delta_S=  2*h*s[random_number1][random_number2];
    delta_S+= 2*J*s[random_number1][random_number2]*(s[(random_number1+1)%N][random_number2]+s[random_number1][(random_number2+1)%N]);
    if (random_number1==0)
    {
        delta_S+= 2*J*s[random_number1][random_number2]*s[N-1][random_number2];
    }
    else
    {
        delta_S+= 2*J*s[random_number1][random_number2]*s[random_number1-1][random_number2];
    }
      if (random_number2==0)
    {
        delta_S+= 2*J*s[random_number1][random_number2]*s[random_number1][N-1];
    }
    else
    {
        delta_S+= 2*J*s[random_number1][random_number2]*s[random_number1][random_number2-1];
    }

    if( random_number_comparison< exp(-delta_S/T))
    {   
        s[random_number1][random_number2] *= (-1);
    }
}
double magnitization_analytic(double J)//analytic way to get magnitization depending on J
{
   if(J>0.440686793509772 )
   {
        return pow(1-(1/pow(sinh(2*J),4)),1/8);
   }
   else
   {
       return 0;
   }
}
double energie_per_site(double J)//energie per site analytically
{
   return (-J*cosh(2*J)*(1+(2/M_PI)*(2*tanh(2*J)*tanh(2*J)-1)* elliptische_integral(4*tanh(2*J)*tanh(2*J)/(cosh(2*J)*cosh(2*J)))));
}

//variables initialize
long double magnetization;
long double magnetization_analytic;
long double measurement[30000];
double  partition_function;
double energie_per_site_numerical;
double energie_per_site_analytical;

double T;
int N;
double h;
long double Z;
double error_h;
double error_N;
double error_gesamt;
int count;




int main()
{
    double measurement_average;
    int times_of_measurement=10000;// this is the number how often we want to measure
    T=1;
    double delta_h = 10e-6;//this is the infintesimal for our derivativ
    
    N=20;
    
    int s[20][20];//Spin-Configuration-Array

     fstream f;//this let us open/write in external Data
    f.open("h_variation.dat", ios::out);
    for (int k = 0; k < 30; k++)
    {
        magnetization=0;
        h=-1.0+2.0*k/29;

        for (int j = 0; j < N; j++)
    {
        for(int i=0;i<N;i++)//first we set all spins to s=-1
        {
          s[i][j]=-1;
        }
    }

    double J=0.3;
    for (int l = 0; l < times_of_measurement; l++)//for-loop for the number of measurements
    {
        for (int n = 0; n < N*N; n++)//for-loop for the sweep over one lattice
        {
          energiedifference(s, N, T,J,h);
        }
    measurement[l]= T/(N*N)*(((-H(s,h+delta_h,N,J)/T)-(-H(s,h-delta_h,N,J)/T))/(2*delta_h));// all measurements stored
    magnetization+= measurement[l];// add all measurements
    }

        f << h  << ' ' << magnetization/times_of_measurement <<endl;
    
    }
    

    
    
    f.close();

    //J variation
    h=0;
    f.open("J_variation.dat", ios::out);
    for (int k = 0; k < 60; k++)
    {
        magnetization=0;
        J=0.25+1.75*k/59;

        for (int j = 0; j < N; j++)
    {
        for(int i=0;i<N;i++)//first we set all spins to s=-1
        {
          s[i][j]=-1;
        }
    }

    
    for (int l = 0; l < times_of_measurement; l++)
    {
        for (int n = 0; n < N*N; n++)
        {
          energiedifference(s, N, T,J,h);
        }
    measurement[l]= T/(N*N)*(((-H(s,h+delta_h,N,J)/T)-(-H(s,h-delta_h,N,J)/T))/(2*delta_h));// derivative numerically
    magnetization+= measurement[l];
    }
        energie_per_site_numerical= H(s,h,N,J)/T;
        energie_per_site_numerical= energie_per_site_numerical/(N*N);// we want per site value of the energie
        energie_per_site_analytical= energie_per_site(J);
        f << J << ' '<< 1/J  << ' ' << fabs(magnetization/times_of_measurement) << ' ' << magnitization_analytic(J)<< ' ' << energie_per_site_numerical <<' ' << energie_per_site_analytical<< endl;
    
    }
    f.close();//we have to close Nconst.dat
    

return 0;//if nothing goes wrong we get 0 in the end
}