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
double elliptische_integral_second(double m)//Elliptische Intrgral numerical. I did vary the number of points and saw that it is stable
{
    double integral_value =0;
    double t;
    for (int i = 0; i < 10000; i++)
    {   
        t= (double)i/10000;
        integral_value += (double)1/10000* sqrt((1-m*t*t)/(1-t*t));
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

long double H_square(int s[][20], double h, int N, double J)// Hamiltonian square
{
    long double value=0;
    long double value_average=0;
    for (int j = 0; j < N-1; j++)
    {
        
       for (int i=0; i<N-1; i++)
    {
        
         value+=-J*s[i][j]*s[(i+1)%N][j];
         value+=-J*s[i][j]*s[i][(j+1)%N];
    if (i==0)
    {
        value+=-J*s[i][j]*s[(N-1)][j];
    }
    else
    {
        value+=-J*s[i][j]*s[(i-1)][j];
    }
      if (j==0)
    {
        value+=-J*s[i][j]*s[i][(N-1)];
    }
    else
    {
        value+=-J*s[i][j]*s[i][(j-1)];
    }
    value_average+= value*value/4;
    value=0;
    }  
    }
    return value_average;
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

double extra_homework_analytic(double J)
{
    double kappa= 2*sinh(2*J)/(cosh(2*J)*cosh(2*J));
    return 4*J*J/(M_PI*tanh(2*J)*tanh(2*J))*(elliptische_integral(kappa*kappa)-elliptische_integral_second(kappa*kappa)-(1-tanh(2*J)*tanh(2*J))*(M_PI/2*(2*tanh(2*J)*tanh(2*J)-1)*elliptische_integral(kappa*kappa)));
}

//variables initialize
long double magnetization[5];
long double magnetization_analytic;
double  partition_function;
double energie_per_site_numerical[5];
double energie_per_site_analytical;
double measurement_average;
double epsilon_square[5];

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

    double C[5];
    int times_of_measurement=2000;// this is the number how often we want to measure
    int times_of_consider_measurment=times_of_measurement-1000;// we dont want the first values because they have a higher error
    T=1;
     double J=0.3;
    double delta_h = 10e-6;//this is the infintesimal for our derivativ
    
    int s[20][20];//Spin-Configuration-Array

     fstream f;//this let us open/write in external Data
    f.open("h_variation.dat", ios::out);

    
    
    for (int k = 0; k < 100; k++)//times of measurement
    {
        
        h=-1.0+2.0*k/99;
    for (int o = 0; o < 5; o++)//for different N
    {
        magnetization[o]=0;
        N=4+o*4;
        for (int j = 0; j < N; j++)
    {
        for(int i=0;i<N;i++)//first we set all spins to s=-1. cold start
        {
          s[i][j]=-1;
        }
    }

   
    for (int l = 0; l < times_of_measurement; l++)//for-loop for the number of sweeps
    {
        for (int n = 0; n < N*N; n++)//for-loop for the sweep over one lattice
        {
          energiedifference(s, N, T,J,h);
        }
    if(l>times_of_measurement-times_of_consider_measurment)
        {
            magnetization[o]+= (((-H(s,h+delta_h,N,J)/T)-(-H(s,h-delta_h,N,J)/T))/(2*delta_h))/(N*N);
        }
    }

    
    }
        f << h  << ' ' << magnetization[0]/times_of_consider_measurment << ' '<< magnetization[1]/times_of_consider_measurment <<' '<<magnetization[2]/times_of_consider_measurment<<' '<< magnetization[3]/times_of_consider_measurment<<' '<<magnetization[4]/times_of_consider_measurment<<endl;
    }

    

    
    
    f.close();

    //J variation
    h=0;
    f.open("J_variation.dat", ios::out);

    
    for (int k = 0; k < 500; k++)//times of measurement
    {
        
        J=0.25+1.75*k/499;
        for (int o = 0; o < 5; o++)//for different N
    {
        magnetization[o]=0;
        N=4+4*o;
        for (int j = 0; j < N; j++)
    {
        for(int i=0;i<N;i++)//first we set all spins to s=-1. cold start
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
        if(l>times_of_measurement-times_of_consider_measurment)
        {
    magnetization[o]+= T/(N*N)*(((-H(s,h+delta_h,N,J)/T)-(-H(s,h-delta_h,N,J)/T))/(2*delta_h));// derivative numerically
        }
    }
        energie_per_site_numerical[o]= H(s,h,N,J)/(N*N);//energy per site
        epsilon_square[o]=H_square(s,h,N,J)/(N*N);
        //cout<< epsilon_square<<' '<<  energie_per_site_numerical[o]<< endl;
        energie_per_site_analytical= energie_per_site(J);//thermodynamic limit
        C[o]=(N*N)*(epsilon_square[o]-energie_per_site_numerical[o]*energie_per_site_numerical[o]);
    }
    //cout<< epsilon_square[4]<<' '<<  energie_per_site_numerical[4]<< endl;
        
       f << J << ' '<< 1/J  << ' ' << fabs(magnetization[0]/times_of_consider_measurment)<<' ' <<fabs(magnetization[1]/times_of_consider_measurment)<<' '<<fabs(magnetization[2]/times_of_consider_measurment)<<' '<<fabs(magnetization[3]/times_of_consider_measurment)<<' '<< fabs(magnetization[4]/times_of_consider_measurment) << ' ' << magnitization_analytic(J)<< ' ' << energie_per_site_numerical[0] <<' ' << energie_per_site_numerical[1]<<' ' << energie_per_site_numerical[2]<<' ' << energie_per_site_numerical[3]<<' ' << energie_per_site_numerical[4]<<' ' << energie_per_site_analytical<<' ' << (magnetization[4]/times_of_consider_measurment)<<  ' ' <<C[0]/(J*J)<<' ' <<C[1]/(J*J)<<' ' <<C[2]/(J*J)<<' ' <<C[3]/(J*J)<< ' ' <<C[4]/(J*J)<<' '<< N*N*extra_homework_analytic(J)/(J*J)<< endl;
    }
 
    
    f.close();//we have to close Nconst.dat
    

return 0;//if nothing goes wrong we get 0 in the end
}