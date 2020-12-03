
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
double  autocorrelation(double array[], double expected_value, int tau, int arraylenght, int startvalue)
{
    int count=0;
    double C=0;
    double C_0=0;
    for (int i = 0; i < arraylenght-startvalue-tau-1; i++)
    {
          C+=(array[i+startvalue]-expected_value)*(array[i+startvalue+tau]*expected_value);
          count+=1;  
    }
    for (int i = 0; i < arraylenght-startvalue; i++)
    {
        C_0+=(array[i+startvalue]-expected_value)*(array[i+startvalue]*expected_value);
    }
    
    C= C*(arraylenght-startvalue)/(count*C_0);
    return C;
}


int N;
int N_md;
double h;
double J;
double array_p_phi[2];
double array[12800], array_spare[12800];
double correlation_function[8000];

int main()
{
// now we want to calculate the long range ising modell
fstream f;//this let us open/write in external Data
    f.open("comparison.dat", ios::out);

    h=0.5;
    N=5;
    J=0.1/N;

    int success_rate;
    double magnitization_per_site[2], error[2], error_m[2];
    double magnitization_chain[10000];
    double p_start, phi_start, random_number,Hstart,Hend, expactation_value;

    default_random_engine generator;// we need this to get an random number from a normal deistribution
    normal_distribution<double> distribution(0,1);

    for (int n = 0; n < 2; n++)
    {
      N_md=4+96*n;
    
    success_rate=0;
    array_p_phi[1]= 0.5;
    
    for (int i = 0; i < 12800; i++)
    {
    random_number=((double)rand() / (double)RAND_MAX);
    array_p_phi[0]= distribution(generator);
    phi_start =array_p_phi[1];//start values stored
    p_start=array_p_phi[0];
    leapfrog(array_p_phi, N_md, J, h, N);//leapfrog algortihm
    Hstart=Hamiltionian(p_start,phi_start, J, N,  h);
    Hend=Hamiltionian(array_p_phi[0],array_p_phi[1], J, N,  h);

    if(random_number<=exp(-(Hend-Hstart)))//calculate the diffrence and choose if we want to accept it
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



    if (n==0)
    {
    for (int i = 0; i < 12800; i++)
    {
        array_spare[i]=array[i];// so we have both arrays for N=4 and N=100
    }
    }
    cout << "sucess rate: "<< success_rate<< " for "<< (4+96*n)<<endl;
    }
    for (int i = 0; i < 8000; i++)
    {
        magnitization_per_site[0]=magnitization(array_spare[i+4800],h);
        magnitization_per_site[1]=magnitization(array[i+4800],h);
         f<< i << ' '<< array_spare[i+4800]<< ' '<< array[i+4800]<< ' '<< magnitization_per_site[0]<< ' '<< magnitization_per_site[1] << endl;
    } 

f.close();


    f.open("correlation.dat", ios::out);
    expactation_value=0;
    for (int i = 0; i < 8000; i++)
    {
        magnitization_chain[i]=magnitization(array[i+4800],h);//calculate the magnitization
        expactation_value+=magnitization(array[i+4800],h);//calculate the error of phi
    }
    expactation_value=expactation_value/8000;
    for (int i = 0; i < 8000; i++)
    {
        correlation_function[i]=autocorrelation(magnitization_chain, expactation_value, i, 8000, 0);
        f<< i<< ' ' << correlation_function[i]<< endl;
    }

    f.close();


    f.open("blocking.dat", ios::out);
    int b, blocks,arraylenght;
    arraylenght=8000;
    double correlation_blocking[4000][6];
    double standard_error[4000][6];
    double expectation_value_blocking=0;

    for (int k = 0; k < 6; k++)
    {
    b= pow(2,k+1);// so we get 2,4,8,16,32,64
    blocks= arraylenght/b;// this is how many blocks we get

    double blocking_array[blocks]={0};

    for (int i = 0; i < blocks; i++)
    {
        for (int j = 0; j < b; j++)// we create the blocks
        {
            blocking_array[i]+=magnitization(array[i*b+j],h);
        }
        blocking_array[i]=blocking_array[i]/b;
        expectation_value_blocking+=blocking_array[i];
    }
    expectation_value_blocking=expectation_value_blocking/blocks;
    
    for (int i = 1; i < blocks; i++)//calculate the autocorrelation and standard error
    {
        standard_error[i][k]= (expectation_value_blocking-blocking_array[i])*(expectation_value_blocking-blocking_array[i]);// standard deviation
        standard_error[i][k]=sqrt(standard_error[i][k]/i);
        standard_error[i][k]= standard_error[i][k]/(sqrt(blocks));// given the formula on the sheet
        correlation_blocking[i][k]=autocorrelation(blocking_array, expactation_value, i, blocks, (int)1000*1/(b*b*b));// the last term is for thermalization. so with higher b we dont throw away so many values because we have less and the thermalization is shorter.
    }
    for (int i = blocks; i < 4000; i++)
    {
        correlation_blocking[i][k]=0;
        standard_error[i][k]=0;
    }
    }
    for (int i = 1; i < 4000; i++)
    {
        f<< i<< ' ' << correlation_blocking[i][0]<< ' ' << standard_error[i][0]<<' '<< correlation_blocking[i][1]<< ' ' << standard_error[i][1]<< ' ' << correlation_blocking[i][2]<< ' ' << standard_error[i][2]<< ' ' << correlation_blocking[i][3]<< ' ' << standard_error[i][3]<< ' ' << correlation_blocking[i][4]<< ' ' << standard_error[i][4]<< ' ' << correlation_blocking[i][5]<< ' ' << standard_error[i][5]<< endl;
    }

    f.close();

    f.open("bootstrap.dat", ios::out);
    int N_bs=1500;

    double magnitization_bootstrap_average[N_bs]={0};
    double bootstrap_delta[N_bs]={0};
    double magnitization_N_bs[N_bs]={0};
    for (int j = 0; j < N_bs; j++)
    {
        for (int i = 0; i < 1000; i++)// randomly sample magnitizations
        {
          magnitization_bootstrap_average[j]+=magnitization(array[4800+rand() % 8000],h);
        }
        magnitization_bootstrap_average[j]=magnitization_bootstrap_average[j]/1000;
    }

    for (int i = 0; i < N_bs; i++)
    {
        for (int j = 0; j <= i; j++)// so we get an array with i different avarages for the magnitization
        {
            magnitization_N_bs[i]+=magnitization_bootstrap_average[j];
        }
        magnitization_N_bs[i]=magnitization_N_bs[i]/(i+1);
       
    } 
    for (int i = 1; i < N_bs; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            bootstrap_delta[i]+=(magnitization_N_bs[i]-magnitization_N_bs[j])*(magnitization_N_bs[i]-magnitization_N_bs[j]); // we use standard deviation to get the bootstrap
        }
        bootstrap_delta[i]=sqrt(bootstrap_delta[i]/(i));

       
    }
    
    for (int i = 1; i < N_bs; i++)
    {
        f<< i <<' ' <<bootstrap_delta[i]<< endl;
    }

    f.close();
    

}