#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <random>

using namespace std;


double hamilton(double u[][3],double a, int N_max, int level)
{
    double H_a=0;
    int N=N_max/pow(2,level);
    for (int i = 0; i < N; i++)
    {
        H_a+=(u[i][level]-u[i-1][level])*(u[i][level]-u[i-1][level])/a;
    }
    return H_a;
}

void fine_to_coarse(double u[][3], int N_max, int level)
{
    int N=N_max/pow(2,level);
    for (int i = 0; i < N/2; i++)
    {
        u[i][level+1]=u[2*i][level];
    }
    
}
void coarse_to_fine(double u[][3], int N_max, int level)
{
    int N=N_max/pow(2,level);
    for (int i = 0; i < N; i++)
    {
        if (i %2==0)
        {
            u[i][level]=u[i][level]+u[i/2][level+1];
        }
        else
        {
            u[i][level]=u[i][level]+(u[(i-1)/2][level+1]+u[(i+1)/2][level+1])/2;
        }
    }
    
}

double hamilton_coarse(double u[][3], double phi[][3], double a, int N_max, int level)
{
    int N=N_max/pow(2,level);
    double H_a=0;
    for (int i = 0; i < N/2; i++)
    {
        H_a+=(u[i][level+1]-u[i-1][level+1])*(u[i][level+1]-u[i-1][level+1])/a;
    }
    for (int i = 0; i < N/2-1; i++)
    {
        H_a+=2*a*phi[i][level+1]*u[i][level+1];
    }
    return H_a;
}

double magnitization(double u[][3], double a, int N_max, int level)
{
    int N=N_max/pow(2,level);
    double magnitization=0;
    for (int i = 0; i < N-1; i++)
    {
        magnitization+=1.0/N*u[i][level];
    }
    return magnitization; 
}
double magnitization_square(double u[][3], double a, int N_max, int level)
{
    int N=N_max/pow(2,level);
    double magnitization=0;
    for (int i = 0; i < N-1; i++)
    {
        magnitization+=1.0/N*u[i][level]*u[i][level];
    }
    return magnitization; 
}
void sweep(double u[][3], double a_max, int N_max, int level, double delta)
{
    default_random_engine generator;// we need this to get an random number from a normal deistribution
    normal_distribution<double> distribution(0,1);
    int N=N_max/pow(2,level);
    double a=a_max*pow(2,level);
    
    double random_number,H_start,H_end;
    int target_site;
    double u_changed;
    for(int i = 0; i < N-1; i++)
    {
        random_number=((double)rand() / (double)RAND_MAX);
        H_start=hamilton(u,a,N,level);
        target_site=1+rand() % (N-2);
        u_changed=u[target_site][level];
        u[target_site][level]=u[target_site][level]+distribution(generator)*delta;
        H_end=hamilton(u,a,N,level);
        if(random_number>exp(-(H_end-H_start)))
        {
            u[target_site][level]=u_changed;
        }
    }
}
void phi_calculation(double phi[][3], int N_max, int level)
{
    int N=N_max/pow(2,level);
    for (int i = 0; i < N/2; i++)
    {
        phi[i][level+1]=(phi[2*i-1][level]+2*phi[2*i][level]+phi[2*i+1][level])/4;
    }
    
}
void multigrid(double u[][3], double phi[][3], double a, int N_max, int level, int sweep_lvl[], int gamma, int max_level, double delta)
{
    int N=N_max/pow(2,level);
    if (level=! max_level)
    {
    for (int i = 0; i < sweep_lvl[level]; i++)
    {
        sweep(u,a,N,level,delta);
    }
    phi_calculation(phi,N,level);
    fine_to_coarse(u,N,level);
    for (int i = 0; i < gamma; i++)
    {
    u[64][level+1]={0};
    multigrid(u, phi, a, N, level-1, sweep_lvl, gamma, max_level, delta);
    }
    coarse_to_fine(u,N,level);
    }
    for (int i = 0; i < sweep_lvl[level]; i++)
    {
        sweep(u,a,N,level,delta);
    }
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

double delta,a,L;
double H_start,H_end, random_number, u_changed;
int target_site,N;

int main()
{
    default_random_engine generator;// we need this to get an random number from a normal deistribution
    normal_distribution<double> distribution(0,1);
    double u[64][3]={0};

    fstream f;//this let us open/write in external Data
    f.open("Metropolis_Hasting.dat", ios::out);
    delta=2;
    N=64;
    L=64;
    double magnitization_per_site[10000];
    double magnitization_per_site_square[10000];
    double energy_per_site[10000];
    double expactation_value;
    double correlation_function[8000];
    a= L/N;

    for (int j = 0; j < 10000; j++)
    {
    for (int i = 0; i < N-1; i++)
    {
        random_number=((double)rand() / (double)RAND_MAX);
        H_start=hamilton(u,a,N,0);
        target_site=1+rand() % (N-2);
        u_changed=u[target_site][0];
        u[target_site][0]=u[target_site][0]+distribution(generator)*delta;
        H_end=hamilton(u,a,N,0);
        if(random_number>exp(-(H_end-H_start)))
        {
            u[target_site][0]=u_changed;
        }
    }
        magnitization_per_site[j]=magnitization(u,a,N,0);
        magnitization_per_site_square[j]=magnitization_square(u,a,N,0);
        energy_per_site[j]= 1.0/N*hamilton(u,a,N,0);
    }

    for (int i = 0; i < 10000; i++)
    {
         f<< i << ' '<< magnitization_per_site[i]<< ' '<< magnitization_per_site_square[i]<< ' '<< energy_per_site[i]<< endl;
    }
    
    
    f.close();

    f.open("Multigrid_gamma1.dat", ios::out);
    
    int gamma;
    double phi[64][3]={0};
    u[64][3]={0};
    int sweep_lvl[3]={4,2,1};
    int max_level=2;
    gamma=1;
    expactation_value=0;
    for (int n = 0; n < 10000; n++)
    {
        multigrid(u,phi,a,N,0, sweep_lvl ,gamma,max_level, delta);
        magnitization_per_site_square[n]=magnitization_square(u,a,N,0);
    }
 
    for (int i = 0; i < 8000; i++)
    {
        expactation_value+=magnitization_per_site_square[i+2000];//calculate the error of phi
    }
    expactation_value=expactation_value/8000;
    for (int i = 0; i < 8000; i++)
    {
        correlation_function[i]=autocorrelation(magnitization_per_site_square, expactation_value, i, 8000, 0);
    }

    for (int i = 0; i < 10000; i++)
    {
         f<< i << ' '<< magnitization_per_site_square[i]<< ' '<< correlation_function[i] << endl;
    }
    f.close();


    f.open("Multigrid_gamma2.dat", ios::out);
    phi[64][3]={0};
    u[64][3]={0};
    max_level=2;
    gamma=2;
    expactation_value=0;
    for (int n = 0; n < 10000; n++)
    {
        multigrid(u,phi,a,N,0, sweep_lvl ,gamma,max_level, delta);
        magnitization_per_site_square[n]=magnitization_square(u,a,N,0);
    }
    for (int i = 0; i < 8000; i++)
    {
        expactation_value+=magnitization_per_site_square[i+2000];//calculate the error of phi
    }
    expactation_value=expactation_value/8000;
    for (int i = 0; i < 8000; i++)
    {
        correlation_function[i]=autocorrelation(magnitization_per_site_square, expactation_value, i, 8000, 0);
    }

    for (int i = 0; i < 10000; i++)
    {
         f<< i << ' '<< magnitization_per_site_square[i]<< ' '<< correlation_function[i] << endl;
    }
    f.close();

}