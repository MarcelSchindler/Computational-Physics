#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <random>
#include <ctime>
#include <chrono>

using namespace std;


double hamilton(double u[][3], double phi[][3],double a_max, int N_max, int level)//Hamilton operator. for lower lvls we see that we just have to exchange a with 2a and N->N/2
{
    double H_a=0;
    int N=N_max/pow(2,level);
    double a=a_max*pow(2,level);
    for (int i = 0; i < N; i++)
    {
        H_a+=(u[i][level]-u[i-1][level])*(u[i][level]-u[i-1][level])/a;
    }
    for (int i = 0; i < N/2-1; i++)
    {
        H_a+=a*phi[i][level]*u[i][level];
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
    for (int i = 1; i < N-1; i++)
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


double magnitization(double u[][3], double a_max, int N_max, int level)
{
    int N=N_max/pow(2,level);
    double a=a_max*pow(2,level);
    double magnitization=0;
    for (int i = 0; i < N-1; i++)
    {
        magnitization+=1.0/N*u[i][level];
    }
    return magnitization; 
}
double magnitization_square(double u[][3], double a_max, int N_max, int level)
{
    int N=N_max/pow(2,level);
    double a=a_max*pow(2,level);
    double magnitization=0;
    for (int i = 0; i < N-1; i++)
    {
        magnitization+=1.0/N*u[i][level]*u[i][level];
    }
    return magnitization; 
}
double get_normal_random()
{
    default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
    normal_distribution<double> distribution(0,1);
    //random_device rd;
    //default_random_engine generator(rd());
    return distribution(generator);
}
void sweep(double u[][3], double phi[][3],double a_max, int N_max, int level, double delta)//here is one sweep programmed
{
    //default_random_engine generator;// we need this to get an random number from a normal deistribution
    //normal_distribution<double> distribution(0,1);
    int N=N_max/pow(2,level);//we need to adjust N and a to the given lvl
    double a=a_max*pow(2,level);

    double H_start,H_end, random_number, u_changed;
    int target_site;
    for (int i = 0; i < N-1; i++)
        {
        random_number=((double)rand() / (double)RAND_MAX);
        H_start=hamilton(u,phi,a,N,level);
        target_site=1+(rand() % (N-2));
        u_changed=u[target_site][level];
        u[target_site][level]=u[target_site][level]+get_normal_random()*delta;
        H_end=hamilton(u,phi,a,N,level);
        if(random_number>exp(-(H_end-H_start)))
        {
            u[target_site][level]=u_changed;
        }
        }
}
void phi_calculation(double phi[][3] ,double u[][3], double a ,int N_max, int level)
{
    int N=N_max/pow(2,level);
    phi[0][level+1]=(phi[N-1][level]+2*phi[0][level]+phi[1][level])/4;
    phi[0][level+1]+=a*a*(u[N-1][level]-u[N-2][level]+u[0][level]-u[2][level]-u[N-3][level]+u[N-1][level]);
    for (int i = 1; i < N/2; i++)
    {
        phi[i][level+1]=(phi[2*i-1][level]+2*phi[2*i][level]+phi[(2*i+1 %N)][level])/4;
        phi[i][level+1]+=a*a*(u[2*i-1][level]-u[2*i-2][level]+u[2*i][level]-u[(2*i+2 )%N][level]-u[2*i-3][level]+u[2*i-1][level]);
    }
}
void multigrid(double u[][3], double phi[][3], double a, int N, int level, int sweep_lvl[], int gamma, int max_level, double delta)//this is the multigrid cycle
{
    for (int i = 0; i < sweep_lvl[level]; i++)//first we sweep over the given lvl
    {
        sweep(u, phi,a,N,level,delta);
    }
    if (level< max_level)//if we are not on the coarsest lvl we do the following
    {
       
        phi_calculation(phi,u, a,N,level);//calculate new phi
        for (int i = 0; i < gamma; i++)
        {
            for (int j = 0; j < 64; j++)//put u(2a) \equiv 0
            {
                u[j][level+1]=0;
            }
        multigrid(u, phi, a, N, level+1, sweep_lvl, gamma, max_level, delta);//here is the recursice part
        }
        coarse_to_fine(u,N,level);//get the new fine values
    }
    for (int i = 0; i < sweep_lvl[level]; i++)//again sweep over
    {
        sweep(u,phi, a,N,level,delta);
    }
}
double  autocorrelation(double array[], double expected_value, int tau, int arraylenght, int startvalue)//autokorallation from last week
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
    //default_random_engine generator;// we need this to get an random number from a normal deistribution
    //normal_distribution<double> distribution(0,1);
    double u[64][3]={0};
    double phi[64][3]={0};

    fstream f;//this let us open/write in external Data
    f.open("Metropolis_Hasting.dat", ios::out);
    delta=2;
    N=64;
    L=64;
    double magnitization_per_site[15000];
    double magnitization_per_site_square[15000];
    double energy_per_site[15000];
    double expactation_value;
    double correlation_function[13000];
    
    a= L/N;

    for (int j = 0; j < 15000; j++)// we do 15000 measurements
    {
        sweep(u,phi,a,N,0,delta);
        magnitization_per_site[j]=magnitization(u,a,N,0);
        magnitization_per_site_square[j]=magnitization_square(u,a,N,0);
        energy_per_site[j]= 1.0/N*hamilton(u,phi,a,N,0);
    }

    for (int i = 0; i < 15000; i++)
    {
         f<< i << ' '<< magnitization_per_site[i]<< ' '<< magnitization_per_site_square[i]<< ' '<< energy_per_site[i]<< endl;
    }
    
    
    f.close();

    f.open("Multigrid_gamma1.dat", ios::out);
    
    int gamma;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 64; j++)
        {
            phi[j][i]=0;
            u[j][i]=0;
        }  
    }
    int sweep_lvl[3]={4,2,1};
    int max_level=2;
    gamma=1;
    expactation_value=0;
    for (int n = 0; n < 15000; n++)
    {
        multigrid(u,phi,a,N,0, sweep_lvl ,gamma,max_level, delta);// here again we do 15k measurements
        magnitization_per_site_square[n]=magnitization_square(u,a,N,0);
        magnitization_per_site[n]=magnitization(u,a,N,0);
    }
 
    for (int i = 0; i < 13000; i++)//autokorellation 
    {
        expactation_value+=magnitization_per_site_square[i+2000];
    }
    expactation_value=expactation_value/13000;
    for (int i = 0; i < 13000; i++)
    {
        correlation_function[i]=autocorrelation(magnitization_per_site_square, expactation_value, i, 15000, 0);
    }

    for (int i = 0; i < 15000; i++)
    {
         f<< i << ' '<< magnitization_per_site_square[i]<< ' '<< correlation_function[i]<< ' '<< magnitization_per_site[i] << endl;
    }
    f.close();


    // And do the same for gamma=2
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 64; j++)
        {
        phi[j][i]=0;
        u[j][i]=0;
        }
        
    }
    max_level=2;
    gamma=2;
    expactation_value=0;
    for (int n = 0; n < 15000; n++)
    {
        multigrid(u,phi,a,N,0, sweep_lvl ,gamma,max_level, delta);
        magnitization_per_site_square[n]=magnitization_square(u,a,N,0);
        magnitization_per_site[n]=magnitization(u,a,N,0);
        
    }
    for (int i = 0; i < 13000; i++)
    {
        expactation_value+=magnitization_per_site_square[i+2000];
    }
    expactation_value=expactation_value/13000;
    for (int i = 0; i < 13000; i++)
    {
        correlation_function[i]=autocorrelation(magnitization_per_site_square, expactation_value, i, 15000, 0);
    }
    f.open("Multigrid_gamma2.dat", ios::out);
    for (int i = 0; i < 15000; i++)
    {
         f<< i << ' '<< magnitization_per_site_square[i]<< ' '<< correlation_function[i] << ' '<< magnitization_per_site[i]<< endl;
    }
    f.close();

}