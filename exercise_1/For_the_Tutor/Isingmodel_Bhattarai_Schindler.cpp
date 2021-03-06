#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>

using namespace std;

long double H(int s[], double h, int N)// Hamiltonian
{
    long double value=0;
       for (int i=0; i<N-1; i++)
{
   value+=-s[i]*s[i+1]-h*s[i];
}
    value += -s[N-1]*s[0]-h*s[N-1];
    return value;
}

long double z(double h, int N, double T)//function for the partition function
{
    int s[N];//Spin-Configuration-Array
    long double Z=0;
    for(int i=0;i<N;i++)//first we set all spins to s=-1
    {
        s[i]=-1;
    }

    for(int j=0;j<pow(2,N);j++)// We have 2^N possibilitys for our Spin configurations
    {
        for(int l=0; l<N;l++)// this is like counting in binary System. With the exception that we have -/+1, not 0/1. After increment we check if we have a carry flag 
        {     
            if(s[l]==3)
            {
                s[l]=-1;
                s[l+1]+=2;
            }
        }
        Z+=exp(-H(s,h,N)/T);
        s[0]+=2;
    }
    return Z;
}
long double z_analytic(double h, int N, double T)//analytic way to get Z
{
    long double Z= pow(exp(1/T)*(cosh(h/T)+sqrt(sinh(h/T)*sinh(h/T)+exp(-4/T))),N);
    Z+= pow(exp(1/T)*(cosh(h/T)+sqrt(sinh(h/T)*sinh(h/T)+exp(-4/T))),N);
    return Z;
}
long double z_analyticexact(double h, double T)//analytic way to get Z for N->00
{ 
    long double Z= exp(1/T)*cosh(h/T)+sqrt(exp(2/T)*sinh(h/T)*sinh(h/T)+exp(-2/T));
    return Z;
}
//variables initialize
long double magnetization;
long double magnetization_analytic;
long double magnetization_analytic_exact;
double T;
int N;
double h;
long double Z;
double error_h;
double error_N;
double error_gesamt;





int main()
{
    T=1;//We let the Temperature T const= 300 
    double delta_h = 10e-6;//this is the infintesimal for our derivativ
    fstream f;//this let us open/write in external Data

   
    //N=const
    f.open("Nconst.dat", ios::out);//Open Nconst.dat
     N=20;//first we start with n=const
    for(int m=0;m<30;m++){//We do it for 30 different h for -1 to 1
    h=-1.0+2.0*m/29;
    magnetization= T/N*((log(z(h+delta_h,N,T))-log(z(h-delta_h,N,T)))/(2*delta_h));//formular on the sheet with numerical derivativ.

    //Error analysis. Only in h dependancy because we saw that for different N the magnitization per Spin does not change
    error_h= T/N*((log(z(h*h+delta_h,N,T))-log(z(h*h-delta_h,N,T)))/(2*delta_h));
    error_h+= -(T/N*((log(z(h+delta_h,N,T))-log(z(h-delta_h,N,T)))/(2*delta_h))*T/N*((log(z(h+delta_h,N,T))-log(z(h-delta_h,N,T)))/(2*delta_h)));
    error_h=sqrt(fabs(error_h));
    

    magnetization_analytic=T/N*((log(z_analytic(h+delta_h,N,T))-log(z_analytic(h-delta_h,N,T)))/(2*delta_h));//analytic formular on the sheet
    magnetization_analytic_exact=T*((log(z_analyticexact(h+delta_h,T))-log(z_analyticexact(h-delta_h,T)))/(2*delta_h));//analytic formular on the sheet

    //Write in Nconst.dat. endl makes a new line
    f << h <<' '<< magnetization << ' '<< error_h<<' '<< magnetization_analytic<< ' '<< magnetization_analytic_exact << endl;
    }
    f.close();//we have to close Nconst.dat

    //h=const
    f.open("hconst.dat", ios::out);//open hconst.dat
    h=0.5;//h is no const
    for(int m=1;m<21;m++){// we do it vor N=1...20
    N=m;
    magnetization= T/N*((log(z(h+delta_h,N,T))-log(z(h-delta_h,N,T)))/(2*delta_h));//formular on the sheet with numerical derivativ.
    magnetization_analytic=T/N*((log(z_analytic(h+delta_h,N,T))-log(z_analytic(h-delta_h,N,T)))/(2*delta_h));//analytic formular on the sheet

    //Error analysis. Only in h dependancy because we saw that for different N the magnitization per Spin does not change
    error_h= T/N*((log(z(h*h+delta_h,N,T))-log(z(h*h-delta_h,N,T)))/(2*delta_h));
    error_h+= -(T/N*((log(z(h+delta_h,N,T))-log(z(h-delta_h,N,T)))/(2*delta_h))*T/N*((log(z(h+delta_h,N,T))-log(z(h-delta_h,N,T)))/(2*delta_h)));
    error_h=sqrt(fabs(error_h));

    //write in hconst.dat
    f << N <<' '<< magnetization << ' '<< error_h<< ' '<< magnetization_analytic << endl;
    }
    f.close();//close hconst

return 0;//if nothing goes wrong we get 0 in the end
}