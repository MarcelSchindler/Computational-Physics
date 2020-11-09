#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>

using namespace std;

double H(int s[], double h, int N)// Hamiltonian
{
    double value=0;
       for (int i=0; i<N-1; i++)
{
   value+=-s[i]*s[i+1]-h*s[i];
}
    value += -s[N-1]*s[0]-h*s[N-1];
    return value;
}

double z(double h, int N, double T)//function for the partition function
{
    int s[N];//Spin-Configuration-Array
    double Z=0;
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
double z_analytic(double h, int N, double T)//analytic way to get Z
{
    double Z= pow(exp(1/T)*(cosh(h/T)+sqrt(sinh(h/T)*sinh(h/T)+exp(-4/T))),N);
    Z+= pow(exp(1/T)*(cosh(h/T)+sqrt(sinh(h/T)*sinh(h/T)+exp(-4/T))),N);
    return Z;
}
//variables initialize
double magnetization;
double magnetization_analytic;
double T;
int N;
double h;
double Z;





int main()
{
    T=300.0;//We let the Temperature T const= 300 
    double delta_h = 10e-6;//this is the infintesimal for our derivativ
    fstream f;//this let us open/write in external Data

    //N=const
    f.open("Nconst.dat", ios::out);//Open Nconst.dat
     N=20;//first we start with n=const
    for(int m=0;m<30;m++){//We do it for 30 different h for -1 to 1
    h=-1.0+2.0*m/29;
    magnetization= 1000*T/N*((log(z(h+delta_h,N,T))-log(z(h-delta_h,N,T)))/(2*delta_h));//formular on the sheet with numerical derivativ. Added *1000 so the plots look nice
    magnetization_analytic=1000*T/N*((log(z_analytic(h+delta_h,N,T))-log(z_analytic(h-delta_h,N,T)))/(2*delta_h));//analytic formular on the sheet
    f << h <<' '<< magnetization <<' '<< magnetization_analytic << endl;//Write in Nconst.dat. endl makes a new line
    }
    f.close();//we have to close Nconst.dat

    //h=const
    f.open("hconst.dat", ios::out);//open hconst.dat
    h=0.5;//h is no const
    for(int m=1;m<21;m++){// we do it vor N=1...20
    N=m;
    magnetization= 1000*T/N*((log(z(h+delta_h,N,T))-log(z(h-delta_h,N,T)))/(2*delta_h));//formular on the sheet with numerical derivativ.
    magnetization_analytic=1000*T/N*((log(z_analytic(h+delta_h,N,T))-log(z_analytic(h-delta_h,N,T)))/(2*delta_h));//analytic formular on the sheet
    f << N <<' '<< magnetization <<' '<< magnetization_analytic << endl;//write in hconst.dat
    }
    f.close();//close hconst

return 0;//if nothing goes wrong we get 0 in the end
}