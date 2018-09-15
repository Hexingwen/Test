#include <random>
#include <iostream>
#include <string>
#include <vector>
#define _USE_MATH_DEFINES
#include "math.h"
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <fstream>
#define lamda           0.005
#define pi              3.14159265358979323846
#define Num_file_total  5
#define Num_file        1
#define Number          200


using namespace std;

int main()
{
    int cache_random(int F,double gamma);
    int transmit_random(int F,double gamma);
    int priority_random(double mi[],int M,int k);
    int possion(double Lambda);
    double U_Random();

    srand( (unsigned)time( NULL ) );
    int i,j,k,n,s,l,m;
    long double frac_n = 1;
    long double frac[Number];
    double cff = 1/( double )Num_file_total ;
    double pf[Num_file_total];
    double cf[Num_file_total];
    double d = 50.0;
    double t,C_d;
    double Pa,Q = 0.972653;
    long double Pa_T[1000];
    double gamma = 1.0;
    double sum = 0.0;
    double sum_f = 0.0;

    double N_t = 1;
    double d_1 = 2;
    double P_0 = 5 ;
    double P_1 = 5 ;
    double alpha = 4;

    double I_1,I_2;
    double theta = 3;
    double B = 3.06E-06;
    //double B = 3.2641*pow(10,-45);
    //double B = 5.2950*pow(10,-5);
    double Coverage_Probability = 0.0;
    double Pa_T_sum;

    for(i=0;i<Num_file_total;i++)
    {
        t = 1 / pow(i+1,gamma);
        sum_f = sum_f + t;
    }

    for(i=0;i<Num_file_total;i++)
    {
        pf[i] = (1 / pow(i+1,gamma)) / sum_f;
        cf[i] =  pf[i];
        //cf[i] = 1 / (double)Num_file_total;
        //printf("pf_%d = %f\n",(i+1),pf[i]);
    }
    I_1 = -(1 / alpha) * pow((P_0/N_t),(2/alpha))*1.77245;
    I_2 = (pow(theta,2/alpha) * pow(d_1,2) * pow(P_1/P_0,2/alpha)*pi)/(alpha * sin(2*pi/alpha));
    printf("I_1 = %f\n",I_1);
    printf("I_2 = %f\n",I_2);
    for(n = 0;n < Number;n++)
    {
        if(n==0)
            frac[0] = 1;
        else
            frac[n] = n * frac[n-1];
        //printf("frac[%d] = %lf\n",n,frac[n]);
    }
    for(l = 0;l < Num_file_total;l++)
    {
        double c_p_l[l];
        c_p_l[l] = 0;
        for(n = 1; n < Number; n++)
        {
            double Pa_pr[10],Pa_pr_sum[n],P_sum[n];
            double c_p_n[n];
            int frac_s = 1;
            c_p_n[n] = 0;
            //Pa_pr[1] = 1;
            //printf("n = %d\n",n);
            frac_n = n*frac_n;
            //P_pr_j_n = powl(lamda*pi*pow(d,2),n) / frac * exp(-lamda*pi*pow(d,2)) / n;
            for(s = 1;s <= n; s++)
            {
                double A_t,c_p_s[s];
                double sum_pr = 0.0;
                Pa_T_sum = 0.0;
                //A[s] = 0;
                //A[2] = 1;

                for(j = 1;j <= n;j++)
                {
                    if(j == 1)
                        P_sum[j] = 1.0;
                    else
                    {
                        //Pa_pr_sum[1] = 1.0;
                        Pa_pr_sum[j] = 0.0;
                        //P_sum[1] = 1.0;
                        for(k = 0;k < Num_file_total;k++)
                        {
                        //Pa_pr[k] = pf[k]*( 1 - Pa_pr_sum[j-1] * ((1 -  cf[k]))) * Pa_pr_sum[j-1];
                            //Pa_pr[k] = pf[k]*( 1 - P_sum[j-1] * ((1 -  cf[k])*(1-Q))) * P_sum[j-1];
                            Pa_pr[k] = pf[k]*( 1 - P_sum[j-1] * ((1 -  cf[k])*(1-Q))) * P_sum[j-1];
                            //P_sum[j] = ( 1 - P_sum[j-1] * ((1 -  cff)*(1-Q))) * P_sum[j-1];
                            //printf("Pa_pr[%d] = %f\n",k,Pa_pr[k]);
                            Pa_pr_sum[j] += Pa_pr[k];
                            //printf("Pa_pr_sum[%d] = %f\n",j,Pa_pr_sum[j]);
                        }
                        P_sum[j] = Pa_pr_sum[j];
                        //printf("P_sum[%d] = %f\n",j,P_sum[j]);
                        //
                    }
                    sum_pr += P_sum[j];
                    //printf("sum_pr[%d] = %f\n",j,sum_pr);
                }
                //printf("sum_pr[%d] = %f\n",s,sum_pr);

                for(m = s;m < Number;m++)
                {
                    double Pa_T[m];
                    A_t = sum_pr * (1 - cf[l]) / m;
                    //printf("A_t = %f\n",A_t);
                    Pa_T[m] =  powl(lamda*pi*pow(d,2),m) / frac[m] * exp(-lamda*pi*pow(d,2)) * A_t;
                    Pa_T_sum += Pa_T[m];
                }
                c_p_s[s] = exp(-2*pi*lamda*Pa_T_sum*I_1) * exp(-2*pi*lamda*Pa_T_sum*I_2) * exp(-2*pi*lamda*Pa_T_sum*B);
                //printf("c_p_s[%d] = %f\n",s,c_p_s[s]);
                c_p_n[n] += (c_p_s[s] / n);
                cout << fixed << setprecision(6) << 100 * (l + ( n + (s/n)) / (double)Number ) / Num_file_total << "% has completed" << "\r";

            }

            c_p_l[l] += powl(lamda*pi*pow(d,2),n) / frac[n] * exp(-lamda*pi*pow(d,2)) * c_p_n[n];

        }

        Coverage_Probability += pf[l]*c_p_l[l];
    }

    printf("The Coverage_Probability is %f\n",Coverage_Probability);

    return 0;
}

double U_Random()   /* 产生一个0~1之间的随机数 */
{
    double f;
    f = (float)(rand() % 10000);
    /* printf("%fn",f); */
    return f/10000;
}

//此函数用于产生随机文件
int transmit_random(int F,double gamma)
{
    //mt19937 eng;
    double pf[F],spf[F];
    double Num = U_Random(),sum=0;
    //uniform_real_distribution<double>uniform(0.0,1.0);
     for(int i=0;i<F;i++)
    {
        sum = sum + (1 / pow((i + 1),gamma));
    }
    for(int i=0;i<F;i++)
    {
        pf[i]=(1 / pow((i + 1),gamma))/sum;
    }
    spf[0]=pf[0];
    for(int i=1;i<F;i++)
    {
        spf[i]=spf[i-1]+pf[i];
    }
    //Num = uniform(eng);
    //printf("num1 = %f\n",Num);
    for(int i=0;i<F;i++)
    {
        if(Num<=spf[i])
            {
                return i+1;
            }
    }

}

int cache_random(int F,double gamma)
{
    //mt19937 eng;
    int i;
    double pf[F];
    //int cf[F];
    double Num = U_Random(),sum=0;
    //uniform_real_distribution<double>uniform(0.0,1.0);

    for(int i=0;i<F;i++)
    {
        sum = sum + (1 / pow((i + 1),gamma));
    }
    for(int i=0;i<F;i++)
    {
        pf[i]=(1 / pow((i + 1),gamma))/sum;
    }
    //Num = uniform(eng);
    //printf("num2 = %f\n",Num);
    for(int i=0;i<F;i++)
    {
        if(Num < pf[i])
        {
            return (i + 1);
        }
        else
            return (0);
    }
}

int priority_random(double org[],int M,int k)
{
    //mt19937 eng;
    int i,j;
    int pri[M];
    double *p_mi,*p_org;
    double t,mi[M];
    p_mi=mi;

    //uniform_real_distribution<double>uniform(0.0,1.0);
    for(i = 0; i < M; i++)
    {
        //mi[i] = U_Random();

        mi[i] = org[i];

        for(j = 0; j < i; j++)
        {
            if(mi[i] > mi[j])
            {
                t = mi[i];
                mi[i] = mi[j];
                mi[j] = t;
            }
        }

    }
    /*for(i = 0; i < M; i++)
    {
        printf("m%d = %f\n",i,mi[i]);
    }*/
    for(i = 0; i < M; i++)
    {
        for(j = 0; j < M; j++)
        {
            if(mi[j] == org[i])
                pri[i] = j + 1;
        }
    }
    return pri[k];
}

int possion(double Lambda)  /* 产生一个泊松分布的随机数，Lamda为总体平均数*/
{
        int k = 0;
        long double p = 1.0;
        long double l = exp(-Lambda);  /* 为了精度，才定义为long double的，exp(-Lambda)是接近0的小数*/
        //printf("%.15Lf\n",l);
        while (p>l)
        {
            double u = U_Random();
            p *= u;
            k++;
        }
        return (k-1);
}
