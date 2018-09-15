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
#define Num_file        1
#define Num_file_total  5


using namespace std;

int main()
{
    int cache_random(int F,double gamma);
    int transmit_random(int F,double gamma);
    int priority_random(double mi[],int M,int k);
    int possion(double Lambda);
    double exponential(double Lambda);
    double U_Random();

    srand( (unsigned)time( NULL ) );
    mt19937 eng;
    uniform_real_distribution<double>uniform(0.0,1.0);
    exponential_distribution<double>Exponential(1.0);
	uniform_real_distribution<double>angle(0.0,2*pi);

    int Number_t,Number_r;
    int i,j,k,l,n,s;
    int flag;
    //int tfile[1000];
	//int cfile[1000][10];
    double cf = (double)Num_file / Num_file_total;
    double lamda_T,d = 50.0,R = 15.0;
    double gamma = 1.0;
    int P_Transmission_flag,P_Transmission_flag_total;
    int seed_1 = 10000,seed_2 = 10;
    double ratio = (double)seed_1 / (double)seed_2;
    int Position_Coverage[seed_1];
	int Fading_Coverage[seed_2];
    int Random_tfile;
	int Random_cfile;
    int count_1;
    double result;
    int A_N_t = pi * pow(d,2) * lamda;
    int A_N_r = pi * pow(d,2) * lamda;
    int priority;
	double m[1000];
	double Distance;
    //printf("Number_t = %d\n",Number_t);

    double N_t = 1;
    double d_1 = 2;
    double P_1 = 5;
    double P_0 = 5;
    double alpha = 4;
    double theta = 3;

    double Random_Angle;
    double Fading,Distance_Probability;

    lamda_T = lamda;
    printf("lamda_T = %f\n",lamda_T);

    Random_Angle = angle(eng);
    double T0_X = d_1 * cos(Random_Angle);
	double T0_Y = d_1 * sin(Random_Angle);
	int PositionCoverage = 0;
	int PositionCount = 0;
	double P_T_flag_probability;
	double P_Coverage;
    double T_flag_total_average_total = 0.0;
    double T_flag_total_average;
    for(i = 0; i < seed_1; i++)
    {
        Number_t = possion(A_N_t);
        if(Number_t == 0 || Number_t == 1)
            Number_t = possion(A_N_t);
        Number_r = Number_t;

        int priority_T[Number_t];

        double TX_Position[Number_t],RX_Position[Number_r];
        int Transmission_flag[Number_t];
        double TX_Position_X[Number_t],TX_Position_Y[Number_t],RX_Position_X[Number_r],RX_Position_Y[Number_r];
        int TX_Cache[Number_t],RX_Cache[Number_r],TX_Transmit[Number_t];

        //Number_t = possion(A_N_t);
        //printf("Number_t = %d\n",Number_t);
		for(j = 0; j < Number_t; j++)
        {
            m[j] = U_Random();
            //printf("m = %f\n",m[j]);
        }
        for(j = 0; j < Number_r; j++)
        {
            priority = priority_random(m, Number_t, j);
           // printf("The priority of m%d is %d.\n",j,priority);
            priority_T[j] = priority;
            priority_T[0] = 10;
            //printf("priority[%d] = %d\n",j,priority_T[j]);


			Distance_Probability = uniform(eng);
			RX_Position[j] = R * sqrt(Distance_Probability);
			RX_Position[0] = 0.0;
			RX_Position_X[j] = RX_Position[j] * cos( Random_Angle);
			RX_Position_Y[j] = RX_Position[j] * sin( Random_Angle);
			//printf("RX_Position[%d] = %f\n",j,RX_Position[j]);

			TX_Position_X[j] = RX_Position_X[j] + d_1 * cos( Random_Angle);
			TX_Position_Y[j] = RX_Position_Y[j] + d_1 * sin( Random_Angle);
			TX_Position[j] = sqrt(pow(TX_Position_X[j],2) + pow(TX_Position_Y[j],2));

            Random_cfile = transmit_random(Num_file_total,gamma);
            TX_Cache[j] = Random_cfile;
            TX_Transmit[j] = TX_Cache[j];
            Random_cfile = transmit_random(Num_file_total,gamma);
            RX_Cache[j] = Random_cfile;

			Transmission_flag[j] = 1;
            //printf("Transmission_flag[%d] = %d\n",j,Transmission_flag[j]);

        }
        int FadingCoverage = 0;
        int T_flag_total ;

		double Interference_T0,Interference_R0,Interference_T,Interference;
		double Signal,SIR;
		int T_flag[Number_t],T_flag_X,T_flag_X_total;
		int counter = 0;


		for(k = 0;k < seed_2;k++)
        {
            double T_I = 0;

			Random_cfile = transmit_random(Num_file_total,gamma);///R0存的文件
			//RX_Cache[0] = Random_cfile;
			Random_tfile = transmit_random(Num_file_total,gamma);///T0传的文件

            for (j = 1;j<Number_r;j++)///T0周围的R是否会导致T0不能传
            {
                Fading = Exponential(eng);
                //printf("fading = %f\n",Fading);

                Distance = sqrt(pow(abs(RX_Position_X[j] - TX_Position_X[0]), 2) + pow(abs(RX_Position_Y[j] - TX_Position_Y[0]), 2));
				Interference_T0 = P_1 * Fading / pow(Distance, alpha);

				if (priority_T[j] < priority_T[0] && Transmission_flag[j])
                {
                    ///T0周围的高优先级的R对T0的干扰大于门限或R不缓存T0传的文件
					if (Interference_T0 > N_t && Random_tfile != RX_Cache[j])
					{
						Fading_Coverage[k] = 0;

						// FadingCoverage += Fading_Coverage[k];

						goto Stop;
					}
				}
            }
            counter++;
            T_flag_total = 0;

            int a[Number_t],b[Number_t],b_flag[Number_t],b_flag_total;
            n = 0;
            //priority_T[0] = 1;
            for(l = 0;l < Number_t;l++)
            {
                for(j = 0;j < Number_t;j++)
                {
                    if(priority_T[j] == l + 1)
                        break;
                }
                a[n] = j;
                n++;
            }

            for(n = 0;n < Number_t;n++)
            {
                if(n == 0)
                {
                    b[n] = 1;
                }
                else
                {
                    if(a[n] == 0)
                        b[n] = 1;
                    else
                    {
                        b_flag_total = 0;
                        for(s = 0;s < n;s++)
                        {
                            if(s != 0)
                            {
                                if(a[s] != a[s-1])
                                    Fading = Exponential(eng);///每次生成不同的fading！！！
                            }
                            else
                                Fading = Exponential(eng);
                            //printf("Fading = %f\n",Fading);
                            Distance = sqrt(pow(abs(TX_Position_X[a[n]] - RX_Position_X[a[s]]), 2) + pow(abs(TX_Position_Y[a[n]] - RX_Position_Y[a[s]]), 2));
                            //printf("Distance = %f\n",Distance);
                            Interference_T = P_1 * Fading / pow(Distance, alpha);
                            //printf("Interference_T = %f\n",Interference_T);

                            //if(TX_Transmit[a[n]] == RX_Cache[a[s]] )
                            if(TX_Transmit[a[n]] == RX_Cache[a[s]] || Interference_T <= N_t)
                            {
                                b_flag[s] = 1;
                            }
                            else
                            {
                                b_flag[s] = 0;
                                //b[n+1] = 0;
                                //a[n+1] = a[n];
                            }
                            b_flag_total += b_flag[s];

                        }
                        b[n] = b_flag_total / n;
                    }
                }

                if(b[n] == 0)
                    a[n] = a[n-1];
                //printf("b[%d] = %d\n",n,b[n]);
                //printf("a[%d] = %d\n",n,a[n]);
                T_flag[a[n]] = b[n];
                T_flag_total += T_flag[a[n]];
            }

			Interference = 0;
			for (j = 1;j < Number_t;j++)
			{
			    ///在active的T中，看哪些传的不是R0的缓存，???但对R0的干扰小于N_t???
				if (T_flag[j] == 1 && Transmission_flag[j] == 1 && TX_Transmit[j] != RX_Cache[0])//&&priority_T[j]<priority_T[0])
				{
					Fading = Exponential(eng);
					//printf("Fading = %f\n",Fading);
					Interference_R0 = P_1 * Fading / pow(TX_Position[j], alpha);
                    //printf("TX_Position[%d] = %f\n",j,TX_Position[j]);
				}
				else
                    Interference_R0 = 0;
                if(Interference_R0 <= N_t)
                    Interference += Interference_R0;
                //printf("Interference = %f\n",Interference);
			}
			Fading = Exponential(eng);
			Signal = P_0 * Fading / pow(d_1, alpha);
			//printf("Signal = %f\n",Signal);
			SIR = Signal / Interference;
			//printf("SIR = %f\n",SIR);

			if (SIR >= theta)
				Fading_Coverage[k] = 1;
			else
				Fading_Coverage[k] = 0;


			FadingCoverage += Fading_Coverage[k];
			Stop:

			cout << fixed << setprecision(6) << 100 * (i*seed_1 + k * ratio) / (seed_1* seed_2 * ratio) << "% has completed" << "\r";
        }
        T_flag_total_average = (double)(T_flag_total ) / (Number_t );
        //printf("T_flag_total_average = %f\n",T_flag_total_average);
        T_flag_total_average_total += T_flag_total_average;
        PositionCount += counter;
		PositionCoverage += FadingCoverage;

    }

    cout << endl;
    cout << "PositionCount = " << PositionCount<<endl;
    cout << "PositionCoverage = " << PositionCoverage<<endl;
    //printf("T_flag_total_average_total = %f\n",T_flag_total_average_total);
    P_T_flag_probability = T_flag_total_average_total / (double)seed_1;
    //printf("P_T_flag_probability = %f\n",P_T_flag_probability);
    P_Coverage = (double)PositionCoverage / PositionCount;
    cout << "The Coverage probability when Nt = " << fixed << setprecision(6) <<  N_t << " is " <<  P_Coverage << endl;

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
    double Num = U_Random() * Num_file_total;
    int Num_int = (int)Num;
    return (Num_int + 1);
}

/*
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
*/
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

double exponential(double Lambda)
{
    double pV = 0.0;
    while(true)
    {
        pV = (double)rand()/(double)RAND_MAX;
        if (pV != 1)
        {
            break;
        }
    }
    pV = (-1.0/Lambda)*log(1-pV);
    return pV;
}
