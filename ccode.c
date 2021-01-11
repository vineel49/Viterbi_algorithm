/* Decoding Convolutional Coded QPSK Over AWGN */
/* Generator matrix is G(D)=[1 (1+D^2)/(1+D+D^2)] */
/* overall rate is 1 */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <time.h>

void main()
{
clock_t start, end;
double cpu_time_used;
start = clock();
/*------------------------------------------------------*/    
    /* simulation parameters (overall rate is 1) */ 
float SNR_dB=6; /*SNR per bit in dB */  
int decoding_delay=20; /*decoding delay of the VA */  
int frame_size=pow(10,3); /*frame size */
double noise_var_1D;/* 1D noise variance */ 
long int sim_runs=pow(10,2); // simulation runs
noise_var_1D = (double)2/(2*pow(10,0.1*SNR_dB));
/*------------------------------------------------------*/
/*       QPSK Constellation                             */
int constell_size = 4; // QPSK constellation size
double re_qpsk_sym[4]={1,1,-1,-1},im_qpsk_sym[4]={1,-1,1,-1};
/*------------------------------------------------------*/
/* Trellis for the convolutional encoder               */
int P_State[4][2]={{0,1},{3,2},{1,0},{2,3}}; // Previous states
int P_Ip[4][2]={{0,1},{0,1},{0,1},{0,1}}; // previous inputs
int Ga_Inx[4][2]={{0,3},{1,2},{0,3},{1,2}}; // branch indices for path-metric recursion 
int num_states=4; // number of states in the trellis
/*------------------------------------------------------*/
int fr_cnt,C_BER=0;
for(fr_cnt=0;fr_cnt<sim_runs;fr_cnt++)
{
srand(time(0));
/*source*/
int a[frame_size],i2;
for(i2=0;i2<frame_size;i2++)
a[i2] = rand()%2;

// Convolutional coding (RSC)
int b[2*frame_size],D1=0,D2=0,temp;
for(i2=0;i2<frame_size;i2++)
{
b[2*i2]=a[i2]; //systematic bit
b[2*i2+1]=a[i2]^D1; // parity bit
temp=b[2*i2+1]^D2;
D2=D1;
D1=temp;
}

// QPSK mapping
int re_sym[frame_size],im_sym[frame_size];
for(i2=0;i2<frame_size;i2++)
{
    re_sym[i2] = 1-2*b[2*i2];
    im_sym[i2] = 1-2*b[2*i2+1];
}
/*-------------------------------------------------------------------*/
/* awgn (Marsaglia algorithm)*/
/* https://rosettacode.org/wiki/Statistics/Normal_distribution */
double re_awgn[frame_size],im_awgn[frame_size];
double x,y,rsq,f;
for(i2=0;i2<frame_size;i2++)
{
do {
x = 2.0 * rand() / (double)RAND_MAX - 1.0;
y = 2.0 * rand() / (double)RAND_MAX - 1.0;
rsq = x * x + y * y;
}while( rsq >= 1. || rsq == 0. );
f = sqrt( -2.0 * log(rsq) / rsq );
re_awgn[i2] = sqrt(noise_var_1D)*x*f;
im_awgn[i2] = sqrt(noise_var_1D)*y*f;
}

/* channel output */
double re_chan_op[frame_size],im_chan_op[frame_size];
for(i2=0;i2<frame_size;i2++)
{
re_chan_op[i2] = re_sym[i2]+re_awgn[i2];
im_chan_op[i2] = im_sym[i2]+im_awgn[i2];
}
/*-------------------------------------------------------------*/
/*                        RECEIVER                  */
// Branch metrics for the VA
double BM[constell_size][frame_size],temp1,temp2;
int u,v;
for(u=0;u<constell_size;u++)
{
for(v=0;v<frame_size;v++)
{
temp1=(re_chan_op[v]-re_qpsk_sym[u])*(re_chan_op[v]-re_qpsk_sym[u]);
temp2=(im_chan_op[v]-im_qpsk_sym[u])*(im_chan_op[v]-im_qpsk_sym[u]);
BM[u][v]=temp1+temp2;
}    
}

// Viterbi algorithm (soft-input, hard-output)
int ip=0,a_hat[frame_size-decoding_delay];
int survivor_node[4][frame_size],survivor_ip[4][frame_size],loc=0,c,bk_cnt;
double path_metric[4][frame_size+1];
double temp3,temp4;
for(u=0;u<frame_size;u++)
{
for(v=0;v<num_states;v++)
{
temp3=path_metric[P_State[v][0]][u]+BM[Ga_Inx[v][0]][u]; 
temp4=path_metric[P_State[v][1]][u]+BM[Ga_Inx[v][1]][u];
if (temp3<=temp4)
{
    path_metric[v][u+1]=temp3;
    survivor_node[v][u]=P_State[v][0];
    survivor_ip[v][u]=P_Ip[v][0];
}
else 
{
    path_metric[v][u+1]=temp4;
    survivor_node[v][u]=P_State[v][1];
    survivor_ip[v][u]=P_Ip[v][1];
}
} //  for(v=0;v<num_states;v++) 

// Back tracing
if (u>=(decoding_delay-1))
{
for(c=1;c<num_states;c++)
if (path_metric[c][u+1] < path_metric[loc][u+1])
 loc = c;

for(bk_cnt=1;bk_cnt<=decoding_delay;bk_cnt++)
{
    ip=survivor_ip[loc][u+1-bk_cnt];
    loc=survivor_node[loc][u+1-bk_cnt];
}
a_hat[u-decoding_delay+1]=ip;
} //if (u>=(decoding_delay-1))
} // for(u=0;u<frame_size;u++)

int errors=0,err;
// calculating bit error rate
for (u=0;u<frame_size-decoding_delay;u++)
{
    err=a[u]==a_hat[u]?0:1;
    errors = errors+err;
}
C_BER = C_BER+errors;
} // for fr_cnt
system("cls"); /*clears the screen */
end = clock();
cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
printf("\n Total number of errors is %ld",C_BER);
printf("\n time is %lf sec",cpu_time_used);
}

