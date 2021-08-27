/*
 Antoniou Christodoulos 2641
 Ilias Diamantis 2685
 Tzounas Antonios 2368
 
 compiled and tested on gcc (Ubuntu 9.3.0-17ubuntu1~20.04) 9.3.0 
 and on gcc (Ubuntu 7.5.0-3ubuntu1~18.04) 7.5.0

compilation command:
 gcc provlima2.c -lm

 oi grafikes parastaseis exoun ginei xrisimopoiontas ta output arxeia .dat tou programmatos sto gnuplot
 arxika dialegoume ena katalilo terminal me tin entoli:
 set terminal png size 800,600
 dimiourgoume ena arxeio proorismou gia to paragomeno plot:
 set output "plot_image.png"
 kai telos, gia paradeigma oi grafikes parastaseis gia ta a b c d sto provlima 1 erotima delta me tin methodo tou euler 
 exoun ginei me tis entoles:
 plot "erotima_delta_euler.dat" using 1:2 linetype 7 linecolor 7 lw 2 with lines
 plot "erotima_delta_euler.dat" using 1:3 linetype 7 linecolor 7 lw 2 with lines
 plot "erotima_delta_euler.dat" using 1:4 linetype 7 linecolor 7 lw 2 with lines
 plot "erotima_delta_euler.dat" using 1:5 linetype 7 linecolor 7 lw 2 with lines
*/



#include <stdio.h>
#include <math.h>


#define M 1
#define g 9.81
#define Iz 0.08
#define AM 2685

double fz, Cz;
double h = 0.001;
double Kpz, Kdz, zdes;

// lisi
#define z(t) (-15.169)*exp((-0.234)*t)+0.169*exp((-20.988)*t)+13.425

// sistima
#define q1_d(t,a,b) b
#define q2_d(t,a,b) (((M*g+Kpz*(zdes-a)-Kdz*b)-g*M-Cz*b)/M)

// Euler
#define An1_Euler_d(t0,a0,b0) a0+h*q1_d(t0,a0,b0)
#define Bn1_Euler_d(t0,a0,b0) b0+h*q2_d(t0,a0,b0)

// Veltiomeni Euler
#define An1_Velt_Euler_d(t0,a0,b0) a0+(h/2)*(q1_d(t0,a0,b0)+q1_d(t0+h,a0+h*q1_d(t0,a0,b0),b0+h*q2_d(t0,a0,b0)))
#define Bn1_Velt_Euler_d(t0,a0,b0) b0+(h/2)*(q2_d(t0,a0,b0)+q2_d(t0+h,a0+h*q1_d(t0,a0,b0),b0+h*q2_d(t0,a0,b0)))



int main() {

    // gia na grapsoume ta apotelesmata se arxeia dat
    FILE *output_file;


    double t;
    double an1, a0, bn1, b0;
    double z0;


    Kpz = 5.0;
    Kdz = 15.0+((double)AM/1000.0);
    z0 = 0.0;
    zdes = ((double)AM/200.0);
    Cz = 3.0+((double)AM/5000.0);
    

    // Euler
    a0 = z0;
    b0 = 0.0;

    // neo arxeio
    output_file = fopen("./provlima_2_euler.dat","w+");
    if (output_file != NULL) fprintf(output_file,"%s\n","#   t       A            B");

    for (t=0.0; t<30.001; t+=h) {

        an1 = An1_Euler_d(t,a0,b0);
        bn1 = Bn1_Euler_d(t,a0,b0);
        //printf("%f - %f\n",an1,bn1);

        // grafoume ta dedomena
        if (output_file != NULL) fprintf(output_file,"%f %.10f %.10f\n",t,an1,bn1);

        // epomenes times a b
        a0 = an1;
        b0 = bn1;
    }

    // Veltiomeni Euler
    a0 = z0;
    b0 = 0.0;

    // neo arxeio
    output_file = fopen("./provlima_2_velt_euler.dat","w+");
    if (output_file != NULL) fprintf(output_file,"%s\n","#   t       A            B");

    for (t=0.0; t<30.001; t+=h) {

        an1 = An1_Velt_Euler_d(t,a0,b0);
        bn1 = Bn1_Velt_Euler_d(t,a0,b0);
        //printf("%f - %f\n",an1,bn1);

        // grafoume ta dedomena
        if (output_file != NULL) fprintf(output_file,"%f %.10f %.10f\n",t,an1,bn1);

        // epomenes times a b
        a0 = an1;
        b0 = bn1;
    }

    // lisi
    output_file = fopen("./provlima_2_lisi.dat","w+");
    if (output_file != NULL) fprintf(output_file,"%s\n","#   t       Z");

    for (t=0.0; t<30.001; t+=h) {

        z0 = z(t);

        if (output_file != NULL) fprintf(output_file,"%f %.10f\n",t,z0);
    }

    printf("%f\n", zdes);


    return 0;
}