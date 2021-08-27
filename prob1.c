/*
 Antoniou Christodoulos 2641
 Ilias Diamantis 2685
 Tzounas Antonios 2368
 
 compiled and tested on gcc (Ubuntu 9.3.0-17ubuntu1~20.04) 9.3.0 
 and on gcc (Ubuntu 7.5.0-3ubuntu1~18.04) 7.5.0
 
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

double fz, Tz, Cz, Cy;
double h = 0.001;
double Kpz, Kdz, Kpy, Kdy, zdes, ydes;

// sistima gia erotima vita
#define q1(t,a,b) b
#define q2(t,a,b) ((fz-g*M-Cz*fabs(b)*b)/M)
#define k1(t,c,d) d
#define k2(t,c,d) ((Tz-0.5*Cy*fabs(d)*d)/Iz)

// Euler gia erotima vita
#define An1_Euler(t0,a0,b0) a0+h*q1(t0,a0,b0)
#define Bn1_Euler(t0,a0,b0) b0+h*q2(t0,a0,b0)
#define Cn1_Euler(t0,c0,d0) c0+h*k1(t0,c0,d0)
#define Dn1_Euler(t0,c0,d0) d0+h*k2(t0,c0,d0)

// Veltiomeni Euler gia erotima vita
#define An1_Velt_Euler(t0,a0,b0) a0+(h/2)*(q1(t0,a0,b0)+q1(t0+h,a0+h*q1(t0,a0,b0),b0+h*q2(t0,a0,b0)))
#define Bn1_Velt_Euler(t0,a0,b0) b0+(h/2)*(q2(t0,a0,b0)+q2(t0+h,a0+h*q1(t0,a0,b0),b0+h*q2(t0,a0,b0)))
#define Cn1_Velt_Euler(t0,c0,d0) c0+(h/2)*(k1(t0,c0,d0)+k1(t0+h,c0+h*k1(t0,c0,d0),d0+h*k2(t0,c0,d0)))
#define Dn1_Velt_Euler(t0,c0,d0) d0+(h/2)*(k2(t0,c0,d0)+k2(t0+h,c0+h*k1(t0,c0,d0),d0+h*k2(t0,c0,d0)))

// sistima gia erotima delta
#define q1_d(t,a,b) b
#define q2_d(t,a,b) (((M*g+Kpz*(zdes-a)-Kdz*b)-g*M-Cz*fabs(b)*b)/M)
#define k1_d(t,c,d) d
#define k2_d(t,c,d) (((Kpy*(ydes-c)-Kdy*d)-0.5*Cy*fabs(d)*d)/Iz)

// Euler gia erotima delta
#define An1_Euler_d(t0,a0,b0) a0+h*q1_d(t0,a0,b0)
#define Bn1_Euler_d(t0,a0,b0) b0+h*q2_d(t0,a0,b0)
#define Cn1_Euler_d(t0,c0,d0) c0+h*k1_d(t0,c0,d0)
#define Dn1_Euler_d(t0,c0,d0) d0+h*k2_d(t0,c0,d0)

// Veltiomeni Euler gia erotima delta
#define An1_Velt_Euler_d(t0,a0,b0) a0+(h/2)*(q1_d(t0,a0,b0)+q1_d(t0+h,a0+h*q1_d(t0,a0,b0),b0+h*q2_d(t0,a0,b0)))
#define Bn1_Velt_Euler_d(t0,a0,b0) b0+(h/2)*(q2_d(t0,a0,b0)+q2_d(t0+h,a0+h*q1_d(t0,a0,b0),b0+h*q2_d(t0,a0,b0)))
#define Cn1_Velt_Euler_d(t0,c0,d0) c0+(h/2)*(k1_d(t0,c0,d0)+k1_d(t0+h,c0+h*k1_d(t0,c0,d0),d0+h*k2_d(t0,c0,d0)))
#define Dn1_Velt_Euler_d(t0,c0,d0) d0+(h/2)*(k2_d(t0,c0,d0)+k2_d(t0+h,c0+h*k1_d(t0,c0,d0),d0+h*k2_d(t0,c0,d0)))



int main() {

    // gia na grapsoume ta apotelesmata se arxeia dat
    FILE *output_file;


    double t;
    double an1, a0, bn1, b0, cn1, c0, dn1, d0;
    double z0 = ((double)AM/1000.0);
    double y0 = 0.0;
    Cz = 3.0-((double)AM/5000.0);
    Cy = 5.0-((double)AM/5000.0);


    // erotima vita - 1o zeugari eisodon
    fz = M*g+((double)AM/1000.0);
    Tz = 0.0;

    // Euler
    a0 = z0;
    b0 = c0 = d0 = 0.0;

    // neo arxeio
    output_file = fopen("./erotima_vita_zeugari_1_euler.dat","w+");
    if (output_file != NULL) fprintf(output_file,"%s\n","#   t       A            B            C            D");
    
    for (t=0.0; t<30.001; t+=h) {

        an1 = An1_Euler(t,a0,b0);
        bn1 = Bn1_Euler(t,a0,b0);
        cn1 = Cn1_Euler(t,c0,d0);
        dn1 = Dn1_Euler(t,c0,d0);
        //printf("%f - %f - %f - %f\n",an1,bn1,cn1,dn1);

        // grafoume ta dedomena
        if (output_file != NULL) fprintf(output_file,"%f %.10f %.10f %.10f %.10f\n",t,an1,bn1,cn1,dn1);
        
        // epomenes times a b c d
        a0 = an1;
        b0 = bn1;
        c0 = cn1;
        d0 = dn1;
    }

    // Veltiomeni Euler
    a0 = z0;
    b0 = c0 = d0 = 0.0;

    // neo arxeio
    output_file = fopen("./erotima_vita_zeugari_1_velt_euler.dat","w+");
    if (output_file != NULL) fprintf(output_file,"%s\n","#   t       A            B            C            D");

    for (t=0.0; t<30.001; t+=h) {

        an1 = An1_Velt_Euler(t,a0,b0);
        bn1 = Bn1_Velt_Euler(t,a0,b0);
        cn1 = Cn1_Velt_Euler(t,c0,d0);
        dn1 = Dn1_Velt_Euler(t,c0,d0);
        //printf("%f - %f - %f - %f\n",an1,bn1,cn1,dn1);

        // grafoume ta dedomena
        if (output_file != NULL) fprintf(output_file,"%f %.10f %.10f %.10f %.10f\n",t,an1,bn1,cn1,dn1);

        // epomenes times a b c d
        a0 = an1;
        b0 = bn1;
        c0 = cn1;
        d0 = dn1;
    }


    // erotima vita - 2o zeugari eisodon
    fz = M*g;
    Tz = ((double)AM/10000.0);

    // Euler
    a0 = z0;
    b0 = c0 = d0 = 0.0;

    // neo arxeio
    output_file = fopen("./erotima_vita_zeugari_2_euler.dat","w+");
    if (output_file != NULL) fprintf(output_file,"%s\n","#   t       A            B            C            D");

    for (t=0.0; t<30.001; t+=h) {

        an1 = An1_Euler(t,a0,b0);
        bn1 = Bn1_Euler(t,a0,b0);
        cn1 = Cn1_Euler(t,c0,d0);
        dn1 = Dn1_Euler(t,c0,d0);
        //printf("%f - %f - %f - %f\n",an1,bn1,cn1,dn1);

        // grafoume ta dedomena
        if (output_file != NULL) fprintf(output_file,"%f %.10f %.10f %.10f %.10f\n",t,an1,bn1,cn1,dn1);

        // epomenes times a b c d
        a0 = an1;
        b0 = bn1;
        c0 = cn1;
        d0 = dn1;
    }

    // Veltiomeni Euler
    a0 = z0;
    b0 = c0 = d0 = 0.0;

    // neo arxeio
    output_file = fopen("./erotima_vita_zeugari_2_velt_euler.dat","w+");
    if (output_file != NULL) fprintf(output_file,"%s\n","#   t       A            B            C            D");

    for (t=0.0; t<30.001; t+=h) {

        an1 = An1_Velt_Euler(t,a0,b0);
        bn1 = Bn1_Velt_Euler(t,a0,b0);
        cn1 = Cn1_Velt_Euler(t,c0,d0);
        dn1 = Dn1_Velt_Euler(t,c0,d0);
        //printf("%f - %f - %f - %f\n",an1,bn1,cn1,dn1);

        // grafoume ta dedomena
        if (output_file != NULL) fprintf(output_file,"%f %.10f %.10f %.10f %.10f\n",t,an1,bn1,cn1,dn1);

        // epomenes times a b c d
        a0 = an1;
        b0 = bn1;
        c0 = cn1;
        d0 = dn1;
    }


    // erotima delta
    Kpz = Kpy = 5.0;
    Kdz = 15.0+((double)AM/1000.0);
    Kdy = 20.0;
    z0 = 0.0;
    y0 = ((double)AM/10000.0);
    zdes = ((double)AM/200.0);
    ydes = -((double)AM/3000.0);
    Cz = 3.0+((double)AM/5000.0);
    Cy = 5.0;
    
    // Euler
    a0 = z0;
    c0 = y0;
    b0 = d0 = 0.0;

    // neo arxeio
    output_file = fopen("./erotima_delta_euler.dat","w+");
    if (output_file != NULL) fprintf(output_file,"%s\n","#   t       A            B            C            D");

    for (t=0.0; t<30.001; t+=h) {

        an1 = An1_Euler_d(t,a0,b0);
        bn1 = Bn1_Euler_d(t,a0,b0);
        cn1 = Cn1_Euler_d(t,c0,d0);
        dn1 = Dn1_Euler_d(t,c0,d0);
        //printf("%f - %f - %f - %f\n",an1,bn1,cn1,dn1);

        // grafoume ta dedomena
        if (output_file != NULL) fprintf(output_file,"%f %.10f %.10f %.10f %.10f\n",t,an1,bn1,cn1,dn1);

        // epomenes times a b c d
        a0 = an1;
        b0 = bn1;
        c0 = cn1;
        d0 = dn1;
    }

    // Veltiomeni Euler
    a0 = z0;
    c0 = y0;
    b0 = d0 = 0.0;

    // neo arxeio
    output_file = fopen("./erotima_delta_velt_euler.dat","w+");
    if (output_file != NULL) fprintf(output_file,"%s\n","#   t       A            B            C            D");

    for (t=0.0; t<30.001; t+=h) {

        an1 = An1_Velt_Euler_d(t,a0,b0);
        bn1 = Bn1_Velt_Euler_d(t,a0,b0);
        cn1 = Cn1_Velt_Euler_d(t,c0,d0);
        dn1 = Dn1_Velt_Euler_d(t,c0,d0);
        //printf("%f - %f - %f - %f\n",an1,bn1,cn1,dn1);

        // grafoume ta dedomena
        if (output_file != NULL) fprintf(output_file,"%f %.10f %.10f %.10f %.10f\n",t,an1,bn1,cn1,dn1);

        // epomenes times a b c d
        a0 = an1;
        b0 = bn1;
        c0 = cn1;
        d0 = dn1;
    }

    printf("zdes: %f ydes: %f\n", zdes, ydes);



    return 0;
}