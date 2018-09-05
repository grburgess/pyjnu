#include <quadmath.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <execinfo.h>
#include <unistd.h>
#include <omp.h>



#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_synchrotron.h>
#include <gsl/gsl_errno.h>

#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>


/**  6th of March
  *  Stand alone version
  */

/**  **********************************
  *  **********************************
  *  Constant declaration
  *  **********************************
  *  **********************************
  */
/**  Constant CGS */
#define sigmaT 6.6524e-25
#define c 2.99792458e10
#define c2 8.9875517873681764e20
#define qe 4.8032068e-10
#define me 9.1093897e-28
#define mec2 8.1871111680068e-7
#define mp 1.6749286e-24
#define Msun 1.989e33
#define h 6.6260755e-27
#define hbar 1.05457266e-27
#define kb 1.380657e-16
#define sigmaSB 5.67051e-5
#define ath 7.5657e-15
#define Ggrav 6.67259e-8
#define Bcrit 4.414e13


/** Constant for Cosmology */
#define Hubble 69.6              /// km s^-1 Mpc^-1
#define Mpc_to_cm 3.08568e24               /// km s^-1 Mpc^-1
#define s_to_year 3.17098e-8

/** Constant for kinetic theory */
#define alphaf 0.0072973525664
#define lambdac 2.4263102367e-10
#define hplank 6.6261e-27
#define cyclosynchrotron 2.7992483604657531e6          /// e/(2 pi m_e c)   x B


/**  Constant m_e c^2 */
#define me_u 1.0
#define mmu_u 206.768289933



/** Numbers */
#define sqrt2 1.4142135623730950
#define sqrt3 1.7320508075688773
#define pi    3.1415926535897932
#define pi2   9.8696044010893586





/** Compilation */
/// gcc -Wall -o main main.c `pkg-config --cflags --libs gsl` -lm
double vabs(double x){
    if(x<0.0){
        return -x;
    }
    else{
        return x;
    }

}

double onemexpdexp(double x){

    if(vabs(x)<1e-3){
        return 1.0 + x + 0.5*x*x + 0.166666666666667*x*x*x;
    }
    else{
        return (1.0 - exp(-x))/x;
    }


}


/**  **********************************
  *  **********************************
  *  Code specific declaration
  *     1- structure
  *     2- grid properties
  *  **********************************
  *  **********************************
  */
#define GRID_E 100
#define GRID_P 200

typedef struct particle{


    /** Energy grid */
    double *E;
    double *Eb;

    /** Distribution function, or coefficients */
    double *f;
    double *f1;
    double *f2;

    /** Characteristics of the particle type */
    int dim;
    double emin;
    double emax;

    double mass;
}particle;



/** Grid initialization */
int make_particle(particle *pt, double emin, double emax, int dim, double mass){

    int i;
    double step = exp(log(emax/emin)/((double) dim));
    double sqrt_step = sqrt(step);

    pt->emin = emin;
    pt->emax = emax;
    pt->dim = dim;

    /** Energy grid */
    pt->E = calloc(dim, sizeof(double));
    pt->Eb = calloc(dim+1, sizeof(double));

    /** Distribution function, or coefficients */
    pt->f = calloc(dim, sizeof(double));
    pt->f1 = calloc(dim, sizeof(double));
    pt->f2 = calloc(dim, sizeof(double));

    pt->mass = mass;
    pt->Eb[0] = emin/sqrt_step;
    for(i=0;i<dim;i++){
        pt->E[i] = emin*pow(step,(double) i);
        pt->Eb[i+1] = pt->E[i]*sqrt_step;
    }


    return 0;
};

int remove_grid(particle *pt){

    free(pt->E);
    free(pt->Eb);
    free(pt->f);
    return 0;

}

/** Initialize the distribution function*/
int init_pl(particle *pt, double emin, double emax, double p, double normalization){

    int i;
    for(i=0; i<pt->dim; i++){
        if(pt->Eb[i] > emin && pt->Eb[i+1]< emax){
            pt->f[i] = normalization*( pt->Eb[i]/((p-1)*pow(pt->Eb[i],p)) - pt->Eb[i+1]/((p-1)*pow(pt->Eb[i+1],p)) );
        }
        else{
            pt->f[i] = 0.0;
        }
    }
    return 0;

};



/**  **********************************
  *  **********************************
  *  Synchrotron function
  *     1- in principle only read from file
  *     2- code for computation of the emissivity is still in place
  *  **********************************
  *  **********************************
  */
/** Renormalized emissivities, see notes*/
double Synchrotron_emissivity(double y, double g){

    double ythreeg2 = y/(3.0*g*g);



    if(ythreeg2 < 705){

        double K43 = gsl_sf_bessel_Knu(4.0/3.0, ythreeg2);
        double K13 = gsl_sf_bessel_Knu(1.0/3.0, ythreeg2);
        double first = K43*K13;
        double second = 0.6*ythreeg2*(K43*K43-K13*K13);
        //printf("ythreeg2 = %g\tresult = %g\t K43 = %g\t K13 = %g\n",ythreeg2, ythreeg2*ythreeg2*(first - second), K43, K13);
        return ythreeg2*ythreeg2*(first - second);

    }
    else{
        return 0;
    }
}

/** Function to be integrated
  * *param = x
  */
double Integral_emissivity(double g, double *param){
    return Synchrotron_emissivity(*param , g);
}


/**  3rd of March 2018
  *  Computation of the opacity
  */
double diff_emissivity(double y, double g){

    double x = y/(3.0*g*g);
    if(x < 705.0){

        double K43 = gsl_sf_bessel_Knu(4.0/3.0, x);
        double K13 = gsl_sf_bessel_Knu(1.0/3.0, x);
        double first = K43*K13;
        double second = 1.2*x*(K43*K43-K13*K13);
        //printf("ythreeg2 = %g\tresult = %g\t K43 = %g\t K13 = %g\n",ythreeg2, ythreeg2*ythreeg2*(first - second), K43, K13);
        return -x*y*2.0/(3.0*g*g*g)*(first - second);
    }
    else{
        return 0.0;
    }

}

double Integral_opacity(double g, double *param){
    return 2.0*g/sqrt(g*g-1.0)*Synchrotron_emissivity( *param, g) + diff_emissivity( *param, g);
    //Synchrotron_emissivity(*param , g, param);
}




/**  3rd of March 2018
  *  Computation of the emissivity
  */
int sync_emissivity(double **syncrate, particle *ph, particle *elec){

    /**  If the file exist, just read.
      */

    int i,j;
    double a;
    FILE *intput;
    intput = fopen("synchrotron_rate.txt", "r");
    while(fscanf(intput, "%d%d%lf", &i,&j,&a) != EOF){
        syncrate[i][j] = a;
    }
    fclose(intput);


    /**  If not, create */
/*
    int i, j;
    double *a;
    double param;
    a = (double *) calloc(2, sizeof(double));

    FILE *output;
    output = fopen("synchrotron_rate.txt", "a");

    double result;
    for(i=0;i<GRID_P;i++){
        for(j=0;j<GRID_E;j++){
            a[0] = elec->Eb[j];
            a[1] = elec->Eb[j+1];;
            param = ph->E[i];
            printf("i = %d\t j = %d\tresult = %g\n",i,j, result);
            result = gauss_kronrod_7_15( a[0], a[1], 1e-4, &Integral_emissivity, &param, 0 );

            syncrate[i][j] = result;
            //printf(" i = %d\t j = %d\t x = %g - %g\t g = %g - %g\tresult = %g\n", i, j, a[0], a[1], b[0], b[1], result);
            fprintf(output, "%d\t%d\t%g\n",i,j,result);
            //getchar();
        }
        //getchar();
    }
    fclose(output);
    printf("result = %g\n",result);


    free(a);*/
    return 0;

}



/**  3rd of March 2018
  *  Computation of the opacity
  */
int opacity(double **opacity_rate, particle *ph, particle *elec){

    /**  If the file exist, just read.
      */

    int i,j;
    double a;
    FILE *intput;
    intput = fopen("tau.txt", "r");
    while(fscanf(intput, "%d%d%lf", &i,&j,&a) != EOF){
        opacity_rate[i][j] = a;
    }
    fclose(intput);


    /**  If not, create */
    /*
    int i, j;
    double *a;
    double param;
    a = (double *) calloc(2, sizeof(double));

    FILE *output;
    output = fopen("tau.txt", "a");

    double result;
    for(i=0;i<GRID_P;i++){
        for(j=0;j<GRID_E;j++){
            a[0] = elec->Eb[j];
            a[1] = elec->Eb[j+1];;
            param = ph->E[i];
            printf("i = %d\t j = %d\tresult = %g\n",i,j, result);
            result = gauss_kronrod_7_15( a[0], a[1], 1e-4, &Integral_opacity, &param, 0 );

            opacity_rate[i][j] = result;
            //printf(" i = %d\t j = %d\t x = %g - %g\t g = %g - %g\tresult = %g\n", i, j, a[0], a[1], b[0], b[1], result);
            fprintf(output, "%d\t%d\t%g\n",i,j,result);
            //getchar();
        }
        //getchar();
    }
    fclose(output);
    printf("result = %g\n",result);
    free(a);
    */
    return 0;
}


/**  **********************************
  *  **********************************
  *  Same for pairs
  *  **********************************
  *  **********************************
  */
double phipair(double x){
    double v = x -1.0;
    double sqrtv = sqrt(v);
    double sqrtx = sqrt(x);
    double w = (sqrtx+sqrtv)/(sqrtx-sqrtv);
    double wp1 = w+1.0;
    double logw = log(w);
    double logw2 = logw*logw;
    double logwp1 = log(w+1.0);
    double logwp12 = logwp1*logwp1;
    double dilog = gsl_sf_dilog(1.0/wp1);
    //printf("x = %g\t dilog = %g\n",x, dilog);
    return (2.0*v*x+1.0)*logw/x - 2.0*(v+x)*sqrtv/sqrtx - logw2 + 2.0*logwp12 + 4.0 * dilog  - pi2/3.0;
}

double Rpp(double x){
    if(x>1.0){
        return (3.0*sigmaT)*phipair(x)/(8.0*x*x);
    }
    else{
        return 0.0;
    }
}


/**  **********************************
  *  **********************************
  *  Inverse Compton scattering
  *     Shall we use the full kernel instead of the Jones approximation?
  *  **********************************
  *  **********************************
  */
double Jones_head_on(double x1, double x, double g){   /// x_1 to x

    double q = x/(4.0*x1*g*g*(1.0-x/g));
    if(q>0.25/g/g && q<=1.0){
        double fqG = 2.0*q*log(q)+(1.0+2.0*q)*(1.0-q)+0.5*(16.0*x1*g*q*x1*g*q)*(1.0-q)/(1.0+4.0*x1*g*q);
        return fqG/(x1*g*g);
    }
    else{
        return 0.0;
    }
}


int main (){

    clock_t begin = clock();


    /** Define the physics  here */
    /// Region parameters, taken from Saugé and Henri (2004)
    double Bfield = 0.1;
    double R = 2.1e15;
    double n = 13.5; //1e-10/(R*sigmaT);

    double nuB = 2.7992483604657531e6*Bfield;
    double A = 4.0*pi*sqrt3*qe*qe*nuB/c;

    double tauCst = R/(8.0*pi*me*nuB*nuB);


    particle photon;
    particle electron;
    double emin = 0.01;
    double emax = 1e30;
    make_particle(&photon, emin, emax, GRID_P, 0.0);
    make_particle(&electron, 5.0, 1e8, GRID_E,0.0);
    init_pl(&electron, 1.1e3, 2.5e6, 2.5, n);


    int i,j,k;
    double **syncrate;
    syncrate = (double **) malloc(GRID_P*sizeof(double *));
    for(i = 0 ; i< GRID_P; i++){
        syncrate[i] = (double *) calloc(GRID_E, sizeof(double));
    }
    sync_emissivity(syncrate, &photon, &electron);

    double **opacity_rate;
    opacity_rate = (double **) malloc(GRID_P*sizeof(double *));
    for(i = 0 ; i< GRID_P; i++){
        opacity_rate[i] = (double *) calloc(GRID_E, sizeof(double));
    }
    opacity(opacity_rate, &photon, &electron);

    FILE *distrib;
    distrib = fopen("electron.txt","a");
    for(i = 0;i< GRID_E;i++){
        fprintf(distrib, "%g\t%g\n", electron.E[i],electron.f[i]);
    }
    fclose(distrib);


    distrib = fopen("photon_abs.txt","a");
    for(i = 0; i< GRID_P;i++){
        double summ = 0.0;
        for(j = 0;j< GRID_E;j++){
            summ += A*syncrate[i][j]*electron.f[j];
        }
        photon.f[i] = summ;
    }
    /**
      *  Note j = h epsilon n / 4 pi = h*h*nu n/4 pi mec2  => n = 4 pi j mec2 /(h^2 nu)
      */
    for(i = 0; i <GRID_P ;i++){
        double summtau = 0.0;
        double summIC = 0.0;
        double summPP = 0.0;
        for(j = 0;j< GRID_E;j++){
            summtau += opacity_rate[i][j]*electron.f[j];
            for(k = 0; k< GRID_P;k++){
                //(4.0*pi*mec2/(h*h*photon.E[k]*nuB))
                double result = photon.E[i] < photon.E[k]*electron.Eb[j]/5.0 ? 0 : electron.f[j]*(electron.Eb[j+1]-electron.Eb[j])*Jones_head_on(h*photon.E[k]*nuB/mec2, h*photon.E[i]*nuB/mec2, electron.E[j])*
                        photon.f[k]*(4.0*pi/(h*c*nuB*photon.E[k]))*nuB*(photon.Eb[k+1]-photon.Eb[k])/mec2;
                summIC += result;
                //printf("prefac = %g\tsummIC = %g\n", (4.0*pi*mec2/(h*h*photon.E[k]*nuB)), summIC);
                //if(summIC != summIC){
                //    getchar();
                //}
            //printf("summIC = %g\n", summIC);
            //    if(summIC > 0.0){
            //        printf("j = %d\t k = %d\tEk = %g\tEi = %g\tsumm = %g\n",j,k, photon.E[k], photon.E[i], result);
            //    }

            }
            /*if(summIC > 0.0){
                getchar();
            }*/

            summPP += Rpp(photon.E[i]*photon.E[j]*h*h*nuB*nuB/mec2/mec2)*photon.f[j]*4.0*pi*mec2/(h*h*photon.E[j]*nuB)*nuB*(photon.Eb[j+1]-photon.Eb[j]);
        }

        printf("summIC = %g\tphoton energy = %g\tsumPP = %g\n", summIC, h*photon.E[i]/mec2, summPP*R);
        //getchar();
        summIC = summIC*h*c*(photon.E[i]*nuB)*0.75*c*sigmaT/(4.0*pi);
        summtau = A*tauCst*summtau/(photon.E[i]*photon.E[i]);
        photon.f2[i] = photon.f[i] + summIC;
        photon.f1[i] = photon.f[i];
        photon.f[i] = (pi*R*R*25*25*25/(134.1*Mpc_to_cm*134.1*Mpc_to_cm))*photon.f[i] * onemexpdexp(summtau);//*onemexpdexp(summPP*R);
        photon.f2[i] = (pi*R*R*25*25*25/(134.1*Mpc_to_cm*134.1*Mpc_to_cm))*photon.f2[i] * onemexpdexp(summtau);//*onemexpdexp(summPP*R);


        //photon.f[i] = summ;
        fprintf(distrib,"%g\t%g\t%g\t%g\n", photon.E[i]*nuB*25.0,R*photon.f[i], R*photon.f2[i], summPP);

    }
    printf("nuB = %g\n", nuB);
    fclose(distrib);

    for(i = 0 ; i< GRID_P; i++){
        free(syncrate[i]);
        free(opacity_rate[i]);
    }
    free(syncrate);
    free(opacity_rate);

    remove_grid(&photon);
    remove_grid(&electron);



    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("time spent = %g\n", time_spent);

	return 0;
}
