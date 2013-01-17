/* Program to run RK5 on a Duffing oscillator
 * File:   duffing.cpp
 * Author: Arjun
 *Run on an Intel Core i7 machine
 * gcc 4.6.1 compiler, NetBeans IDE
 */


#include<cstdio>
#include<cmath>
#include<conio.h>
#include<cstdlib>
/*----All commented out printfs were used for debugging purposes----*/
/*---------------------------------RK5------------------------------*/
double ode5(double yold[], double ynew[], double &h, double t, 
		//&h - Call by reference - Increments h in main
        int numberofequations, void frhs(double[], double, double[])) {
    int i;
    double *k1, *k2, *k3, *k4, *k5, *k6, *temp, *ynewstar; //*yerror;  //pointers to arrays

    //Values from Numerical Recipes

    double c2 = 0.2, c3 = 0.3, c4 = 0.8, c5 = 8.0 / 9.0, //Constants
            a21 = 0.2, //Coefficients
            a31 = 3.0 / 40.0, a32 = 9.0 / 40.0,
            a41 = 44.0 / 45.0, a42 = -56.0 / 15.0, a43 = 32.0 / 9.0,
            
            a51 = 19372.0 / 6561.0, a52 = -25360.0 / 2187.0, 
            a53 = 64448.0 / 6561.0, a54 = -212.0 / 729.0,
            
            a61 = 9017.0 / 3168.0, a62 = -355.0 / 33.0, a63 = 46732.0 / 5247.0, 
            a64 = 49.0 / 176.0, a65 = -5103.0 / 18656.0,
            
            a71 = 35.0 / 384.0, a73 = 500.0 / 1113.0, a74 = 125.0 / 192.0, 
            a75 = -2187.0 / 6784.0, a76 = 11.0 / 84.0;
    /*Error value coefficients*/ 
    double e1 = 71.0 / 57600.0, e3 = -71.0 / 16695.0, e4 = 71.0 / 1920.0, 
            e5 = -17253.0 / 339200.0, e6 = 22.0 / 525.0, e7 = -1.0 / 40.0;
    double *error, *scale; //Scale now
    double hmin = 1.0e-10; 
    //Minimum size to avoid getting stuck in an infinite loop
    
    double tolerable_error = 1.0e-7; //Error limit
    double error_min = 1.0e-6; //Error limit for scale
    double max_error_i; //To find maximum error
    k1 = new double[10 * numberofequations]; 
    //integer*numberofequations where the integer=# of arrays being used
    
    if (NULL == k1) {
        printf("Cannot allocate k1 in ode5() \n");
        return (0);
    }
    k2 = k1 + numberofequations;
    k3 = k2 + numberofequations;
    k4 = k3 + numberofequations;
    k5 = k4 + numberofequations;
    k6 = k5 + numberofequations;
    temp = k6 + numberofequations;
    ynewstar = temp + numberofequations; //For the embedded 4th order
    error = ynewstar + numberofequations; //For the error array
    scale = error + numberofequations; //For the scale array
    frhs(yold, t, k1); // Get the RHS
    //printf("\n Calling frhs (yold,t,k1) \n");


repeat_current_step: //goto label
    //printf("h = %f\n",h);
    double h_remember = h;
    for (i = 0; i < numberofequations; i++)//Step 1
    {
        temp[i] = yold[i] + a21 * h * k1[i];
    }
    frhs(temp, t + c2*h, k2); //Step 2
    //printf("\n Calling frhs (yold,t+c2*h,k2) \n");
    for (i = 0; i < numberofequations; i++) {
        temp[i] = yold[i] + h * (a31 * k1[i] + a32 * k2[i]);
    }
    //printf("\n Calling frhs (yold,t+c3*h,k3) \n");
    frhs(temp, t + c3*h, k3); //Step 3
    for (i = 0; i < numberofequations; i++) {
        temp[i] = yold[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
    }
    //printf("\n Calling frhs (yold,t+c4*h,k4) \n");
    frhs(temp, t + c4*h, k4); //Step 4
    for (i = 0; i < numberofequations; i++) {
        temp[i] = yold[i] + h * (a51 * k1[i] + a52 * k2[i] + 
                a53 * k3[i] + a54 * k4[i]);
    }
    //printf("\n Calling frhs (yold,t+c5*h,k5) \n");
    frhs(temp, t + c5*h, k5); //Step 5
    for (i = 0; i < numberofequations; i++) {
        temp[i] = yold[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] 
                                + a64 * k4[i] + a65 * k5[i]);
    }
    //printf("\n Calling frhs (yold,t+c6*h,k6) \n");
    frhs(temp, t + h, k6); //Step 6
    for (i = 0; i < numberofequations; i++) {
        ynew[i] = yold[i] + h * (a71 * k1[i] + a73 * k3[i] + a74 * k4[i] 
                                + a75 * k5[i] + a76 * k6[i]);
    }

    //printf("\n Final \n");
    frhs(ynew, t + h, ynewstar); //Final

    for (i = 0; i < numberofequations; i++)//to set up the scale
        scale[i] = fabs(yold[i]) + fabs(h * k1[i]) + fabs(error_min);

    for (i = 0; i < numberofequations; i++)//Estimate error
    {
        error[i] = (h * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] 
                + e5 * k5[i] + e6 * k6[i] + e7 * ynewstar[i])) / scale[i];
        //printf("Error(%d) = %g\n",i,error[i]);
    }

    max_error_i = fabs(error[0]);
    for (i = 0; i < numberofequations; i++) {
        if (fabs(error[i]) > max_error_i) {
            max_error_i = fabs(error[i]);
        }
    }//Finding the maximum error

    if (max_error_i > tolerable_error) {
        h = h / 5;

        if (h < hmin) exit(0);

        goto repeat_current_step;
    }

    if (max_error_i < tolerable_error) {
        h = h * (pow(tolerable_error / max_error_i, 0.2));
    }
    //printf("h value in ode5 = %f\n",h);

    //getch();
    //exit(0);

    delete[]k1; //Free memory
    return (h_remember);
}//end ode5
//end rk5
const int N = 2; /*For Duffing oscillator, SHM and Mathieu system*/

//evaluate the ODE RHS
//Change for each problem
/*Old systems for testing*/
void myrhs1(double yold[], double t, double f[])//Harmonic oscillator
{
    f[1] = -yold[0]; //Harmonic oscillator with omega = 1 //f[1] = y'(t)
    f[0] = yold[1]; //f(0) = y(t)
}

void myrhs(double yold[], double t, double f[])
//Chaotic system - Mathieu equation
{
    f[1] = -(cos(2 * t) * yold[0]);
    f[0] = yold[1];
}

/*--------------Duffing Oscillator-------------*/
void duffing(double yold[], double t, double f[]) {
    /*Parameters :
     * alpha    -       a
     * beta     -       b
     * damping coefficient      -       mu
     * mass     -       m
     * periodicity coefficient  -       A
     * frequency        -       w       (omega)
     */


    double b = 1, a = 1, mu = 0.1, A = 0.35, w = 1.4, m = 1;
    f[1] = (-b * yold[0] - a * yold[0] * yold[0] * yold[0] 
                - mu * yold[1] + A * cos(w * t)) / m;
    
    f[0] = yold[1]; //f(0) = y(t)
}

//Main program

int main() {
    int istep; //nstep ;
    double yold[N], ynew[N];
    double t, t_initialize = 0.0; //initial value of t
    double tmax;
    double h = 1; //step size
    int n_period = 100;
    double h_increment;

    FILE *fp, *fp1;
    t = t_initialize;
    double pi = 4 * atan(1);
    double w = 1.4;
    tmax = 2 * pi * n_period / w;

    yold[0] = 0; //y(0)=1
    yold[1] = 0; //y'(0)=0
    /* To open a new output file and give it a title */
    fp = fopen("Duffing.dat", "w+");
    if (NULL == fp) {
        printf(" Unable to open the file \n");
        return ( 0);

    }
    fp1 = fopen("PoincareDuffing.dat", "w+");
    if (NULL == fp1) {
        printf(" Unable to open the file \n");
        return ( 0);

    }
    //printf("\n Starting process \n\n");
    int time = 0;

    for (istep = 0; t < tmax; istep++) {
        
        //printf("istep = %d \n",istep);
        h_increment = ode5(yold, ynew, h, t, N, duffing); 
        //Uses call by reference
        
        //printf("t,ynew[0],ynew[1]=%f,%f,%f\n",t,ynew[0],ynew[1]);
        fprintf(fp, "%f %f %f\n", t, ynew[0], ynew[1]);

       for (int i = 0; i < N; i++) 
       {
            yold[i] = ynew[i];
       }
       
        //Print out ynew and copy ynew-->yold
        t = t + h_increment;
        //printf("h in main = %f\n",h);
    }

    /*For the Poincare section
     * Set up a loop with tmin = 0 and tmax = one period ( 1*T)
     * Perform RK5 and get the values
     * Loop it such that tmin = 1*T and tmax = 2T 
     * Perform RK5 and get the values
     * Keep doing this for as many time periods as needed
     */


    double tbound = 2 * pi / w;
    /*Set the initial conditions again*/
    t = t_initialize;
    yold[0] = 0; //y(0)=1
    yold[1] = 1; //y'(0)=0
    int count = 0;

    while (t < tmax) {
        for (istep = 0; t < tbound; istep++) 
        {
            //printf("istep = %d \n",istep);
            h_increment = ode5(yold, ynew, h, t, N, duffing);
            //printf("t,ynew[0],ynew[1]=%f,%f,%f\n",t,ynew[0],ynew[1]);

            for (int i = 0; i < N; i++) 
            {
                yold[i] = ynew[i];
            }
                 
            t = t + h_increment;
          
        }
        fprintf(fp1, "%f %f %f\n", t, ynew[0], ynew[1]);
        //printf("%d\t%f\n",count,tbound);
        count = count + 1;
        tbound = tbound + 2 * pi / w;
    }
    return (0);
    //end
}


