/* Lorenz and Rossler systems
 * File:   lorenz.cpp
 * Author: Arjun
 *
 * Created on May 1, 2012, 8:45 PM
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
    double *k1, *k2, *k3, *k4, *k5, *k6, *temp, *ynewstar; 
    //*yerror;  //pointers to arrays
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
    double hmin = 1.0e-10; //Minimum size to avoid infinite loop
    double tolerable_error = 1.0e-6; //Error limit
    double error_min = 1.0e-2; //Error limit for scale
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
        temp[i] = yold[i] + h * (a51 * k1[i] + a52 * k2[i] 
                        + a53 * k3[i] + a54 * k4[i]);
    }
    //printf("\n Calling frhs (yold,t+c5*h,k5) \n");
    frhs(temp, t + c5*h, k5); //Step 5
    for (i = 0; i < numberofequations; i++) {
        temp[i] = yold[i] + h * (a61 * k1[i] + a62 * k2[i] 
                                    + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
    }
    //printf("\n Calling frhs (yold,t+c6*h,k6) \n");
    frhs(temp, t + h, k6); //Step 6
    for (i = 0; i < numberofequations; i++) {
        ynew[i] = yold[i] + h * (a71 * k1[i] + a73 * k3[i] 
                        + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
    }

    //printf("\n Final \n");
    frhs(ynew, t + h, ynewstar); //Final

    for (i = 0; i < numberofequations; i++)//to set up the scale
        scale[i] = fabs(yold[i]) + fabs(h * k1[i]) + fabs(error_min);

    for (i = 0; i < numberofequations; i++)//Estimate error
    {
        error[i] = (h * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + 
                e5 * k5[i] + e6 * k6[i] + e7 * ynewstar[i])) / scale[i];
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
const int N = 3; /*For Lorenz system*/

void lorenz(double x[], double t, double dx[]) {
    double R = 28, B = 8.0 / 3.0, Sigma = 10.0;

    dx[0] = Sigma * (-x[0] + x[1]);
    dx[1] = R * x[0] - x[1] - x[0] * x[2];
    dx[2] = -B * x[2] + x[0] * x[1];
}

void rossler(double x[], double t, double dx[]) {
    double a = 0.2, b = 0.2, c = 10.0;

    dx[0] = -x[1] - x[2];
    dx[1] = x[0] + a * x[1];
    dx[2] = b + x[0] * x[2] - c * x[2];
}
//Main program
int main() {
    int istep; //nstep ;
    double x[N], dx[N];
    double x1[N], dx1[N];
    
    //To check for sign change in y for Poincare section
    double flagx[N], flagt, flagdx[N]; 



    /* x[N] -->  x , y , z
     * dx[N] -->  dx/dt , dy/dt , dz/dt 
     */

    double t, t_initialize = 0.0; //initial value of t
    double tmax;
    double h = 0.5; //step size

    double h_increment;

    FILE *fp, *fp1, *fp2;
    t = t_initialize;

    tmax = 5e2;

    /*Lorenz plot*/
    x[0] = 0.0; //x0
    x[1] = 1.0; //y0
    x[2] = 0.0; //z0

    /*Rossler plot*/
    x1[0] = 0.0; //x0
    x1[1] = 0.0; //y0
    x1[2] = 0.0; //z0



    /* To open a new output file and give it a title */
    fp = fopen("Lorenz.dat", "w+");
    if (NULL == fp) {
        printf(" Unable to open the file \n");
        return ( 0);

    }

    fp1 = fopen("PoincareLorenz.dat", "w+");
    if (NULL == fp1) {
        printf(" Unable to open the file \n");
        return ( 0);
    }

    fp2 = fopen("Rossler.dat", "w+");
    if (NULL == fp2) {
        printf(" Unable to open the file \n");
        return ( 0);

    }
    //printf("\n Starting process \n\n");
    for (istep = 0; t < tmax; istep++) {

        for (int i = 0; i < N; i++) {
            flagx[i] = x[i];
            flagdx[i] = dx[i]; //Change to x1,dx1 for Rossler
            flagt = t;
        }
        //printf("istep = %d \n",istep);
        h_increment = ode5(x, dx, h, t, N, lorenz); //Uses call by reference
        //printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 
        //        t , x[0], x[1], x[2],dx[0] , dx[1] , dx[2]);
        fprintf(fp, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 
                        t, x[0], x[1], x[2], dx[0], dx[1], dx[2]);





        h_increment = ode5(x1, dx1, h, t, N, rossler); //Uses call by reference
        fprintf(fp2, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 
                        t, x1[0], x1[1], x1[2], dx1[0], dx1[1], dx1[2]);



        for (int i = 0; i < N; i++) {
            x[i] = dx[i];
            x1[i] = dx1[i];
        }
        
        /*For Poincare section where y changes sign (crosses zero)*/
        if (((flagx[1] < 0) && (x[1] > 0)) || ((flagx[1] > 0) && (x[1] < 0))) 
        {
            //Change to x1,dx1 for Rossler
            fprintf(fp1,"%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
                    0.5 * (t + h_increment + flagt),
                    0.5 * (flagx[0] + x[0]),
                    0.5 * (flagx[1] + x[1]), 0.5 * (flagx[2] + x[2]),
                    0.5 * (flagdx[0] + dx[0]), 0.5 * (flagdx[1] + dx[1]),
                    0.5 * (flagdx[2] + dx[2]));
        }
        //Print out ynew and copy ynew-->yold
        t = t + h_increment;
        //printf("h in main = %f\n",h);
    }
    return (0);
    //end
}



