#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <time.h>

#define pi 3.141592653589793238462643383279502884197169399375105820974

double* aloc(int size){
    double* vector = (double*) calloc(size,sizeof(double));
    return vector;
}

double f(double x){
    return 1/(pow((pow(x,2)),(1.0/3.0)));
    //return 1.0/ (sqrt(4.0-pow(x,2)));
}

double X(double xi,double xf,double ak){
    return ((xi+xf)/2.0)+((xf-xi)/2.0)*ak;
}

double xs_simple(double s,double a, double b){
    return ((a+b)/2.0)+(((b-a)/2.0)*tanh(s));
}

double xs_double(double s,double a,double b){
    return ((a+b)/2.0)+(((b-a)/2.0)*tanh((pi/2.0) * sinh(s)));
}

double exp_simple(double s,double a,double b){
    return (f(xs_simple(s,a,b)))*(((b-a)/2.0)*(1/pow(cosh(s),2)));
}

double exp_double(double s,double a, double b){
    return (f(xs_double(s,a,b)))*(((b-a)/2.0)*((pi/2.0)*(cosh(s)/pow((cosh((pi/2.0)*sinh(s))),2))));
}

double GaussL4(double a,double b,double lim_inf,double lim_sup,short int simple){
    double w1 = (1.0/2.0)-(sqrt(5.0/6.0))/6.0;
    double w2 = (1.0/36.0)*(18+sqrt(30.0));
    double c = -(sqrt((3+(2*sqrt(6.0/5.0)))/7.0));
    double d = -(sqrt((3-(2*sqrt(6.0/5.0)))/7.0));
    double e = -d;
    double g = -c;

    switch (simple){
    case 1: 
        return ((b-a)/2.0)*((exp_simple(X(a, b, c),lim_inf,lim_sup))*w1+(exp_simple(X(a, b, d),lim_inf,lim_sup))*w2+(exp_simple(X(a, b, e),lim_inf,lim_sup))*w2+(exp_simple(X(a, b, g),lim_inf,lim_sup))*w1);
    break;

    default: 
        return ((b-a)/2.0)*((exp_double(X(a, b, c),lim_inf,lim_sup))*w1+(exp_double(X(a, b, d),lim_inf,lim_sup))*w2+(exp_double(X(a, b, e),lim_inf,lim_sup))*w2+(exp_double(X(a, b, g),lim_inf,lim_sup))*w1);
    break;
    }

  }


double integrate(int c,double lim_inf,double lim_sup,short int flag_exp){
    /*
    Arguments:
        c - integration limits of the numeric formula
        lim_inf - inferior limit of the integral(a)
        lim_sup - superior limit of the integral(b)
        flag_exp - if 1: Simple Exponential,else: double exponential
    */   
    
    double tolerance=10E-4;
    double result=0;
    double aux = 0,aux_c = 0;
    double difference=1;
    int l = 0,n = 0,a = 0,b = 0;

    printf("\n\nc = %d",c);
    aux_c=result;    
    result=0;
    difference=1;
    n=0;
    a=-c;                       
    b=c;
    while(difference>tolerance){
        l = 0;
        aux = result;
        result = 0;

        while(l<pow(2,n)){
            result += GaussL4(a+((l*(b-a))/pow(2,n)), a+((l+1)*((b-a)/pow(2,n))),lim_inf,lim_sup,flag_exp);
            ++l;  
            }    

        difference=fabs(result-aux);
        ++n;
        printf("\niteration %d = %f",n,result);
        //if (n==30) break;  //Depending on the c and the function(f), you might have to use this  
    }
    return result;
}


int main (){
    short int method = 0;
    double lim_inf = 0,lim_sup = 0;
    int c = 0;
    double result;
    printf("[1]Simple Exponential \n[0]Double Exponential\n");
    scanf("%hd", &method);

    if (method > 1 || method < 0){
        printf("Invalid Digit\n");
        return 0;
    }

    printf("Enter the inferior limit(a)\n");
    scanf("%lf",&lim_inf);

    printf("Enter the superior limit(b)\n");
    scanf("%lf",&lim_sup);

    printf("Enter the c\n");
    scanf("%d",&c);
        
    clock_t t = clock();

    result=integrate(c,lim_inf,lim_sup,method);
    printf("\nResult=%f",result);

    t = clock() - t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC; 
    printf("\nc = %d - time taken =  %.5lf segundos \n",c,time_taken / 3.0);

    return 0;
}
