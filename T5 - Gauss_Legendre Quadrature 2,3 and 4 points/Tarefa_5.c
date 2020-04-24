#include <stdio.h>
#include <stdlib.h> 
#include <math.h>

double f(double x){
    double a = ((sin(2*x)) + (4*pow(x,2)) + (3*x));
    return pow(a,2);
}

double X(double xi, double xf, double ak){
    return ((xi+xf)/2)+((xf-xi)/2)*ak;
}

double GaussL2(double a, double b){
    double c = -sqrt(1.0/3.0);
    double d = sqrt(1.0/3.0);
    return ((b-a)/2)*(f(X(a, b, c))+f(X(a, b, d)));
}

double GaussL3(double a, double b){
    double w1 = 5.0/9.0;
    double w2 = 8.0/9.0;
    double c = -sqrt(3.0/5.0);
    double d = 0;
    double e = sqrt(3.0/5.0);
    return ((b-a)/2)*((f(X(a, b, c))*w1)+(f(X(a, b, d))*w2)+(f(X(a, b, e))*w1));
}

double GaussL4(double a, double b){
    /*printf("a chegou=%f\n",a);
    printf("b chegou=%f\n",b);
    double dja = 5/6;
    printf("sqrt(5/6) =%f\n",sqrt(5.0/6.0));
    printf("sqrt(5/6)/6 =%f\n",sqrt(5/6)/6);
    */
    double w1 = ((1.0/2.0)-(sqrt(5.0/6.0)/6));
    double w2 = (1.0/36.0)*(18+sqrt(30.0));
    double c = -(sqrt((3+(2*sqrt(6.0/5.0)))/7));
    double d = -(sqrt((3-(2*sqrt(6.0/5.0)))/7));
    double e = -d;
    double g = -c;
    /*printf("primeira parte %f\n",((b-a)/2)*(f(X(a, b, c))*w1));
    printf("(b-a)/2 %f\n",((b-a)/2));
    printf("X(a, b, c) %f\n",X(a, b, c));
    printf("f(X(a, b, c) %f\n",f(X(a, b, c)));
    printf("f(X(a, b, c) *w1 %f\n",f(X(a, b, c))*w1);
    printf("w1=%f",w1);
    */


    return ((b-a)/2)*(f(X(a, b, c))*w1+f(X(a, b, d))*w2+f(X(a, b, e))*w2+f(X(a, b, g))*w1);
}


double integrate(short int degree,int a,int b){
    /*
    Arguments:
        degree - Legendre's polynomial degree 
        a e b - integration limits
    */
    
   if (degree > 4 || degree <2){
       printf("degree %d not implemented\n",degree);
       return 0 ;
   }

    double tolerance=10E-6;
    double result=0;
    double difference=1;
    double h = 0;
    int l = 0;
    double aux = 0;
    double xi = 0;
    int n = 0;
    
    while(difference>tolerance){
        l = 0;
        aux = result;
        result = 0;

        while(l<pow(2,n)){
            xi = (a+(l*(b-a))/pow(2,n));
            
            if (degree==2)
                    result += GaussL2(a+(l*(b-a)/pow(2,n)), a+((l+1)*((b-a)/pow(2,n))));
            else if(degree==3)
                    result += GaussL3(a+(l*(b-a)/pow(2,n)), a+((l+1)*((b-a)/pow(2,n))));
            else
                    result += GaussL4(a+(l*(b-a)/pow(2,n)), a+((l+1)*((b-a)/pow(2,n))));
            l++;    
            }    

        difference=abs(result-aux);
        n++;
        printf("\niteration %d = %f\n",n,result);
    }
    return result;
}

int main (){
    short int degree = 0;
    int a = 0;
    int b = 0;


    printf("Enter the Legendre's polynomial degree(2,3 or 4)\n");
    scanf("%hd", &degree);

    printf("\nEnter the a\n");
    scanf("%d", &a);

    printf("\nEnter the b\n");
    scanf("%d", &b);

    double result = integrate(degree,a,b);
    printf("Result = %f \n",result);

    return 0;
}