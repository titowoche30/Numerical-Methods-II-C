#include <stdio.h>
#include <stdlib.h> 
#include <math.h>

double int_closed_degree1(double a,double b,double h){
    return (h/2) * (a+b);
}

double int_open_degree1(double a,double b,double h){
    return ((3*h)/2) * (a+b);
}

double int_closed_degree2(double a, double b, double c, double h){
    return (h/3) * (a + 4*b + c);
}  

double int_open_degree2(double a,double b,double c,double h){
    return (4/3)*h*(2*a - b + 2*c);
}
                    
double int_closed_degree3(double a, double b, double c, double d,double h){
    return (3/8)*h * (a + 3*b + 3*c + 1*d);
}

double int_open_degree3(double a,double b, double c, double d, double h){
    return (5/24)*h * (11*a + b + c + 11*d);
}

double int_closed_degree4(double a, double b, double c, double d, double e, double h){
    return ((2*h)/45) * (7*a + 32*b + 12*c + 32*d + 7*e);
}

double int_open_degree4(double a,double b, double c, double d,double e, double h){
    return (3/5)*h*( (11/2)*a - 7*b + 13*c - 7*d + (11/2)*e);
}

double f(double x){
    double a = ((sin(2*x)) + (4*pow(x,2)) + (3*x));
    return pow(a,2);
}

double integrate(short int degree,short int closed,int a,int b){
    /*
    Arguments:
        degree - substitution polynomial degree 
        closed - boolean value, 1 = closed e 0 = open
        a e b - integration limits
    */
    
   if (degree > 4 || degree <1){
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
    int n=0;
    
    while(difference>tolerance){
        l = 0;
        aux = result;
        result = 0;
        
        if (closed == 1)
             h = ((b-a)/pow(2,n)) / degree;
        else
             h = ((b-a)/pow(2,n)) / (degree +2);
        

        while(l<pow(2,n)){
            xi = (a+(l*(b-a))/pow(2,n));
            
            if (closed == 1){
                if (degree==4)
                    result += int_closed_degree4(f(xi),f(xi+h),f(xi+2*h),f(xi+3*h),f(xi+4*h),h);
                else if (degree==3)
                    result += int_closed_degree3(f(xi),f(xi+h),f(xi+2*h),f(xi+3*h),h);
                else if (degree==2)
                    result += int_closed_degree2(f(xi),f(xi+h),f(xi+2*h),h);
                else
                    result += int_closed_degree1(f(xi),f(xi+h),h);
            }else{
                if (degree==4)
                    result += int_open_degree4(f(xi+h),f(xi+2*h),f(xi+3*h),f(xi+4*h),f(xi+5*h),h);
                else if (degree==3)
                    result += int_open_degree3(f(xi+h),f(xi+2*h),f(xi+3*h),f(xi+4*h),h);
                else if (degree==2)
                    result += int_open_degree2(f(xi+h),f(xi+2*h),f(xi+3*h),h);
                else
                    result += int_open_degree1(f(xi+h),f(xi+2*h),h);
            }    
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
    short int closed = 0;
    int a = 0;
    int b = 0;


    printf("Enter the substitution polynomial degree(1,2,3 or 4)\n");
    scanf("%hd", &degree);

    printf("\nEnter the approach\n[0]:Open\n[1]:Closed\n");
    scanf("%hd", &closed);

    printf("\nEnter the a\n");
    scanf("%d", &a);

    printf("\nEnter the b\n");
    scanf("%d", &b);

    double result = integrate(degree,closed,a,b);
    printf("Result = %f \n",result);

    return 0;
}