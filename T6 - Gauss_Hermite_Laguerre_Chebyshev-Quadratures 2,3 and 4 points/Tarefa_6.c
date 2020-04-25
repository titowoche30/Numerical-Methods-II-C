#include <stdio.h>
#include <stdlib.h> 
#include <math.h>

#define pi 3.141592653589793238462643383279502884197169399375105820974

double f_her(double x){
    return pow((x+1),2);
    //return cos(x);
    //return pow(x,2);
}

double f_lag(double x){
    //return pow((x-1),3);
    //return pow(x,3);
    return x*(pow(x,2)+4);
}

double f_cheb(double x){
    //return x+1
    //return  (2*x) / (pow((pow(x,2)+1),2));
    return x/ (pow((pow(x,2)+4),3.0/2.0));
}

double her(int n){
    /*
    Gauss-Hermite
    For functions of type: exp(-x^2)* f(x)
    a = -inf , b = +inf

    Argumentos
        f - funcao lambda ou normal
        n - grau do polinomio de Hermite
    */
   
    double x[n];
    double w[n];
    double result = 0;

    for (int i = 0; i < n; i++){
        x[i] = 0;
        w[i] = 0;
    }
    

    if (n==2){
        x[0] = -1/sqrt(2.0);
        x[1] = -x[0];
        w[0] = sqrt(pi)/2.0;
        w[1] = w[0];
    }
    else if (n==3){
        x[0] = -(sqrt(3.0/2.0));
        x[1] = 0;
        x[2] = -x[0];
        w[0] = sqrt(pi)/6.0;
        w[1] = (2*sqrt(pi))/3.0;
        w[2] = w[0];
    }
    else{
        x[0] = -(sqrt(3.0/2.0 + sqrt(3.0/2.0)));
        x[1] = -(sqrt(3.0/2.0 - sqrt(3.0/2.0)));
        x[2] = -x[1];
        x[3] = -x[0];
        w[0] = sqrt(pi)/(4*(3+sqrt(6.0)));
        w[1] = -sqrt(pi)/(4*(sqrt(6.0)-3));
        w[2] = w[1];
        w[3] = w[0];
    }

    for (int i = 0; i < n; i++){
        result+=f_her(x[i])*w[i];
    }
    
    return result;
}
        
double lag(int n){
    /*
    Gauss-Laguerre
    For functions of type: exp(-x)* f(x)
    a = 0 , b = +inf

    Argumentos
        f - funcao lambda ou normal
        n - grau do polinomio de Laguerre
    */
    
    double x[n];
    double w[n];
    double result = 0;

    for (int i = 0; i < n; i++){
        x[i] = 0;
        w[i] = 0;
    }

    if (n==2){
        x[0] = 2-sqrt(2.0);
        x[1] = 2+sqrt(2.0);
        w[0] = (1.0/4.0)*(2+sqrt(2.0));
        w[1] = (1.0/4.0)*(2-sqrt(2.0));
    }
    else if (n==3){
        x[0] = 0.4157745568;
        x[1] = 2.2942803603;
        x[2] = 6.2899450829;
        w[0] = 0.7110930099;
        w[1] = 0.2785177336;
        w[2]=  0.0103892565;
    }
    else{
        x[0] = 0.32255;
        x[1] = 1.7458;
        x[2] = 4.5366;
        x[3] = 9.3951;
        w[0] = 0.603115;
        w[1] = 0.357347;
        w[2] = 0.0388894;
        w[3] = 0.000539278;
    }

    for (int i = 0; i < n; i++){
        result+=f_lag(x[i])*w[i];
    }
    
    return result;   
}

double cheb(int n){
    /*
    Gauss-Chebyshev
    For functions of type: 1/sqrt(1-x^2)* f(x)

    a = -1 , b = +1

    Argumentos
        f - funcao lambda ou normal
        n - grau do polinomio de Chebyshev
    */
    
    double x[n];
    double w[n];
    double result = 0;

    for (int i = 0; i < n; i++){
        x[i] = 0;
        w[i] = 0;
    }

    if (n==2){
        x[0] = -1/sqrt(2.0);
        x[1] = -x[0];
        w[0] = pi/2.0;
        w[1] = w[0];
    }
    else if (n==3){
        x[0] = -(sqrt(3.0)/2.0);
        x[1] = 0;
        x[2] = -x[0];
        w[0] = pi/3.0;
        w[1] = w[0];
        w[2]=  w[0];
    }
    else{
        x[0] = -(sqrt(2.0 + sqrt(2.0))/2.0);
        x[1] = -(sqrt(2.0 - sqrt(2.0))/2.0);
        x[2] = -x[1];
        x[3] = -x[0];
        w[0] = pi/4.0;
        w[1] = w[0];
        w[2] = w[0];
        w[3] = w[0];
    }

    for (int i = 0; i < n; i++){
        result+=f_cheb(x[i])*w[i];
    }
    
    return result; 
}

int main (){
    short int method = 0;
    int a = 0;
    int b = 0;
    double resultado = 0;

    printf("[1] - Hermite (-inf,+inf) \n[2] - Laguerre (0,+inf) \n[3] - Chebyshev (-1,1)\n");
    scanf("%hd", &method);

    if (method > 3 || method < 1){
        printf("Invalid Digit\n");
        return 0;
    }

    if (method==1)
        for (int i = 1; i <= 3; i++){
            printf("Grau %d\n",i);
            resultado=her(i+1);
            printf("Resultado= %f\n",resultado);
        }

    else if (method==2)
        for (int i = 1; i <= 3; i++){
            printf("Grau %d\n",i);
            resultado=lag(i+1);
            printf("Resultado= %f\n",resultado);
        }
    else
       for (int i = 1; i <= 3; i++){
            printf("Grau %d\n",i);
            resultado=cheb(i+1);
            printf("Resultado= %f\n",resultado);
        }        
                    
 
    return 0;
}