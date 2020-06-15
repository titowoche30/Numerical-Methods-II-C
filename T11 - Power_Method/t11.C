#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <time.h>

//Aloca vetor double com tamanho size
double* aloc(int size){
    double* vector = (double*) calloc(size,sizeof(double));
    return vector;
}

//Printa matriz de double
void printMf(double* vetor,int n,int esp){
    //Vetor, num de elementos, numero de elementos numa linha(num de colunas)
    for(int i = 1; i <= n; i++){
        printf("%f ",vetor[i-1]);
        if(i%esp == 0) printf("\n");
    }
    printf("\n");
}

double vecVec(double* vector1, double* vector2,int m,int n){
    //m = linhas do 1,n = linhas do 2

    if (m != n){ printf("Dimensões incompatíveis dos vetores"); exit(0); }

    double result = 0;
    for(int i = 0; i < m; i++)
        result += vector1[i] * vector2[i];
    
    return result;
}

//Multiplicacao matriz com vetor
double* matVec(double* matrix,double* vetor,int m,int n,int p){
    //matriz,vetor,linhas da matriz,colunas da matriz,linhas do vetor
    if (n != p){ printf("Colunas da matriz != Linhas do vetor"); exit(0); }

    double* b = aloc(m);
    for(int i = 0; i < m; i++)
        for(int j = 0; j < n; j++)
            b[i] = b[i] + matrix[i*n +j] * vetor[j];
    
    return b;
}

double norm(double* vector,int n){
    double norm_ = 0;

    for(int i = 0; i < n; i++)
        norm_+=vector[i]*vector[i];
    
    return sqrt(norm_);
}

double* vectorDiv(double* vector,double scaler,int l){
    double* new_vector = aloc(l);

    for (int i = 0; i < l; i++)
        new_vector[i] = vector[i] / scaler;
    
    return new_vector;
}

double* potenciaRegular(double* A,int l, double* valor){
    double* v = aloc(l);
    double* vetor = aloc(l);

    for (int i = 0; i < l; i++)
        v[i] = 1;
    
    double tolerancia = 1.E-10;
    double lambda = 0, lambda_velho = 0, dif = 1;
    

    while (dif > tolerancia){
        lambda_velho = lambda;
        vetor = vectorDiv(v,norm(v,l),l);
        v = matVec(A,vetor,l,l,l);
        lambda = vecVec(vetor,v,l,l);
        dif = fabs((lambda - lambda_velho)/lambda);
    }

    valor[0] = lambda;
    return vetor;
}

int main (){
    double matriz1[9] = {5.,2.,1.,2.,3.,1.,1.,1.,2.};
    double matriz2[9] = {-14.,1.,-2.,1.,-1.,1.,-2.,1.,-11.};
    double matriz3[25] = {40.,8.,4.,2.,1.,8.,30.,12.,6.,2.,4.,12.,20.,1.,2.,2.,6.,1.,25.,4.,1.,2.,2.,4.,5.};

    int l = 25;
    int n = sqrt(l);

    double* vetor = aloc(n);
    double* valor = aloc(1);
    vetor = potencia_regular(matriz3,n,valor);
    printf("Autovetor \n");
    printMf(vetor,n,1);
    printf("Autovalor \n");
    printf("%f \n",valor[0]);

    return 0;
}