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

//Produto interno de 2 vetores
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

//Norma-2 de vetor
double norm(double* vector,int n){
    //vetor, tamanho

    double norm_ = 0;

    for(int i = 0; i < n; i++)
        norm_+=vector[i]*vector[i];
    
    return sqrt(norm_);
}

//Divisao de vetor(ou matriz) por escalar
double* vectorDiv(double* vector,double scaler,int l){
     //vetor,escalar,tamanho
    double* new_vector = aloc(l);

    for (int i = 0; i < l; i++)
        new_vector[i] = vector[i] / scaler;
    
    return new_vector;
}

//Método da potencia para calculo de autovalores e autovetores
void potenciaRegular(double* A,int l, double* valor,double *vetor){
    //Matriz,num de linhas(ou col), variaveis que vão ser armazenados os auto
    
    double* v = aloc(l);
    double* aux = aloc(l);
 //   printf("Endereço do vetor que chegou = %p \n",vetor);

    for (int i = 0; i < l; i++) v[i] = 1;
    
    double tolerancia = 1.E-10;
    double lambda = 0, lambda_velho = 0, dif = 1;
    
    while (dif > tolerancia){
        lambda_velho = lambda;
      //  vetor = vectorDiv(v,norm(v,l),l);
        aux = vectorDiv(v,norm(v,l),l);
        for (int i = 0; i < l; i++)                 //Pra n mudar o endereço do vetor
            vetor[i] = aux[i];
    //    printf("Endereço do vetor no laco = %p \n",vetor);
        v = matVec(A,vetor,l,l,l);
        lambda = vecVec(vetor,v,l,l);
        dif = fabs((lambda - lambda_velho)/lambda);
    }

    *valor = lambda;                                //Value no endereço de valor agora é lambda
}

int main (){
    double matriz1[9] = {5.,2.,1.,2.,3.,1.,1.,1.,2.};
    double matriz2[9] = {-14.,1.,-2.,1.,-1.,1.,-2.,1.,-11.};
    double matriz3[25] = {40.,8.,4.,2.,1.,8.,30.,12.,6.,2.,4.,12.,20.,1.,2.,2.,6.,1.,25.,4.,1.,2.,2.,4.,5.};
    int l = (int) sizeof(matriz3)/sizeof(matriz3[0]);
    int n = sqrt(l);

    double* auto_vetor = aloc(n);
    double* auto_valor = aloc(1);

    potenciaRegular(matriz3,n,auto_valor,auto_vetor);
    printf("Autovetor \n");
    printMf(auto_vetor,n,1);
    printf("Autovalor \n");
    printf("%f \n",*auto_valor);

    return 0;
}