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
            b[i] += matrix[i*n +j] * vetor[j];
    
    return b;
}

//Multiplicacao de matrizes
void matMult(double* matrix0,double* matrix1,double* b,int l0, int c0,int l1,int c1){
    //Linha e coluna da matriz 0 e da matriz1 
    if (c0 != l1){ printf("Colunas da matriz0 != Linhas da matriz1"); exit(0); }

    for(int i = 0; i < l0; i++)
        for(int j = 0; j < c1; j++)
            for(int k = 0; k < c0; k++)
                b[i*c1 +j] += matrix0[i*c0 + k] * matrix1[k*c1 +j];

          
   // return b;
}

//Norma 2 de vetor
double norm(double* vector,int n){
    //vetor, tamanho
    double norm_ = vecVec(vector,vector,n,n);    
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

//Maior elemento de um vetor(ou matriz)
double max(double* arr, int n, int *line) {  
    //vetor, tamanho, linha do maior
    int max_ = fabs(arr[0]); 
    int line_ = 0;
    for (int i = 1; i < n; i++) 
        if (fabs(arr[i]) > fabs(max_)) {
            max_ = fabs(arr[i]); 
            line_ = i;
        }
    *line = line_;
    return max_; 
} 

//Retorna uma coluna duma matriz
double* get_column(double *A, int l,int k){
    //Matriz,num. de linhas,qual coluna
    double* column = aloc(l);
    for (int i = 0; i < l; i++){
        column[i] = A[i*l + k];
    }
    
    return column;
}

//Troca linhas duma matriz(inplace)
void switch_lines(double *A, int m, int k, int l){
    //Matriz, linhas m e k, tamanho da linha

 //   printf("Antes\n");
 //   printMf(A,l*l,l);
    double* aux = aloc(l);
    double* lM = aloc(l);
    double* lK = aloc(l);

    for (int j = 0; j < l; j++)
        aux[j] = A[m*l + j];            //aux = a
        
    for (int j = 0; j < l; j++){
        A[m*l + j] = A[k*l + j] ;  // a = b
        A[k*l + j] = aux[j]    ;   // b = aux
    }
    
 //   printf("Depois\n");
//    printMf(A,l*l,l);

}

//Eliminacao de Gauss com pivoteamento parcial. Guarda L e U na A e guarda as mudancas de linhas
double* LU(double* A,int l){
    //Matriz, num. de linhas(ou colunas)
    double* trocas = aloc(l);
    int* m =  (int*) calloc(1,sizeof(int));

    short int flag = 0;
    
    double pivo = 0;
    //do livro comeca no 1
    for (int k = 0; k < l-1; k++){
        
        pivo = max(get_column(A,l,k),l,m);
        if (pivo == 0){
            flag = 1;
            trocas[k] = 0;
            printf("A é singular");
        }
        
        else{
            trocas[k] = *m;
     
            if (*m != k){
                printf("Trocou %d com %d\n",*m,k);
                switch_lines(A,*m,k,l);
            }
            
            for (int i = k; i < l; i++)
                A[i*l + k] = A[i*l + k] / A[k*l + k];
             for (int j = k; j < l; j++)                // li <- li - m*lj
                for (int i = k; i < l; i++)
                    A[i*l+j] = A[i*l+j] - (A[i*l+k] * A[k*l+j]);                
        }
        
    }
    if (A[((l-1)*l)+(l-1)] == 0){
        flag = 1;
        trocas[l] = 0;
    }
    else{
        trocas[l] = l;
    }

 /*   printf("trocas\n");
    printMf(trocas,l,1); */

    return trocas;
}


/*
//Método da potencia para calculo de autovalores e autovetores
void potenciaInverso(double* A,int l, double* valor,double *vetor){
    //Matriz,num de linhas(ou col), variaveis que vão ser armazenados os auto
    
    double* v = aloc(l);
    double* aux = aloc(l);
 //   printf("Endereço do vetor que chegou = %p \n",vetor);

    for (int i = 0; i < l; i++) v[i] = 1;
    
    double tolerancia = 1.E-10;
    double lambda = 0, lambda_velho = 0, dif = 1;
    
   // A = LU(A)                  //Parte superior de A temm o U, inferior tem o L 

    while (dif > tolerancia){
        lambda_velho = lambda;
      //  vetor = vectorDiv(v,norm(v,l),l);
        aux = vectorDiv(v,norm(v,l),l);
        for (int i = 0; i < l; i++)                 //Pra n mudar o endereço do vetor
            vetor[i] = aux[i];
    //    printf("Endereço do vetor no laco = %p \n",vetor);
    //    v = luSolver(A,vetor);
        lambda = vecVec(vetor,v,l,l);
        dif = fabs((lambda - lambda_velho)/lambda);
    }

    *valor = lambda;                                //Value no endereço de valor agora é lambda
}
/*
def potencia_inverso(A):
    v = np.ones((A.shape[0],1))
    tolerancia = 1E-10
    
    LU,piv = sci.lu_factor(A)               #LU tem o U na parte superior e o L na inferior
    
    lambda_ = 0
    dif = 1

    while dif > tolerancia:
        lambda_velho=lambda_
        x1 = v / np.linalg.norm(v)
        v = sci.lu_solve((LU, piv), x1)
        lambda_ = x1.T.dot(v)
        dif = abs((lambda_-lambda_velho)/lambda_)

        #print('autovetor = \n ',x1)
        #print('autovalor = ',lambda_,'\n')
        
    lambda_ = 1/lambda_
    return x1,lambda_

def potencia_deslocamento(A,desloc):
    A = A - (desloc*np.eye(A.shape[0]))

    vetor,valor = potencia_inverso(A)
    valor += desloc
    
    return vetor,valor
*/

int main (){
    double matriz1[9] = {5.,2.,1.,2.,3.,1.,1.,1.,2.};
    double matriz2[9] = {-14.,1.,-2.,1.,-1.,1.,-2.,1.,-11.};
    double matriz3[25] = {40.,8.,4.,2.,1.,8.,30.,12.,6.,2.,4.,12.,20.,1.,2.,2.,6.,1.,25.,4.,1.,2.,2.,4.,5.};
    int l = (int) sizeof(matriz3)/sizeof(matriz3[0]);
    int n = sqrt(l);

    double* auto_vetor = aloc(n);
    double* auto_valor = aloc(1);
    printf("Antes\n");
    printMf(matriz3,l,n);
    
    double *trocas = aloc(n);
    trocas = LU(matriz3,n);
    printf("Depois\n");
    
    printMf(matriz3,l,n);

//    switch_lines(matriz3,0,4,n);
//    printMf(matriz3,l,5);
//    printMf(get_column(matriz3,5,2),5,1);
 //   printf("%f",max(get_column(matriz3,n,2),l));
 /*   potenciaRegular(matriz3,n,auto_valor,auto_vetor);
    vetor = potenciaInverso(matriz3,n,valor);  
    printf("Autovetor \n");
    printMf(auto_vetor,n,1);
    printf("Autovalor \n");
    printf("%f \n",*auto_valor);*/

    return 0;
}