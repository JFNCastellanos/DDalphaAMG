#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <complex.h>
#include <time.h>
#include <math.h>
#include <assert.h>

typedef double complex complex_double;

#define Nx 4
#define Ny 4
#define Nz 4
#define Nt 4
#define N 12*Nx*Ny*Nz*Nt
#define S 4
#define C 3

complex_double gammaMat[4][4][4];  //Gamma matrices
double conf[4*Nx*Ny*Nz*Nt][18]; //SU(3) links

int dim[4] = {Nx, Ny, Nz, Nt};
int IndX[N], IndY[N], IndZ[N], IndT[N], IndS[N], IndC[N];
double m0=1, h=1; //Default values for mass and lattice spacing
int n_mu[4], n_pmu[4];
int Identity[4][4] = {{1, 0, 0, 0},
                         {0, 1, 0, 0},
                         {0, 0, 1, 0},
                         {0, 0, 0, 1}};

//Gamma matrices definition
void gammas_init(){
    //gamma0
    gammaMat[0][0][0] = 0.0; gammaMat[0][0][1] = 0.0; gammaMat[0][0][2] = 0.0; gammaMat[0][0][3] = 1.0*_Complex_I;
    gammaMat[0][1][0] = 0.0; gammaMat[0][1][1] = 0.0; gammaMat[0][1][2] = 1.0*_Complex_I; gammaMat[0][1][3] = 0.0;
    gammaMat[0][2][0] = 0.0; gammaMat[0][2][1] = -1.0*_Complex_I; gammaMat[0][2][2] = 0.0; gammaMat[0][2][3] = 0.0;
    gammaMat[0][3][0] = -1.0*_Complex_I; gammaMat[0][3][1] = 0.0; gammaMat[0][3][2] = 0.0; gammaMat[0][3][3] = 0.0;

    //gamma1
    gammaMat[1][0][0] = 0.0; gammaMat[1][0][1] = 0.0; gammaMat[1][0][2] = 0.0; gammaMat[1][0][3] = -1.0;
    gammaMat[1][1][0] = 0.0; gammaMat[1][1][1] = 0.0; gammaMat[1][1][2] = 1.0; gammaMat[1][1][3] = 0.0;
    gammaMat[1][2][0] = 0.0; gammaMat[1][2][1] = 1.0; gammaMat[1][2][2] = 0.0; gammaMat[1][2][3] = 0.0;
    gammaMat[1][3][0] = -1.0; gammaMat[1][3][1] = 0.0; gammaMat[1][3][2] = 0.0; gammaMat[1][3][3] = 0.0;

    //gamma2
    gammaMat[2][0][0] = 0.0; gammaMat[2][0][1] = 0.0; gammaMat[2][0][2] = 1.0*_Complex_I; gammaMat[2][0][3] = 0.0;
    gammaMat[2][1][0] = 0.0; gammaMat[2][1][1] = 0.0; gammaMat[2][1][2] = 0.0; gammaMat[2][1][3] = -1.0*_Complex_I;
    gammaMat[2][2][0] = -1.0*_Complex_I; gammaMat[2][2][1] = 0.0; gammaMat[2][2][2] = 0.0; gammaMat[2][2][3] = 0.0;
    gammaMat[2][3][0] = 0.0; gammaMat[2][3][1] = 1.0*_Complex_I; gammaMat[2][3][2] = 0.0; gammaMat[2][3][3] = 0.0;

    //gamma3
    gammaMat[3][0][0] = 0.0; gammaMat[3][0][1] = 0.0; gammaMat[3][0][2] = 1.0; gammaMat[3][0][3] = 0.0;
    gammaMat[3][1][0] = 0.0; gammaMat[3][1][1] = 0.0; gammaMat[3][1][2] = 0.0; gammaMat[3][1][3] = 1.0;
    gammaMat[3][2][0] = 1.0; gammaMat[3][2][1] = 0.0; gammaMat[3][2][2] = 0.0; gammaMat[3][2][3] = 0.0;
    gammaMat[3][3][0] = 0.0; gammaMat[3][3][1] = 1.0; gammaMat[3][3][2] = 0.0; gammaMat[3][3][3] = 0.0;
}

//Indices definition
void indices_init(){
    int count = 0;
    for(int t = 0; t < Nt; t++){
        for(int z = 0; z < Nz; z++){
            for(int y = 0; y < Ny; y++){
                for(int x = 0; x < Nx; x++){
                    for(int s=0; s < S; s++){
                        for(int c=0; c < C; c++){  
                            IndX[count] = x; IndY[count] = y; IndZ[count] = z; IndT[count] = t; 
                            IndS[count] = s; IndC[count] = c;
                            count += 1;
                        }
                    }
                }
            }
        }
    }
}

//Modulus function
static inline int mod(int a, int b){
    int r = a % b;
    return r < 0 ? r + b : r;
}

//Kronecker delta for ints
static inline int kDelta(int i, int j){
    if(i == j) return 1;
    else return 0;
}

//Kronecker delta for arrays
static inline int kDeltaArr(int *i, int *j){
    for(int k = 0; k < 4; k++){
        if(i[k] != j[k]) return 0;
    }
    return 1;
}

//Returns a number for the coordinate (t,z,y,x) in the mu-direction
static inline int Coord(int *n, int mu){
    return mu + 4*n[3] + 4*Nx*n[2] + 4*Nx*Ny*n[1] + 4*Nx*Ny*Nz*n[0];
}

void periodic_boundary(int *n_vec, int mu){
    for(int i = 0; i < 4; i++){
        n_mu[i] = mod((n_vec[i]-Identity[mu][i]), dim[i]); //n-hat{mu} with periodic boundary
        n_pmu[i] = mod((n_vec[i]+Identity[mu][i]), dim[i]); //n+hat{mu} with periodic boundary
    }
}

complex_double DiracEntry(int n, int m){
    /*
    *alfa, beta: 0,1,2,3 (spin indices)
    *c, d: 0,1,2 (color indices)
    *n, m: 0,1,...,Nx Ny Nz Nt (volume)
    *m0: bare mass
    *h: lattice spacing
    */
    int n_vec[4]={IndX[n], IndY[n], IndZ[n], IndT[n]};
    int m_vec[4]={IndX[m], IndY[m], IndZ[m], IndT[m]};
    int alfa = IndS[n], beta = IndS[m]; //Spin indices
    int a = IndC[n], b = IndC[m]; //Color indices
    complex_double Dnm = (m0+4)/h * kDelta(alfa,beta) * kDelta(a,b) * kDeltaArr(n_vec, m_vec);
    
    for (int mu=0; mu<4; mu++){
        periodic_boundary(n_vec,mu); //n_mu and n_pmu
        int delta1 = kDeltaArr(n_pmu,m_vec), delta2 = kDeltaArr(n_mu,m_vec);
        if (delta1 != 0 || delta2 != 0){
            //SU3[i,j] = U[2*(i+j*3)] + 1j*U[2*(i+j*3)+1]
            //U_mu(n)[a,b]
            complex_double Unab = conf[Coord(n_vec,mu)][2*(a+b*3)] + _Complex_I*conf[Coord(n_vec,mu)][2*(a+b*3)+1];
            //U*_m(n-hat{mu})[b,a]
            complex_double Un_mu_ba =  conf[Coord(n_mu,mu)][2*(b+a*3)] - _Complex_I*conf[Coord(n_mu,mu)][2*(b+a*3)+1];
            Dnm += -1.0 / (2.0*h) * (1.0*Identity[alfa][beta]-gammaMat[mu][alfa][beta]) * Unab * delta1;
            Dnm += -1.0 / (2.0*h) * (1.0*Identity[alfa][beta]+gammaMat[mu][alfa][beta]) * Un_mu_ba * delta2;
        }
    }
    
    return Dnm;
}


int main(){
    printf("Nx=%d, Ny=%d, Nz=%d, Nt=%d\n", Nx, Ny, Nz, Nt);
    gammas_init();
    indices_init();
    //----------Read configuration----------//
    FILE* ptr_conf;
    char File[100];
    char line[1000];
    sprintf(File, "%dx%dx%dx%d_random.txt",Nt,Nz,Ny,Nx); //Open the binary file
    ptr_conf = fopen(File, "r");
    if (!ptr_conf) {
        printf("Unable to open file!");
        return 1;
    }
    int count = 0;
    while (fgets(line, sizeof(line), ptr_conf) && count < N) {
        char *ptr = line;
        //if (count > 1){
            for (int j = 0; j < 18; j++) {
                conf[count][j] = strtod(ptr, &ptr);
            }
        //}
        count++;
    }
    fclose(ptr_conf);
    //------------------------------------//

    //----------Write Dirac matrix to file----------//
    clock_t start, end;
    start = clock();
    FILE *ftxt = NULL;
    FILE *fout = NULL;
    char s[100];
    sprintf( s, "%dx%dx%dx%d_DiracMatrix.txt", Nt,Nz,Ny,Nx); //This one is to compare with python (remove later)
    assert( ( ftxt = fopen( s, "wb" ) ) != NULL ); 
    sprintf( s, "%dx%dx%dx%d_DiracMatrix", Nt,Nz,Ny,Nx);
    assert( ( fout = fopen( s, "wb" ) ) != NULL ); 

    complex_double Dnm;
    //Loop over the matrix entries
    for(int i =0;i<N;i++){
        for(int j=0; j<N; j++){
            Dnm = DiracEntry(i,j);
            //Print row, column, real part and imaginary part of the non-zero elements    
            if ( sqrt(creal(Dnm)*creal(Dnm) + cimag(Dnm)*cimag(Dnm)) > 1e-12){
                fprintf(ftxt, "%-15d%-15d%-30.17g%-30.17g\n", i,j,creal(Dnm),cimag(Dnm));
                double real_part = creal(Dnm); double imag_part = cimag(Dnm); 
                fwrite(&i, sizeof(int), 1, fout);
                fwrite(&j, sizeof(int), 1, fout);
                fwrite(&real_part, sizeof(double), 1, fout);
                fwrite(&imag_part, sizeof(double), 1, fout);
                fflush(0);
            }
        }
    }
    fclose( ftxt ); 
    fclose( fout );
    end = clock();
    double time_taken = (double)(end - start) / (double)(CLOCKS_PER_SEC);
    printf("Time taken to assemble the Dirac matrix: %g secs\n", time_taken);
    //------------------------------------//
    return 0;

}
    