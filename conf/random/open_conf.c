#include <stdio.h>
#include <stdlib.h>

struct Lattice {
    int Nx, Ny, Nz, Nt;
};


int main() {
    int Nx, Ny, Nz, Nt;
    double plaq;
    FILE* ptr_conf;
    double SU3[18];
    char File[500];
    char line[1024];
    sprintf(File, "4x4x4x4_random");

    ptr_conf = fopen(File, "r");
    if (!ptr_conf) {
        printf("Unable to open file!");
        return 1;
    }
    struct Lattice vol;
    
	for (int i = 0; i <= 100; i++){
        if (i == 0){fread(&vol, sizeof(struct Lattice), 1, ptr_conf); 
        printf("Nx,Ny,Nz,Nt = %d,%d,%d,%d\n", vol.Nx, vol.Ny, vol.Nz, vol.Nt);}
        else if (i==1){fread(&plaq, sizeof(double), 1, ptr_conf);}
        else{
		    fread(&SU3, sizeof(double), 18, ptr_conf);
            for (int j = 0; j < 18; j++) {
                printf("%-30.17g ", SU3[j]);
            }
            printf("\n");
        }
	}
	fclose(ptr_conf);

    return 0;
}