#include<stdio.h>
#include <stdlib.h>

// Dirac Matrix /
struct DMatrix{
	int row, column;
	double real, imag;
};

int main(int argc, char** argv) {
	FILE* ptr_Dmatrix;
	struct DMatrix D;
	int Nx, Ny, Nz, Nt;
	Nx = atoi(argv[1]);
	Ny = atoi(argv[2]);
	Nz = atoi(argv[3]);
	Nt = atoi(argv[4]);
	char File[100];
	sprintf( File, "%dx%dx%dx%d_DiracMatrix", Nt,Nz,Ny,Nx);

	ptr_Dmatrix = fopen(File, "rb");
	if (!ptr_Dmatrix) {
		printf("Unable to open file!");
		return 1;
	}

	//Print the first 100 non-zero entries//
	printf("%-5s%-15s%-30s%-30s\n", "Row", "Column", "Real part", "Imaginary part");
	for (int i = 1; i <= 100; i++){
		fread(&D, sizeof(struct DMatrix), 1, ptr_Dmatrix);
		printf("%-5d%-15d%-30.15g%-30.15g\n", D.row, D.column, D.real, D.imag);
	}
	fclose(ptr_Dmatrix);

	return 0;
}