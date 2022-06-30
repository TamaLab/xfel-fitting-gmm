/*************************************************/
// this program is a part of GMM_DIFRACT_TOOLS. 
// Written by T. Nagai 
/*************************************************/
#include<stdio.h>
#include<math.h>
#define DIM 3

int rotate(double mu[][DIM], double sig[][DIM][DIM],int nmu, double matrix[DIM][DIM]){
	int i,j,k;
	double mu_tmp[nmu][DIM];
	double sigma_tmp[nmu][DIM][DIM];
	for(i=0; i<nmu;i++){
		for(j=0;j<DIM;j++){
			mu_tmp[i][j]=0.0;
		}
	}
	for(i=0; i<nmu;i++){
		for(j=0;j<DIM;j++){
			for(k=0;k<DIM;k++){
				mu_tmp[i][j]+=matrix[j][k]*mu[i][k];
			}
		}
		for(j=0;j<DIM;j++){
			mu[i][j]=mu_tmp[i][j];
		}
	}
	int m;
	for(i=0; i<nmu;i++){
		for(j=0;j<DIM;j++){
			for(k=0;k<DIM;k++){
				sigma_tmp[i][j][k]=0.0;
			}
		}
		for(j=0;j<DIM;j++){
			for(k=0;k<DIM;k++){
				for(m=0;m<DIM;m++){
					sigma_tmp[i][j][k]+=sig[i][j][m]*matrix[k][m];
				}
			}
		}
		for(j=0;j<DIM;j++){
			for(k=0;k<DIM;k++){
				sig[i][j][k]=sigma_tmp[i][j][k];
				sigma_tmp[i][j][k]=0.0;
			}
		}
		for(j=0;j<DIM;j++){
			for(k=0;k<DIM;k++){
				for(m=0;m<DIM;m++){
					sigma_tmp[i][j][k]+=matrix[j][m]*sig[i][m][k];
				}
			}
		}
		for(j=0;j<DIM;j++){
			for(k=0;k<DIM;k++){
				sig[i][j][k]=sigma_tmp[i][j][k];
			}
		}
	}
	return 0;
}


//
//void main(){
//	int nmu=2;
//	double pi[2]={0.5,0.5};
//	double mu[2][DIM]={{0.0,-10.1,0.0},{0.0, 10.1,0.0}};
//	double sigma[2][DIM][DIM]={{{3.3,0.2,0.},{0.2,0.1,0.},{0.0,0.0,0.1}},{{3.3,0.0,0.},{0.0,0.1,0.3},{0.0,0.3,0.1}}};
//
//	//double matrix[DIM][DIM]={{0.7071,-0.7071,0.},{0.7071,0.7071,0.},{0.0,0.0,1.0}};
//	double matrix[DIM][DIM]={{0.0,-1.0,0.},{1.00,0.00,0.},{0.0,0.0,1.0}};
//
//
//	//printf("%f %f %f\n", mu[0][0], mu[0][1], mu[0][2]);
//	//printf("%f %f %f\n", mu[1][0], mu[1][1], mu[1][2]);
//	int i;
//	printf("====\n");
//	for(i=0;i<DIM;i++){
//		printf("%f %f %f\n", sigma[0][i][0], sigma[0][i][1],sigma[0][i][2]);
//	}
//	printf("====\n");
//	for(i=0;i<DIM;i++){
//		printf("%f %f %f\n", sigma[1][i][0], sigma[1][i][1],sigma[1][i][2]);
//	}
//	rotate(pi, mu, sigma,nmu, matrix);
//
//	for(i=0;i<DIM;i++){
//		printf("%f %f %f\n", matrix[i][0], matrix[i][1], matrix[i][2]);
//	}
//	//printf("%f %f %f\n", mu[0][0], mu[0][1], mu[0][2]);
//	//printf("%f %f %f\n", mu[1][0], mu[1][1], mu[1][2]);
//	printf("====\n");
//	for(i=0;i<DIM;i++){
//		printf("%f %f %f\n", sigma[0][i][0], sigma[0][i][1],sigma[0][i][2]);
//	}
//	printf("====\n");
//	for(i=0;i<DIM;i++){
//		printf("%f %f %f\n", sigma[1][i][0], sigma[1][i][1],sigma[1][i][2]);
//	}
//}
//
//
