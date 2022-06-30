/*************************************************/
// this program is a part of GMM_DIFRACT_TOOLS. 
// Written by T. Nagai 
/*************************************************/
#include<stdio.h>
#include<math.h>
#include <complex.h>
#define DIM 3


double get_abs_struct_factor(double kx, double ky, double kz,double pi[], double mu[][DIM], double sigma[][DIM][DIM],int nmu){
	int i,j,k;
	//printf("%f, %f\n", kx,ky);
	//for(i=0; i<nmu;i++){
	//	printf("for muid=%d\n",i);
	//	printf("%f  %f  %f\n",mu[i][0],mu[i][1],mu[i][2]);
	//	printf("\n");
	//}
	//for(i=0; i<nmu;i++){
	//	printf("sigma for muid=%d\n",i);
	//	for(j=0; j<DIM;j++){
	//		printf("%f  %f  %f\n",sigma[i][j][0],sigma[i][j][1],sigma[i][j][2]);
	//	}
	//	printf("\n");
	//}
	double S[DIM];
	S[0]=2.0*M_PI*kx;
	S[1]=2.0*M_PI*ky;
	S[2]=2.0*M_PI*kz;

	double _Complex f=0; // Complex Intensity
	double s_mu=0;
	double s_sig_s=0;
	
	for(i=0; i<nmu;i++){
		s_mu=0;
		s_sig_s=0;
		for(j=0;j<DIM;j++){
			s_mu+=S[j]*mu[i][j];
		}
		for(j=0;j<DIM;j++){
			for(k=0;k<DIM;k++){
				s_sig_s+=S[j]*sigma[i][j][k]*S[k];
			}
		}
		f+=pi[i]*cexp(I*s_mu)*exp(-s_sig_s*0.5);
		//printf("f's real=%f  imag=%f\n", creal(f),cimag(f));
	}
	double F=cabs(f);
	if(F<pow(10,-320)){F=pow(10,-320);}

	return F;
}


#include<omp.h>
void get_abs_struct_factor_all_grid_point(double* diff, int num_grid, double det_width, double det_dist, double wave_length, double pi[], double mu[][DIM], double sigma[][DIM][DIM],int nmu){
	int i,j,k;
	int x_index, y_index;
	double _Complex f=0; // Complex Intensity
	double F;            // for absolute Intensity
	double s_mu;
	double s_sig_s;
	double S[3];
	double dx,dy;
	double kx,ky,kz;
	double k_inc=1.0/wave_length;
        //printf("I am using this\n");

	#pragma omp parallel for schedule(guided) private(kx,ky,kz,dx,dy,i,j,k,y_index, f, F,s_mu, s_sig_s,S) 
	for(x_index=0; x_index<num_grid; x_index++){
		for(y_index=0; y_index<num_grid; y_index++){
			dx=2.0*det_width/(double)(num_grid-1)*x_index-det_width;
			dy=2.0*det_width/(double)(num_grid-1)*y_index-det_width;
			kx=(k_inc*dx/sqrt(det_dist*det_dist+dx*dx+dy*dy));
			ky=(k_inc*dy/sqrt(det_dist*det_dist+dx*dx+dy*dy));
			kz=(k_inc-sqrt(k_inc*k_inc-kx*kx-ky*ky));
			//fprintf(stderr,"debug: %d %d %10.7f %10.7f %10.7f %10.7f %10.7f ; %10.8f\n",x_index,y_index,dx/1e10,dy/1e10,kx,ky,kz,(det_dist*det_dist+dx*dx+dy*dy));
			S[0]=2.0*M_PI*kx;
			S[1]=2.0*M_PI*ky;
			S[2]=2.0*M_PI*kz;
			f=0; 
			for(i=0; i<nmu;i++){
				s_mu=0;s_sig_s=0;
				for(j=0;j<DIM;j++){
					s_mu+=S[j]*mu[i][j];
				}
				for(j=0;j<DIM;j++){
					for(k=0;k<DIM;k++){
						s_sig_s+=S[j]*sigma[i][j][k]*S[k];
					}
				}
				f+=pi[i]*cexp(I*s_mu)*exp(-s_sig_s*0.5);
			}
			F=cabs(f);
			if(F<pow(10,-320)){F=pow(10,-320);}
			diff[x_index*num_grid+y_index]=F;
		}
	}
}


