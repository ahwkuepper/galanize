/***************************************************************************
 *   Copyright (C) 2010-2016 by Andreas Kuepper                            *
 *   ahwkuepper@gmail.com                                                  *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
/***********************
* particle[][0] = mass *
* particle[][1] = x    *
* particle[][2] = y    *
* particle[][3] = z    *
* particle[][4] = vx   *
* particle[][5] = vy   *
* particle[][6] = vz   *
* particle[][7] = phif *
* particle[][8] = mgc  *
* particle[][9] = rh   *
************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>
#include<sys/stat.h>
#include "galanize.h"


int main (int argc, const char * argv[]) {
	
	int i,j;
	
	clock_t t1, t2;                       //start stop watch
	t1 = clock(); 
	
	if (argc>1) {
		
		//get scale mass and radius from command line
		double Mscale, Rscale, tscale, vscale;
		if (argc >= 3) {
			Mscale = atof(argv[2]);  //in Msun
			Rscale = atof(argv[3]);  //in pc
			tscale = 14.94*pow(pow(Rscale,3)/Mscale,0.5);   //in Myr 
			vscale = Rscale/tscale;  //in pc/Myr ~ km/s
		} else {
			Mscale = 1.0;  //N-body units
			Rscale = 1.0;
			tscale = 1.0;
			vscale = 1.0;
		}		
		
		//open file
		FILE *dat;
		char *fname=calloc(100,sizeof(char));
		strcpy (fname,argv[1]); 
		strcat (fname,".POS");
		dat=fopen(fname,"rb");
		if (dat==NULL) {
			printf("\nFile %s not found!\n\n", fname);
			return 0;
		} 
		
		printf("\n##########################################################\n");
		printf("Analyzing model %s \n", fname);		
		
		
		long lSize, lPos;
		fseek (dat , 0 , SEEK_END);      // obtain file size
		lSize = ftell (dat);
		rewind (dat);
		
		if (lSize) {
			
			printf("File size: %li Bytes\n",lSize);
			printf("##########################################################\n");
			
			
			//main data array
			int columns = 10;
			double **particle;
			particle = (double **)calloc(NMAX,sizeof(double *));
			for (j=0;j<NMAX;j++){
				particle[j] = (double *)calloc(columns,sizeof(double));
				if (particle[j] == NULL) {
					printf("\nMemory allocation failed!\n\n");
					return 0;
				}
			}			
			
			//temporary array
			double *temp;
			temp = (double *)calloc(NMAX,sizeof(double));
			
			double time;
			int nbody;
			
			do {
				
				fread (&time,sizeof(double),1,dat);//time
				time *= tscale;
				
				fread (&nbody,sizeof(int),1,dat);//N
				
				fread (temp,sizeof(double),nbody,dat);//mass
				for (i=0;i<nbody;i++) particle[i][0] = Mscale*temp[i];
				
				fread (temp,sizeof(double),nbody,dat);//mgc
				for (i=0;i<nbody;i++) particle[i][8] = temp[i];
				
				fread (temp,sizeof(double),nbody,dat);//rh
				for (i=0;i<nbody;i++) particle[i][9] = temp[i];
				
				fread (temp,sizeof(double),nbody,dat);//x
				for (i=0;i<nbody;i++) particle[i][1] = Rscale*temp[i];
				
				fread (temp,sizeof(double),nbody,dat);//y
				for (i=0;i<nbody;i++) particle[i][2] = Rscale*temp[i];
				
				fread (temp,sizeof(double),nbody,dat);//z
				for (i=0;i<nbody;i++) particle[i][3] = Rscale*temp[i];
				
				fread (temp,sizeof(double),nbody,dat);//vx
				for (i=0;i<nbody;i++) particle[i][4] = vscale*temp[i];
				
				fread (temp,sizeof(double),nbody,dat);//vy
				for (i=0;i<nbody;i++) particle[i][5] = vscale*temp[i];
				
				fread (temp,sizeof(double),nbody,dat);//vz
				for (i=0;i<nbody;i++) particle[i][6] = vscale*temp[i];
				
				fread (temp,sizeof(double),nbody,dat);//phif		
				for (i=0;i<nbody;i++) particle[i][7] = temp[i];
				
				
				
				
				//CoM correction   //Be careful when SMBH is implemented
				printf("\nApplying centre-of-mass correction.\n");
				double cmr[7];
				for (j=0; j<7; j++) cmr[j] = 0.0;
				
				for (j=0; j<nbody; j++) {
					for (i=1;i<7;i++) 
						cmr[i] += particle[j][0]*particle[j][i]; 
				} 
				
				for (j=0; j<nbody; j++) {
					for (i=1;i<7;i++)
						particle[j][i] -= cmr[i];
				}
				
				
				
				
				//here goes the analysis part
				
				
				//total mass
				double M = 0.0;
				for (j=0;j<nbody;j++) M += particle[j][0];
				if (Mscale == 1.0) printf("\nTotal mass: %g [NBODY]\n",M);
				else printf("\nTotal mass: %g [Msun]\n",M);
				
			      	//triaxiality+Lagrangeradii
				printf("Attention: Determination of Lagrangeradii assumes equal mass particles \n");
				double ha = 0.0, hb = 0.0, hc = 0.0, LAGR = 0;
				int NLAGR;
				NLAGR = 0.01*nbody;
				triax(particle, NLAGR, &ha, &hb, &hc, &LAGR, nbody);
				printf("Main Axis at 1.00 percent Lagrangeradius: %f\t ha = %f\thb = %f\thc = %f\n", LAGR, ha, hb, hc);
				NLAGR = 0.02*nbody;
				triax(particle, NLAGR, &ha, &hb, &hc, &LAGR, nbody);
				printf("Main Axis at 2.00 percent Lagrangeradius: %f\t ha = %f\thb = %f\thc = %f\n", LAGR, ha, hb, hc);
				NLAGR = 0.05*nbody;
				triax(particle, NLAGR, &ha, &hb, &hc, &LAGR, nbody);
				printf("Main Axis at 5.00 percent Lagrangeradius: %f\t ha = %f\thb = %f\thc = %f\n", LAGR, ha, hb, hc);
				NLAGR = 0.1*nbody;
				triax(particle, NLAGR, &ha, &hb, &hc, &LAGR, nbody);
				printf("Main Axis at 10.0 percent Lagrangeradius: %f\t ha = %f\thb = %f\thc = %f\n", LAGR, ha, hb, hc);
				NLAGR = 0.25*nbody;
				triax(particle, NLAGR, &ha, &hb, &hc, &LAGR, nbody);
				printf("Main Axis at 25.0 percent Lagrangeradius: %f\t ha = %f\thb = %f\thc = %f\n", LAGR, ha, hb, hc);
				NLAGR = 0.5*nbody;
				triax(particle, NLAGR, &ha, &hb, &hc, &LAGR, nbody);
				printf("Main Axis at 50.0 percent Lagrangeradius: %f\t ha = %f\thb = %f\thc = %f\n", LAGR, ha, hb, hc);
				NLAGR = 0.6*nbody;
				triax(particle, NLAGR, &ha, &hb, &hc, &LAGR, nbody);
				printf("Main Axis at 60.0 percent Lagrangeradius: %f\t ha = %f\thb = %f\thc = %f\n", LAGR, ha, hb, hc);
				NLAGR = 0.7*nbody;
				triax(particle, NLAGR, &ha, &hb, &hc, &LAGR, nbody);
				printf("Main Axis at 70.0 percent Lagrangeradius: %f\t ha = %f\thb = %f\thc = %f\n", LAGR, ha, hb, hc);
				NLAGR = 0.8*nbody;
				triax(particle, NLAGR, &ha, &hb, &hc, &LAGR, nbody);
				printf("Main Axis at 80.0 percent Lagrangeradius: %f\t ha = %f\thb = %f\thc = %f\n", LAGR, ha, hb, hc);
				NLAGR = 0.9*nbody;
				triax(particle, NLAGR, &ha, &hb, &hc, &LAGR, nbody);
				printf("Main Axis at 90.0 percent Lagrangeradius: %f\t ha = %f\thb = %f\thc = %f\n", LAGR, ha, hb, hc); 
				
				
				//radial profile
				double Rh2D, Rh3D;
				radial_profile(particle, nbody, M, Rscale, &Rh2D, &Rh3D); 
				
				//velocity distribution
 				beta(particle, Rscale, nbody); //NEW
				
				
				lPos = ftell(dat);
				//printf("%li bytes of %li\n",lPos,lSize);
				
			} while (lPos+1<=lSize);
			
			free(temp);
			
			for (j=0;j<NMAX;j++) free (particle[j]);
			free(particle);
			
			
		} else {
			
			printf("File %s is empty!\n",fname);
			printf("##########################################################\n");		
			
		}
		
		fclose(dat);
		
	} else {
		
		printf("\nUsage: \n\n > galanize modelname [Mscale Rscale]\n\n");		
		return 0;
		
	}
	
	t2 = clock();             
	printf("\nElapsed time: %g sec\n",(double)(t2-t1)/CLOCKS_PER_SEC); //write elapsed time to STDOUT
	
    return 0;
}



void radial_profile(double **star, int N, double M, double Rscale, double *Rh2D, double *Rh3D) {
	int i, j;
	*Rh2D = 0.0;
	*Rh3D = 0.0;
	double Mtemp;
	double **rarray;
	rarray = (double **)calloc(N,sizeof(double *));
	for (j=0;j<N;j++){
		rarray[j] = (double *)calloc(2,sizeof(double));
		if (rarray[j] == NULL) {
			printf("\nMemory allocation failed!\n");
			return;
		}
	} 
	double **rarray2D;
	rarray2D = (double **)calloc(N,sizeof(double *));
	for (j=0;j<N;j++){
		rarray2D[j] = (double *)calloc(2,sizeof(double));
		if (rarray2D[j] == NULL) {
			printf("\nMemory allocation failed!\n");
			return;
		}
	} 
	
	for (j=0; j<N; j++) {
		rarray[j][0] = sqrt(star[j][1]*star[j][1]+star[j][2]*star[j][2]+star[j][3]*star[j][3]);
		rarray[j][1] = star[j][0];
		rarray2D[j][0] = sqrt(star[j][1]*star[j][1]+star[j][2]*star[j][2]);//x-y
		rarray2D[j][1] = star[j][0];
	}
	shellsort_reverse(rarray, N, 2);	
	shellsort_reverse(rarray2D, N, 2);	
	
	int noofradialbins = 100;//20
	double rprofile[noofradialbins][7];
	double rmax = 10.0*Rscale; //pc
	double rmin = 0.1*Rscale;
	double stepsize;
	stepsize = (log10(rmax)-log10(rmin))/(noofradialbins-1);
	
	for (j=0;j<noofradialbins;j++) {
		rprofile[j][0] = pow(10.0, log10(rmin) + stepsize*j); //radius
		if (j == 0)
			rprofile[j][1] =  4.0/3.0*PI*pow(rprofile[j][0],3); //volume
		else
			rprofile[j][1] = 4.0/3.0*PI*(pow(rprofile[j][0],3) - pow(rprofile[j-1][0],3)); //volume
		rprofile[j][2] = 0; //number
		rprofile[j][3] = 0; //mass
		if (j == 0)
			rprofile[j][4] =  PI*pow(rprofile[j][0],2); //area
		else
			rprofile[j][4] = PI*(pow(rprofile[j][0],2) - pow(rprofile[j-1][0],2)); //area
		rprofile[j][5] = 0; //number (2D)
		rprofile[j][6] = 0; //mass (2D)
	}
	
	j = 0; i = 0; Mtemp = 0.0;
	while ((j < noofradialbins) && (i < N)) {
		if (rarray[i][0] < rprofile[j][0]) {//3D binning in astrophysical units
			rprofile[j][2] += 1.0;
			rprofile[j][3] += rarray[i][1];
			Mtemp += rarray[i][1];
			if ((Mtemp>=0.5*M) && !(*Rh3D)) *Rh3D = rarray[i][0];
			i++;
		} else {
			j++;
		}
	}
	
	j = 0; i = 0; Mtemp = 0.0;
	while ((j < noofradialbins) && (i < N)) {
		if (rarray2D[i][0] < rprofile[j][0]) {//2D binning in astrophysical units
			rprofile[j][5] += 1.0;
			rprofile[j][6] += rarray2D[i][1];
			Mtemp += rarray2D[i][1];
			if ((Mtemp>=0.5*M) && !(*Rh2D)) *Rh2D = rarray2D[i][0];
			i++;
		} else {
			j++;
		}
	}
	
	//print radial density profile to screen	
	printf("\nRadial density profile:\n\n#  R [pc]   N [1/pc^3]      M [Msun/pc^3]      N [1/pc^2]      M [Msun/pc^2]\n");
	for (j=0;j<noofradialbins;j++) printf("%9.4f  %14.4f  %14.8f  %14.4f  %14.8f\n",rprofile[j][0],rprofile[j][2]/rprofile[j][1],rprofile[j][3]/rprofile[j][1],rprofile[j][5]/rprofile[j][4],rprofile[j][6]/rprofile[j][4]);
	
	for (j=0;j<N;j++) free (rarray[j]);//NEW
	free(rarray);//NEW

	for (j=0;j<N;j++) free (rarray2D[j]);//NEW
	free(rarray2D);//NEW
	
}

void triax(double **particle, int NLAGR, double *ha, double *hb, double *hc, double *LAGRANGE, int nbody) {
	double r,lagr;
	int i;
	int nrot = 0;
	
	
	
	double *dis;
	dis = (double *)calloc(nbody,sizeof(double));
	for (i=0;i<nbody;i++) {
		dis[i] = sqrt(particle[i][1]*particle[i][1]+particle[i][2]*particle[i][2]+particle[i][3]*particle[i][3]);
	}
	
	
	shellsort_reverse_1d(dis, nbody);
	//printf("Lagrangeradii  %8.4f\n",dis[NLAGR]);  
	lagr = dis[NLAGR];
	*LAGRANGE = lagr;
	
	float **am_r;
	am_r = (float **)calloc(4, sizeof(float *));
	for (i=0;i<4;i++) am_r[i] = (float *)calloc(4,sizeof(float));
	
	float **v;
	v = (float **)calloc(4, sizeof(float *));
	for (i=0;i<4;i++) v[i] = (float *)calloc(4,sizeof(float));
	
	float d[4];
	
	am_r[1][1] =  0.0;
	am_r[1][2] =  0.0; 
	am_r[1][3] =  0.0;
	am_r[2][1] =  0.0;
	am_r[2][2] =  0.0;
	am_r[2][3] =  0.0;
	am_r[3][1] =  0.0;
	am_r[3][2] =  0.0;
	am_r[3][3] =  0.0;
	int count = 0;
	
	for (i=0;i<nbody;i++) {
		if (particle[i][0] > 0.0) { 
			r = sqrt(particle[i][1]*particle[i][1]+particle[i][2]*particle[i][2]+particle[i][3]*particle[i][3]);
			if (r < lagr) {
				count++;
				am_r[1][1] += particle[i][2]*particle[i][2]+particle[i][3]*particle[i][3]; 
				am_r[1][2] -= particle[i][1]*particle[i][2];
				am_r[1][3] -= particle[i][1]*particle[i][3];
				am_r[2][1] -= particle[i][2]*particle[i][1];
				am_r[2][2] += particle[i][1]*particle[i][1]+particle[i][3]*particle[i][3];
				am_r[2][3] -= particle[i][2]*particle[i][3];
				am_r[3][1] -= particle[i][3]*particle[i][1];
				am_r[3][2] -= particle[i][3]*particle[i][2];
				am_r[3][3] += particle[i][1]*particle[i][1]+particle[i][2]*particle[i][2];
			}
		} 
	}
	
	jacobi(am_r, 3, d, v, &nrot);
	eigsrt(d, v, 3);                               
	
	//Comparison to homogenous ellipsoid	
	*hc = sqrt(5.0/2.0*(d[2]+d[3]-d[1]))/sqrt(5.0/2.0*(d[1]+d[2]-d[3]));
	*hb = sqrt(5.0/2.0*(d[3]+d[1]-d[2]))/sqrt(5.0/2.0*(d[1]+d[2]-d[3]));
	*ha = sqrt(5.0/2.0*(d[1]+d[2]-d[3]))/sqrt(5.0/2.0*(d[1]+d[2]-d[3]));    
	
	
	for (i=0;i<4;i++) free(am_r[i]);
	free(am_r);
	
	for (i=0;i<4;i++) free(v[i]);
	free(v);
	
	free(dis); //NEW
	
}

void beta(double **particle, double scalesize, int nbody) {
	int i,j,NR;
	int l,k,var; //for Bootstrapping
	double er1,er2; //for Bootstrapping
	int JMAX = 40;
	double theta,phi,sr2,sp2,st2;
	double *dis,beta1 = 0.0 ,beta2 = 0.0;
	
	printf("\n");
	printf("Anisotropy parameter beta Method1: 1-(st^2+sp^2)/(2*sr^2) Method2: 1-(st^2)/(sr^2)  \n");
	printf("Attention: Assumes nonrotating model, Error evaluated with Bootstrapping method \n");
      	printf("\n");
	dis = (double *)calloc(nbody,sizeof(double));
	for (i=0;i<nbody;i++) {
		dis[i] = sqrt(particle[i][1]*particle[i][1]+particle[i][2]*particle[i][2]+particle[i][3]*particle[i][3]);
	}
	
	
	for (j=0;j<JMAX;j++) { //Main Loop    
		sr2 = 0.0;
		sp2 = 0.0;
		st2 = 0.0;
		NR  = 0;
		for (i=0;i<nbody;i++) { //Loop over all particles
			if ((dis[i] < pow(10.0,-3.0+0.1*(j+1))*scalesize) && (dis[i] > pow(10.0,-3.0+0.1*j)*scalesize)) {
			    theta = acos(particle[i][3]/dis[i]);
                       	    phi = atan2(particle[i][2],particle[i][1]);
			    NR  += 1;
			    sr2 += pow(particle[i][4]*sin(theta)*cos(phi) + particle[i][5]*sin(theta)*sin(phi) + particle[i][6]*cos(theta),2.0);        
		            sp2 += pow(particle[i][4]*cos(theta)*cos(phi) + particle[i][5]*cos(theta)*sin(phi) - particle[i][6]*sin(theta),2.0);
			    st2 += pow(particle[i][5]*cos(phi) - particle[i][4]*sin(phi),2.0);
			}		
		}
		beta1 = 1.0 - (st2+sp2)/(2.0*sr2); //Method1
		beta2 = 1.0 - (st2)/(sr2); //Method1
		
		//*****Bootstrapping******
		var=0;
		double **aa;
		aa = (double **)calloc(NR, sizeof(double *));
		for (k=0;k<NR;k++) aa[k] = (double *)calloc(6,sizeof(double));
		
		for (i=0;i<nbody;i++){ 	         
			if ((dis[i] < pow(10.0,-3.0+0.1*(j+1))*scalesize) && (dis[i] > pow(10.0,-3.0+0.1*j)*scalesize)) {
				for (l=0;l<6;l++) aa[var][l] = particle[i][l+1];
				var += 1;
			}
		}
		bootstrap(aa, NR, &er1, &er2);
		for (k=0;k<NR;k++) free(aa[k]);
		free(aa);
		//*********END*******/
		printf("r: %11.5f\t  beta[Method1]: %11.5f\t pm %11.5f\t beta[Method2]: %11.5f\t pm %11.5f\t total N in shell %i\n ",
		  pow(10.0,-3.0+0.1*(j+1))*scalesize, beta1, er1, beta2, er2, NR);
		//printf("r: %11.5f\t sr**2 %11.5f\t st**2 %11.5f\n",pow(10.0,-3.0+0.1*(j+1))*scalesize,sr2/NR,(st2+sp2)/NR);
 

	}	
	
	free(dis);//NEW
	
}

void bootstrap(double **aa, int n, double *er1, double *er2) {
    int i,k,l,m,a;
    time_t t;
    double xxx,sr2,sp2,st2;
    double theta,phi,*beta1,*beta2,xdummy; 
    // time(&t);
    srand48((unsigned int)1); //NEW             /* Zufallsgenerator initialisieren */
    int LMAX = 100;
    beta1 = (double *)calloc(LMAX,sizeof(double));
    beta2 = (double *)calloc(LMAX,sizeof(double));
	
	
    for (m=0;m<LMAX;m++){
		double **arr;
		arr = (double **)calloc(n, sizeof(double *));
		for (k=0;k<n;k++) arr[k] = (double *)calloc(6,sizeof(double));
		
		for (i=0;i<n;i++){ 
			do{//NEW
				xxx = (1.*n*drand48());//NEW
				a   = (int)xxx;
			} while (a>=n || a<0);//NEW
			for (l=0;l<6;l++) arr[i][l] = aa[a][l];
		} 
		sr2 = 0.0;
		sp2 = 0.0;
		st2 = 0.0;
		for (i=0;i<n;i++){  //Berechnung von beta 
       	                theta = acos(arr[i][2]/sqrt(arr[i][0]*arr[i][0]+arr[i][1]*arr[i][1]+arr[i][2]*arr[i][2]));
                        phi = atan2(arr[i][1],arr[i][0]);
			sr2 += pow(arr[i][3]*sin(theta)*cos(phi) + arr[i][4]*sin(theta)*sin(phi) + arr[i][5]*cos(theta),2.0);            
			sp2 += pow(arr[i][3]*cos(theta)*cos(phi) + arr[i][4]*cos(theta)*sin(phi) - arr[i][5]*sin(theta),2.0);
			st2 += pow(arr[i][4]*cos(phi) - arr[i][3]*sin(phi),2.0); 
		} 
      		beta1[m] = 1.0 - (st2+sp2)/(2.0*sr2); //Method1
		beta2[m] = 1.0 - (st2)/(1.0*sr2); //Method2
      
		for (k=0;k<n;k++) free(arr[k]);
		free(arr);
    }
    
	
	
    // 2nd Part: Calculate 1 sigma error from Distribution 
	
	double mu1,mu2,sigma1,sigma2;
	
	mu1 = 0.0;
	sigma1 = 0.0;
	mu2 = 0.0;
	sigma2 = 0.0;
	
	for (m=0;m<LMAX;m++) {
		mu1 += beta1[m]/(1.0*LMAX); //NEW
		mu2 += beta2[m]/(1.0*LMAX); //NEW
	}
	for (m=0;m<LMAX;m++){
		sigma1 += (beta1[m]-mu1)*(beta1[m]-mu1);  
		sigma2 += (beta2[m]-mu2)*(beta2[m]-mu2); 
	}
	*er1 = sqrt(sigma1/(LMAX-1.0));//NEW
	*er2 = sqrt(sigma2/(LMAX-1.0));//NEW
	
	free(beta1); //NEW
	free(beta2); //NEW
	
}



///****************************************************************************************
/* ---------- sorting functions ------------ */
void shellsort(double **array, int N, int k) {//largest up
	int i,j,n,o;
	N = N-1;
	double swap[k];
	//guess distance n
	for (n = 1; n <= N/9; n = 3*n+1);
	for (; n > 0; n /= 3) {
		for (i = n; i <= N; i++) {
			for (o=0; o<k; o++) swap[o] = array[i][o];
			for (j = i; ((j >= n) && (array[j-n][0] < swap[0])); j -= n) {
				for (o=0; o<k; o++) array[j][o] = array[j-n][o];
			}
			for (o=0; o<k; o++) array[j][o] = swap[o];
		}
	}
}

void shellsort_reverse(double **array, int N, int k) {//smallest up
	int i,j,o,n;
	N = N-1;
	double swap[k];
	//guess distance n
	for (n = 1; n <= N/9; n = 3*n+1);
	for (; n > 0; n /= 3) {
		for (i = n; i <= N; i++) {
			for (o=0; o<k; o++) swap[o] = array[i][o];
			for (j = i; ((j >= n) && (array[j-n][0] > swap[0])); j -= n) {
				for (o=0; o<k; o++) array[j][o] = array[j-n][o];
			}
			for (o=0; o<k; o++) array[j][o] = swap[o];
		}
	}
}

void shellsort_1d(double *array, int N) {//largest up
	int i,j,n;
	N = N-1;
	double swap;
	//guess distance n
	for (n = 1; n <= N/9; n = 3*n+1);
	for (; n > 0; n /= 3) {
		for (i = n; i <= N; i++) {
			swap = array[i];
			for (j = i; ((j >= n) && (array[j-n] < swap)); j -= n) {
				array[j] = array[j-n];
			}
			array[j] = swap;
		}
	}
}

void shellsort_reverse_1d(double *array, int N) {//smallest up
	int i,j,n;
	N = N-1;
	double swap;
	//guess distance n
	for (n = 1; n <= N/9; n = 3*n+1);
	for (; n > 0; n /= 3) {
		for (i = n; i <= N; i++) {
			swap = array[i];
			for (j = i; ((j >= n) && (array[j-n] > swap)); j -= n) {
				array[j] = array[j-n];
			}
			array[j] = swap;
		}
	}
}


/*----------- Numerical recipes ------------- */


void jacobi(float **am, int n, float d[], float **v, int *nrot) {
	
	int j,iq,ip,i;
	float tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
	
	b=vector(1,n);
	z=vector(1,n);
	
	for (ip=1;ip<=n;ip++) {
		for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	
	for (ip=1;ip<=n;ip++) {
		b[ip]=d[ip]=am[ip][ip];
		z[ip]=0.0;
	}
	
    *nrot=0;
	
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++)
				sm += fabs(am[ip][iq]);
		}
		if (sm == 0.0) {
			free_vector(z,1,n);
			free_vector(b,1,n);
			return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++) {
				g=100.0*fabs(am[ip][iq]);
				if (i > 4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip]) && (float)(fabs(d[iq])+g) == (float)fabs(d[iq])) am[ip][iq]=0.0;
				else if (fabs(am[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((float)(fabs(h)+g) == (float)fabs(h)) t=(am[ip][iq])/h;
					else {
						theta=0.5*h/(am[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*am[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					am[ip][iq]=0.0;
					
					for (j=1;j<=ip-1;j++) {
						ROTATE(am,j,ip,j,iq);
					}
					
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(am,ip,j,j,iq);
					}
					
					for (j=iq+1;j<=n;j++) {
						ROTATE(am,ip,j,iq,j);
					}
					
					for (j=1;j<=n;j++) {
						ROTATE(v,j,ip,j,iq);
					}
					
					++(*nrot);
					
				}
			}
		}
		
		for (ip=1;ip<=n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	
	nrerror("Too many iterations in routine jacobi");
	
}

void eigsrt(float d[], float **v, int n) {
	
	int k,j,i;
	float p;
	for (i=1;i<n;i++) {
		p=d[k=i];
		for (j=i+1;j<=n;j++) if (d[j] >= p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			for (j=1;j<=n;j++) {
				p=v[j][i];
				v[j][i]=v[j][k];
				v[j][k]=p;
			}
		}
	}
	
}

float *vector(long nl, long nh) {
	/* allocate a float vector with subscript range v[nl..nh] */
	
	float *v;
	
	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
	
	return v-nl+NR_END;
	
}

void free_vector(float *v, long nl, long nh) {
	/* free a float vector allocated with vector() */
	
	free((FREE_ARG) (v+nl-NR_END));
	
}

void nrerror(char error_text[]) {
	/* Numerical Recipes standard error handler */
	
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	
	exit(1);
	
}

