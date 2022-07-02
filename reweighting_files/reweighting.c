#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define n_frame 60000
#define n_system 16
#define n_res 1
#define n_space 24 
#define n_gap 15

void CMAPSPL (int dx, double *y, int n, double output[]){
	int i;
	double y2[2 * n_space];
	double u[2 * n_space];
	double pinv;
	double dxinv = 1.0 / dx;
	memset(y2, 0, sizeof(y2));
	memset(u, 0, sizeof(u));
	y2[0] = 0.0;
	u[0] = 0.0;

	for(i = 1;i < n - 1;i++){
		pinv = 1.0 / (y2[i - 1] + 4.0);
		y2[i] = -1.0 * pinv;
		u[i] = ((6 * y[i + 1] - 12 * y[i] + 6 * y[i - 1]) * dxinv * dxinv - u[i - 1]) * pinv;
	}

	y2[n - 1] = 0.0;
	u[n - 1] = 0.0;

	for(i = n - 2;i >= 0;i--){
		y2[i] = y2[i] * y2[i + 1] + u[i];
	}
	//print(y2);	
	for(i = 0;i < n;i++){
		//printf("%.4f\n",y2[i]);
		output[i] = y2[i];
	}
}

void CMAPSPI(int xmin, int dx, double *ya, double *y2a, int x, double *y, double *y1){
	int inx = floor((x - xmin) / dx);
	double a = (xmin + inx * dx + dx - x) / dx;
	double b = (x - xmin - inx * dx) / dx;
	*y = a * ya[inx] + b * ya[inx + 1] + ((a * a * a - a) * y2a[inx] + (b * b * b - b) * y2a[inx 
+ 1]) * (dx * dx) / 6.0;
	*y1 = (ya[inx + 1] - ya[inx]) / dx - (3 * a * a - 1) / 6 * dx * y2a[inx] + (3 * b * b - 1) / 
6 * dx * y2a[inx + 1];
}

void setcmap2(double cmap[n_space][n_space], double ggrd0[n_space][n_space], double ggrd1[n_space][n_space], double ggrd2[n_space][n_space], double ggrd3[n_space][n_space]){
	int i;
	int j;
	int ii = 0;
	int jj = 0;
	int num = n_space;
	int xm = num / 2;
	int dnum = 2 * num;
	int dx = 360 / num;
	int xmin = -360;
	double tgmap[2 * n_space][2 * n_space];

	for(i = 0;i < num;i++){
		for(j = 0;j < num;j++){
			ggrd0[i][j] = cmap[i][j];
		}
	}

	for(i = 0;i < dnum;i++){
		ii = (i + xm) % num;
		for(j = 0;j < dnum;j++){
			jj = (j + xm) % num;
			tgmap[i][j] = cmap[jj][ii];
		}
	}

	/*
	// test hliu
	printf("tgmap:\n");
	for(i = 0; i < dnum; i++){
		for(j = 0; j < dnum; j++){
			printf("%d\t%d\t%.4f\n", i, j, tgmap[i][j]);
		}
	}
	*/
	
	double y2temp[2 * n_space];
	double y2a[2 * n_space][2 * n_space];
	for(j = 0;j < dnum;j++){
		memset(y2temp, 0, sizeof(y2temp));
		CMAPSPL(dx, tgmap[j], dnum, y2temp);
		for (i = 0;i < dnum;i++){
			y2a[j][i] = y2temp[i];
		}
	}

	/*
	// test hliu
	printf("y2a:\n");
	for(i = 0; i < dnum; i++){
		for(j = 0; j < dnum; j++){
			printf("%d\t%d\t%.4f\n", i, j, y2a[i][j]);
		}
	}
	*/

	int phi;
	int psi;
	double a = 0.0;
	double b = 0.0;
	double v = 0.0;
	double v1 = 0.0;
	double v2 = 0.0;
	double v12 = 0.0;
	double yytmp[2 * n_space];
	double y1tmp[2 * n_space];
	double u2[2 * n_space];
	for(ii = xm;ii < xm + num;ii++){
		phi = (ii - xm) * dx - 180;
		for(jj = xm;jj < xm + num;jj++){
			psi = (jj - xm) * dx - 180;

			for(j = 0;j < dnum;j++){
				CMAPSPI(-360, dx, tgmap[j], y2a[j], psi, &a, &b);
				yytmp[j] = a;
				y1tmp[j] = b;
			}
				
			CMAPSPL(dx ,yytmp, dnum, u2);
			CMAPSPI(-360, dx, yytmp, u2, phi, &v, &v1);
			CMAPSPL(dx ,y1tmp, dnum, u2);
			CMAPSPI(-360, dx, y1tmp, u2, phi, &v2, &v12);
				
			ggrd1[jj - xm][ii - xm] = v1 * dx;
			ggrd2[jj - xm][ii - xm] = v2 * dx;
			ggrd3[jj - xm][ii - xm] = v12 * dx * dx;
		}
	}
	
	/*
	// test hliu
	printf("ggrd:\n");
	for(i = 0;i<num;i++){
		for(j=0;j<num;j++){
			printf("%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\n",i,j,ggrd0[i][j],ggrd1[i][j],ggrd2[i][j],ggrd3[i][j]);
		}
	}
	*/
}

double ecmapfast(double ggrd0[n_space][n_space], double ggrd1[n_space][n_space], double ggrd2[n_space][n_space], double ggrd3[n_space][n_space], int IP1, int IP2, int IP1P1, int IP2P1, double TT, double TU, int wt[16][16]){

	int inn = 0;
	int i;
	int j;
	int k;
	double xx;
	double TC[4][4];

	for(i = 0;i < 4;i++){
		for (j = 0;j < 4;j++){

			xx = 0.0;

			xx += wt[inn][0] * ggrd0[IP2][IP1];
			xx += wt[inn][1] * ggrd0[IP2][IP1P1];
			xx += wt[inn][2] * ggrd0[IP2P1][IP1P1];
			xx += wt[inn][3] * ggrd0[IP2P1][IP1];
			xx += wt[inn][4] * ggrd1[IP2][IP1];
			xx += wt[inn][5] * ggrd1[IP2][IP1P1];
			xx += wt[inn][6] * ggrd1[IP2P1][IP1P1];
			xx += wt[inn][7] * ggrd1[IP2P1][IP1];
			xx += wt[inn][8] * ggrd2[IP2][IP1];
			xx += wt[inn][9] * ggrd2[IP2][IP1P1];
			xx += wt[inn][10] * ggrd2[IP2P1][IP1P1];
			xx += wt[inn][11] * ggrd2[IP2P1][IP1];
			xx += wt[inn][12] * ggrd3[IP2][IP1];
			xx += wt[inn][13] * ggrd3[IP2][IP1P1];
			xx += wt[inn][14] * ggrd3[IP2P1][IP1P1];
			xx += wt[inn][15] * ggrd3[IP2P1][IP1];

			inn++;
			TC[i][j] = xx;
		}
	}
	double E = 0.0;
	int II;

	for(II = 3;II > -1;II--){
		E = TT * E + ((TC[II][3] * TU + TC[II][2]) * TU + TC[II][1]) * TU + TC[II][0];
	}
	

	return E;
}

double RmsCmap(double cmap1[n_space][n_space], double cmap2[n_space][n_space]){
	static double rms = 0.0;
	int i;
	int j;
	for(i=0; i<n_space; i++){
		for(j=0; j<n_space; j++){
			rms += pow((cmap1[i][j] - cmap2[i][j]), 2);
		}
	}
	rms = pow(rms/n_space/n_space, 0.5);
	return rms;
}

double SumDouble(double *array, int len){
	double sum = 0.0;
	int i;
	for(i=0; i<len; i++){
		sum += array[i];
	}
	return sum;
}

double SumDoubleWeigted(double *array, double *weight, int len){
	double sum = 0.0;
	int i;
	for(i=0; i<len; i++){
		sum += array[i] * weight[i];
	}
	return sum;
}

// for MergeSort
void merge(double *array, int begin, int mid, int end){
	int i, j, k;
	int n1 = mid - begin + 1;
	int n2 = end - mid;

	double left[n1], right[n2];

	for(i = 0; i < n1; i++){
		left[i] = array[begin+i];
	}
	for(j = 0; j < n2; j++){
		right[j] = array[mid+1+j];
	}

	i = 0;
	j = 0;
	k = begin;
	while (i < n1 && j < n2){
		if(left[i] <= right[j]){
			array[k] = left[i];
			i++;
		}
		else{
			array[k] = right[j];
			j++;
		}
		k++;
	}

	while(i < n1){
		array[k] = left[i];
		i++;
		k++;
	}

	while(j < n2){
	 array[k] = right[j];
		j++;
		k++;
	}
}

void MergeSort(double *array, int begin, int end){
	if(begin < end){
		int mid = begin + (end - begin) / 2;
		//printf("%d\n", mid);

		MergeSort(array, begin, mid);
		MergeSort(array, mid+1, end);

		merge(array, begin, mid, end);
	}
}
// for MergeSort

double *Weighting(double *e_init, double *e_new){
	double kJ2kcal = 4.1858518;
	double invkT = 1.0/(0.0019872041 * 300);
	static double w[n_system*n_frame];
	int i;
	for(i = 0; i<n_system*n_frame; i++){
		w[i] = exp((e_init[i] - e_new[i]) * invkT / kJ2kcal);
	}
	return w;
}

int WeightingCheck(double *w, int begin, int end){
// top x% should have less than y% weight

	int len = end - begin;
	double w_new[len];
	double top_percent = 0.05;
	double cutoff = 0.3;
	int top = ceil(top_percent * len);
	int i;
	
	for(i = 0; i<len; i++){
		w_new[i] = w[begin+i];
	}

	MergeSort(w_new, 0, len-1);
	
	double w_part[top];
	for(i=0; i<top; i++){
		w_part[i] = w_new[len-1-i];
	}
	double c = SumDouble(w_part, top) / SumDouble(w_new, len);
	
	if(c < cutoff){
		return 1;
	}
	else{
		return 0;
	}
}

double *JcoupPred(double *jcoup, double *weight){
	static double jcoup_new[n_system];
	double jcoup_sum = 0.0;
	double w_sum = 0.0;
	int i;
	for(i=0; i<n_system*n_frame; i++){
		jcoup_sum += jcoup[i] * weight[i];
		w_sum += weight[i];
		if(((i+1) % n_frame) == 0){
			jcoup_new[i/n_frame] = jcoup_sum / w_sum;
			jcoup_sum = 0.0;
			w_sum = 0.0;
		}
	}
	return jcoup_new;
}

int main(int argc, char **argv){
	srand((unsigned)time(NULL));
	
	if(argc != 7){
		printf("Usage: %s file_CMAP file_dih file_jcoup file_out n_step w_rmsd\n", argv[0]);
		exit(0);
	}

	char *file_cmap = argv[1];
	char *file_dih = argv[2];
	char *file_jcoup = argv[3];
	char *file_out = argv[4];
	int n_step = atoi(argv[5]);
	double w_rmsd = atof(argv[6]);
	int wt[16][16] = {
		{   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  },
		{   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0  },
		{  -3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0  },
		{   2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0  },
		{   0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  },
		{   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0  },
		{   0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1  },
		{   0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1  },
		{  -3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  },
		{   0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0  },
		{   9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2  },
		{  -6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2  },
		{   2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  },
		{   0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0  },
		{  -6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1  },
		{   4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1  }
	};

	// Step 1. Read files
	// read cmap file	
	FILE *f_cmap = fopen(file_cmap, "r");
	double cmap_trans[n_space][n_space];
	int i = 0;
	int j = 0;
	int k = 0;
	int n1 = -1;
	int n2 = 0;
	int ch;
	char value[20];
	char line[100];

	while((ch = fgetc(f_cmap)) != EOF){
		line[i] = ch;
		i++;
		if(ch == '\n'){
			line[i] = '\0';
			if(line[0] == '!'){
				n1++;
				n2 = 0;
			}
			else if(strcmp(line, "") == 0){
				continue;
			}
			else{
				for(j=0; j<i-1; j++){
					value[k] = line[j];
					k++;
					if(k == 8){
						value[k] = '\0';
						//printf("%d\t%d\t%.4f\n", n1, n2, atof(value));
						k = 0;
						cmap_trans[n1][n2] = atof(value);
						j++;
						n2++;
					}
				}
			}
			i = 0;
		}
	}
	//printf("%.4f\t%.4f\n", cmap[0][0], cmap[23][23]);
	double cmap[n_space][n_space];
	for(i = 0; i < n_space; i++){
		for(j = 0; j < n_space; j++){
			cmap[j][i] = cmap_trans[i][j];
		}
	}
	
	// test hliu
	//printf("cmap:\n");
	//for(i = 0; i < n_space; i++){
	//	for(j = 0; j < n_space; j++){
	//		printf("%d\t%d\t%.4f\n", i, j, cmap[i][j]);
	//	}
	//}
	

	// read dih file
	i = 0;
	j = 0;
	FILE *f_dih = fopen(file_dih, "r");
	double phi[n_system * n_frame];
	double psi[n_system * n_frame];
	while((ch = fgetc(f_dih)) != EOF){
		value[i] = ch;
		i++;
		if(ch == '\t'){
			value[i] = '\0';
			i = 0;
			phi[j] = atof(value);
		}
		else if(ch == '\n'){
			value[i] = '\0';
			i = 0;
			psi[j] = atof(value);
			j++;
		}
	}
	//printf("phi%d\t%d\t%d\t%d\n",&phi[639999], &psi[639999], &phi[-1], &psi[-1]);
	
	

	// test hliu
	//printf("dihedrals:\n");
	//for(i = 0; i < n_system*n_frame; i++){
	//	printf("%d\t%.4f\t%.4f\n", i, phi[i], psi[i]);
	//}

	// read jcoup file
	i = 0;
	j = 0;
	double jcoup[n_system * n_frame];
	FILE *f_jcoup = fopen(file_jcoup, "r");
	while((ch = fgetc(f_jcoup)) != EOF){
		value[i] = ch;
		i++;
		if(ch == '\n'){
			value[i] = '\0';
			i = 0;
			jcoup[j] = atof(value);
			j++;
		}
	}
	
	//printf("%.4f\t%.4f\n", jcoup[0], jcoup[n_system*n_frame-1]);

	// Step2. Initialize
	// initial target
	// target = t1 + w_rmsd * t2
	// t1: rms(jcoup)
	// t2: rms(cmap)
	double jcoup_exp[n_system] = {6.075, 6.495, 7.318, 7.177, 6.931, 7.547, 6.615, 6.991,
						 	   7.495, 7.351, 7.007, 6.974, 7.642, 7.294, 7.050, 7.103};
	double jcoup_avg[n_system];
	double jcoup_best[n_system];
	double *jcoup_pred;
	
	double tmp = 0.0;
	j = 0;
	for(i=0; i<n_system*n_frame; i++){
		tmp += jcoup[i];
		if((i+1) % n_frame == 0){
			tmp /= n_frame;
			jcoup_avg[j] = tmp;
			j++;
			tmp = 0.0;
		}
	}
	double t1 = 0.0;
	double t2 = 0.0;
	double t1_best = 0.0;
	double t2_best = 0.0;
	double target_init;
	double target_old;
	double target_new;
	double target_diff;
	double target_best;
	tmp = 0.0;
	for(i=0; i<n_system; i++){
		tmp += pow((jcoup_avg[i] - jcoup_exp[i]), 2.0);
		//printf("%.4f\t%.4f\t%.4f\n", jcoup_avg[i], jcoup_exp[i], tmp);
	}
	t1 = pow((tmp/n_system), 0.5);

	double cmap_old[n_space][n_space];
	double cmap_new[n_space][n_space];
	double cmap_best[n_space][n_space];
	for(i=0; i<n_space; i++){
		for(j=0; j<n_space; j++){
			cmap_old[i][j] = cmap[i][j];
			cmap_new[i][j] = cmap[i][j];
		}
	}
	t2 = RmsCmap(cmap, cmap_new);
	target_init = t1 + w_rmsd * t2;
	t1_best = t1;
	t2_best = t2;
	target_old = target_init;
	target_new = target_init;
	target_best = target_init;
	//printf("Check target: %.4f %.4f\n", target_init, target_new);
	
	// initial cmap energy
	int ii = 0;
	int ij = 0;
	int iphi1;
	int ipsi1;
	int iframe = 0;
	int pai = 180;
	int IP1[n_system*n_frame][n_res], IP2[n_system*n_frame][n_res], IP1P1[n_system*n_frame][n_res], IP2P1[n_system*n_frame][n_res];
	double TT[n_system*n_frame][n_res], TU[n_system*n_frame][n_res];
	for(i=0; i<n_system*n_frame; i++){
		ii = i % n_res;
		iframe = i / n_res ;
		iphi1 = floor((phi[i] + 180.0) / n_gap);
		//printf("%i\n",iphi1);
		ipsi1 = floor((psi[i] + 180) / n_gap);
		TT[iframe][ij] = (phi[i] - iphi1 * n_gap + 180.0) / n_gap * 1.0;
		TU[iframe][ij] = (psi[i] - ipsi1 * n_gap + 180.0) / n_gap * 1.0;
		
		IP1[iframe][ij] = iphi1 % n_space;
		IP2[iframe][ij] = ipsi1 % n_space;
		IP1P1[iframe][ij] = (IP1[iframe][ij] + 1) % n_space;
		IP2P1[iframe][ij] = (IP2[iframe][ij] + 1) % n_space;
		
		ij = (ij + 1) % n_res;
	}
	
	double g10[n_space][n_space];
	double g11[n_space][n_space];
	double g12[n_space][n_space];
	double g13[n_space][n_space];
	setcmap2(cmap, g10, g11, g12, g13);
	
	double energy_init[n_system*n_frame];
	double energy_new[n_system*n_frame];
	for(i=0; i<n_system*n_frame; i++){
		tmp = 0.0;
		for(j=0; j<n_res; j++){
			tmp += ecmapfast(g10, g11, g12, g13, IP1[i][j], IP2[i][j], IP1P1[i][j], IP2P1[i][j], TT[i][j], TU[i][j], wt);
			//printf("%d\t%d\t%d\t%d\t%.4f\t%.4f\n",IP1[i][j], IP2[i][j], IP1P1[i][j], IP2P1[i][j], TT[i][j], TU[i][j]);
		}
		//printf("%.4f\n", tmp);
		energy_init[i] = tmp;
	}

	// test hliu
	//printf("energy_init:\n");
	//for(i = 0; i < n_system*n_frame; i++){
	//	printf("%d\t%.4f\n", i, energy_init[i]);
	//}

	//exit(0);
	
	// Step 3. MCSA(Mote Carlo Simulated Annealing) calculations
	// condition of MCSA
	double temp0 = 10.0;
	double temp;
	int step = 0;
	double prob = 0.0;
	double prob0;
	double boltz;
	int accepted = 1;

	int qweight;
	int begin;
	int end;
	double *weight_new;
	
	printf("n_frame: %d, n_step: %d, w_rmsd: %.4f, target_init: %.4f\n", n_frame, n_step, w_rmsd, target_init);
	printf("%11s%11s%11s%11s%11s%11s%11s%11s\n",
		   "MC step", "Temp", "P",  "Accepted?", "target", "t_best", "t1", "t2");
	while(step < n_step){
		temp = temp0 * exp(-1.0*(step/(n_step/4.0)));
		// modify cmap parameters
		for(i=0; i<n_space; i++){
			for(j=0; j<n_space; j++){
				cmap_new[i][j] = cmap_old[i][j] + rand() * 0.01 / RAND_MAX - 0.005;
			}
		}
		
		// update cmap energy
		setcmap2(cmap_new, g10, g11, g12, g13);
		for(i=0; i<n_system*n_frame; i++){
			tmp = 0.0;
			for(j=0; j<n_res; j++){
				tmp += ecmapfast(g10, g11, g12, g13, IP1[i][j], IP2[i][j], IP1P1[i][j], IP2P1[i][j], TT[i][j], TU[i][j], wt);
			}
			energy_new[i] = tmp;
		}
		weight_new = Weighting(energy_init, energy_new);
		jcoup_pred = JcoupPred(jcoup, weight_new);		

		qweight = 0;
		for(i=0; i<n_system; i++){
			begin = i * n_frame;
			end = (i+1) * n_frame;
			qweight += WeightingCheck(weight_new, begin, end);
			jcoup_best[i] = jcoup_pred[i];
			//printf("Check %d %.4f %.4f\n", i, jcoup_pred[i], jcoup_exp[i]);
		}

		if(qweight == n_system){
			tmp = 0.0;
			for(i=0; i<n_system; i++){
				tmp += pow((jcoup_pred[i] - jcoup_exp[i]), 2.0);
			}
			
			t1 = pow(tmp/n_system, 0.5);
			t2 = RmsCmap(cmap_new, cmap_old);
			target_new = t1 + w_rmsd * t2;
			target_diff = target_new - target_old;
			
			//printf("Check target_new target_diff: %.4f %.4f\n", target_new, target_diff);
			// Metropolis criterion
			if(step == 0){
				prob = 0.0; // required to prevent overflow
			}
			else{
				boltz = -1.0 * target_diff / (0.0019872041 * temp);
				prob = exp(boltz);
			}
			prob0 = rand() * 1.0 / RAND_MAX;
			
			if(target_diff < 0.0){
				accepted = 1;
				prob = 1.0;
			}
			else if(prob0 < prob){
				accepted = 1;
			}
			else{
				accepted = 0;
			}
		}
		//printf("Check accepted: %d\n", accepted);
		// update
		if(accepted){
			target_old = target_new;
			for(i=0; i<n_space; i++){
				for(j=0; j<n_space; j++){
					cmap_old[i][j] = cmap_new[i][j];
				}
			}
			if(target_new < target_best){
				target_best = target_new;
				t1_best = t1;
				t2_best = t2;
				for(i=0; i<n_system; i++){
					jcoup_best[i] = jcoup_pred[i];
				}
				for(i=0; i<n_space; i++){
					for(j=0; j<n_space; j++){
						cmap_best[i][j] = cmap_new[i][j];
					}
				}
			}
		}
		
		if((step % 1) == 0){
			printf("%11d%11.2f%11.4f%11d%11.4f%11.4f%11.4f%11.4f\n", 
				step, temp, prob, accepted, target_new, target_best, t1, t2);
		}
		step++;
	}
	printf("%11d%11.2f%11.4f%11d%11.4f%11.4f%11.4f%11.4f\n",
		step, temp, prob, accepted, target_new, target_best, t1, t2);

	// output best parameters and predicted jcoupling
	printf("Best predicted jcoupling values are:");
	for(i = 0; i<n_system; i++){
		printf(" %.4f", jcoup_best[i]);
	}
	printf("\n");
	printf("Best t1 and t2 are: %.4f, %.4f\n", t1_best, t2_best);
	int aa;
	FILE *f_out = fopen(file_out, "w");
	for(j=0;j<n_space;j++){
		aa=-180+360*j/n_space;
		fprintf(f_out,"!%d\n",aa);
		for(i=0;i<4;i++){
			fprintf(f_out,"%8.2f%8.2f%8.2f%8.2f%8.2f\n",cmap_best[5*i][j],cmap_best[5*i+1][j],cmap_best[5*i+2][j],cmap_best[5*i+3][j],cmap_best[5*i+4][j]);
		}
		fprintf(f_out,"%8.2f%8.2f%8.2f%8.2f\n\n",cmap_best[20][j],cmap_best[21][j],cmap_best[22][j],cmap_best[23][j]);
	}
	fclose(f_out);
	return 0;
	
}
