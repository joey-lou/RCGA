/* Sample RCGA code to find minimum of simple spherical function 
Jiujiu Lou, 06/28/2018

Objective function:
min f(x1,x2,x3) = x1^2+x2^2+(x3-2.5)^2
x1 in [-1,1]
x2 in [-1,1]
x3 in [-1,1]

Methods proposed: 
tournament selection method 
laplace crossover/ BLX-alpha crossover
power mutation
penalty for out of constraint variables 

Possible modification:
instead of having to store two generations, dissecting the current
generation into two parts and keep portion of the current generation
while updating the rest. This will save computation of copying generations
at every iteration. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
/* constants */
const int PPL_SIZE = 30;						// population size (must be even)
const int VAR_SIZE = 3;							// number of variables
const float CO_PROB = 0.8;						// crossover probability
const float MU_PROB = 0.1;						// mutation probability
const int ITER = 100;							// iteration number
const float TOUR_PERCENT = 0.17;					// percentage of population enters tournaments
const int TOUR_SIZE = TOUR_PERCENT*PPL_SIZE;				// number of individuals enter tournaments
const float P = 10; 							// power of mutation
const float uplimit = 1.0;
const float lowlimit = -1.0;
#define PI 3.14159265358979

/* global variables */
float n_gen[PPL_SIZE][VAR_SIZE];					// new generation 2D array
float o_gen[PPL_SIZE][VAR_SIZE];					// old generation 2D array
float ngen_value[PPL_SIZE];						// function values of new generation
float ogen_value[PPL_SIZE];						// function values of old generation
int parent1, parent2;							// parents selected from old generation
int tour_idx[TOUR_SIZE];						// indices of those enter tournaments
float fworst = INFINITY;						// the worst feval result from previous generation
int total_idx[PPL_SIZE];


/* helper function list */
float sum(float *arr, int size);				// returns sum of given array
int max(float *arr, int size);					// returns index of max in array
int min(float *arr, int size, int *idx_arr);			// return minimum from indexed individuals in new gen
int fit(int idx, int entry);					// check if indexed individual fits constraints
float feval(int idx, int fchoice);				// evaluate function value for indexed individual
float penval(int idx);						// evaluate penalty added to indexed individual
void copygen(void);						// copy new genration to old generation
void updatetour(void);						// update indices in tournament
int roulette(void);
float homogeneitycal(void);
float randomgen(float minimum, float maximum);
/* crossover methods */
void laplacex(float a, float b, int idx);			// laplace crossover
void BLX_alpha(float alpha, int idx);
void arithmcross(int idx, float alpha);

/* mutation methods */
void uniformutate(int idx, float delta);
void powermutate(int idx);					// mutate indexed individual in new gen


/* visualization tools */
void printarrayf(float *arr,int size);
void printarrayi(int *arr,int size);
void printgen(float arr[PPL_SIZE][VAR_SIZE],int idx);
void display(int iter, float homogeneity);

int main(int argc, char* argv[]){
	// important to start with all values within constraint
	if (argc!=5){
		printf("not enough arguments\n");
		return 1;
	}
	// laplace crossover variables:
	float a = 0.00;
	float b = 0.35;
	// BLX-alpha variable:
	float alphaX = 1;

	// example variables:
	float alpha = 0.366;
	float delta_i = 0.01;
	float delta_e = 0.00;

	/* initialize */
	// initialize new generation, uniform distribution within bounds
	int arg1 = atoi(argv[1]);	// crossover method
	int arg2 = atoi(argv[2]);	// mutation method
	int arg3 = atoi(argv[3]);	// selection method
	int arg4 = atoi(argv[4]);	// function choice

	for (int i=0;i<PPL_SIZE;i++){total_idx[i]=i;}
	srand((unsigned)time(NULL));
	for (int i=0;i<PPL_SIZE;i++){
		for (int j=0;j<VAR_SIZE;j++){
			n_gen[i][j] = randomgen(lowlimit,uplimit);
		}
		ngen_value[i] = feval(i,arg4);
	}
	// update fworst, old generation
	fworst = ngen_value[max(ngen_value,PPL_SIZE)];
	copygen();
	
	/* iterate */ 
	int iter = 0;
	float homogeneity = 0.0;
	//display(iter,homogeneity);
	float delta = delta_i;

	while (iter<ITER && homogeneity<0.9){
		for (int i=0;i<PPL_SIZE;i+=2){
			/* selection method */
			switch(arg3){
				case 0:
				updatetour();
				parent1 = min(ogen_value,TOUR_SIZE,tour_idx);
				updatetour();
				parent2 = min(ogen_value,TOUR_SIZE,tour_idx);
				break;
				case 1:
				parent1 = roulette();
				parent2 = roulette();
				//printarrayf(ogen_value,PPL_SIZE);
				//printf("%i,%i,%f,%f\n",parent1,parent2,ogen_value[parent1],ogen_value[parent2] );
				break;
				default:
				updatetour();
				parent1 = min(ogen_value,TOUR_SIZE,tour_idx);
				updatetour();
				parent2 = min(ogen_value,TOUR_SIZE,tour_idx);				
			}

			
			/* crossover method */
			switch(arg1){
				case 2: laplacex(a, b, i);
				break;
				case 1: BLX_alpha(alphaX,i);
				break;
				case 0: arithmcross(i,alpha);
				break;
				default : arithmcross(alphaX,i);
			}

			/* mutation methods */
			switch(arg2){
				case 0: 
				uniformutate(i,delta);
				uniformutate(i+1,delta);
				break;
				case 1:
				powermutate(i);
				powermutate(i+1);
				break;
				default:
				uniformutate(i,delta);
				uniformutate(i+1,delta);												
			}

			ngen_value[i] = feval(i,arg4);
			ngen_value[i+1] = feval(i+1,arg4);
		}
		fworst = ngen_value[max(ngen_value,PPL_SIZE)];
		copygen();
		iter ++;
		delta += (delta_e-delta_i)/ITER;

		homogeneity = homogeneitycal();
		
	}
	display(iter, homogeneity);
	return 0;
}

void display(int iter,  float homogeneity){
	// display current generation scores, min results
	int min_idx = min(ogen_value, PPL_SIZE, total_idx);
	//printarrayf(ngen_value,PPL_SIZE);
	printf("Iter[%i] Homogeneity is %f\n",iter,homogeneity);
	printf("min is %g, solution is: ",ogen_value[min_idx]);
	printgen(o_gen,min_idx);
	//printf("\n");
	return;	
}

void printgen(float arr[PPL_SIZE][VAR_SIZE],int idx){
	/* print individuals in the generation */
	for (int j=0;j<VAR_SIZE;j++){
		printf("%f, ", arr[idx][j]);
	}
	printf("\n");
}

void printarrayf(float *arr,int size){
	/* print a 1D array, sperated by periods */
	for (int i=0;i<size;i++){
		printf("%f, ", arr[i]);
	}
	printf("\n");
	return;
}

void printarrayi(int *arr,int size){
	/* print a 1D array, sperated by periods */
	for (int i=0;i<size;i++){
		printf("%i.", arr[i]);
	}
	printf("\n");
	return;
}

int max(float *arr, int size){
	/* simple max */
	int max = 0;
	for (int i=0;i<size;i++){
		if (arr[i]>arr[max]){
			max = i;
		}
	}
	return max;
}

float sum(float *arr, int size){
	/* sum up the given array */
	float sum = 0.0;
	for (int i=0;i<size;i++){
		sum += arr[i];
	}
	return sum;
}

int min(float *arr, int size, int *idx_arr){
	/* simple min */
	int min = 0;
	for (int i=0;i<size;i++){
		if (arr[idx_arr[i]]<arr[min]){
			min = idx_arr[i];
		}
	}
	return min;
}

int fit(int idx, int entry){
	/* check if constraints are satisfied */
	if (fabs(n_gen[idx][entry])>uplimit){
		return 0;
	}
	return 1;
}

float feval(int idx, int fchoice){
	/* evaluate function for indexed individual */
	float ans = 10.0*VAR_SIZE;
	if (fit(idx,0)&&fit(idx,1)&&fit(idx,2)){
		if (fchoice == 0){
			ans = powf(n_gen[idx][0],2.0)+powf(n_gen[idx][1],2.0)+powf(n_gen[idx][2],2.0);
		}
		else if (fchoice == 1){
			for (int i=0;i<VAR_SIZE;i++){
				ans += powf(n_gen[idx][i],2.0)-10.0*cos(2*PI*n_gen[idx][i]);
			}			
		}

	}
	else {
		ans = fworst + penval(idx);
	}
	//printf("answer = %f, ", ans);
	//printgen(n_gen,idx);

	return ans;
}

float penval(int idx){
	/* return penalty added to out of bound individuals */
	// penalty calculated 2 norm
	float penalty = 0.0;
	for (int i=0;i<VAR_SIZE;i++){
		if (!fit(idx,i)){
			penalty += powf(fabs(n_gen[idx][i])-uplimit,2);
		}
	}
	return penalty;
}

void copygen(void){
	/* copy new generation to old generation */
	for (int i=0;i<PPL_SIZE;i++){
		for (int j=0;j<VAR_SIZE;j++){
			o_gen[i][j] = n_gen[i][j];
		}
		ogen_value[i] = ngen_value[i];
	}
	return;
}

void updatetour(void){
	/* update indices in tournament */
	int fill = 0;
	int idx = 0;
	while (fill<TOUR_SIZE){
		if (randomgen(0,1)<TOUR_PERCENT){
			tour_idx[fill] = idx;
			fill ++;
		}
		idx ++;
		if (idx>=PPL_SIZE){
			idx = 0;
		}
	}
}

int roulette(void){
	/* find an index that corresponds to selected individual */
	float vmax = ogen_value[max(ogen_value,PPL_SIZE)];
	float vsum =  sum(ogen_value,PPL_SIZE);
	vsum = PPL_SIZE*vmax-vsum;
	float randomp = randomgen(0,1);
	float stackp = 0.0;
	int outindex;
	//printf("vsum is %f, randomp is %f\n", vsum,randomp);
	if (vsum == 0.0){
		return 0;
	}
	for (int i=0; i<PPL_SIZE; i++){
		stackp += float(vmax-ogen_value[i])/vsum;
		//printf("stackp is %f\n", stackp);
		if (stackp >= randomp){
			outindex = i;
			return outindex;
		}
	}
}

void laplacex(float a, float b, int idx){
	/*a & b are laplace distribution parameters, 
	idx is index for first offspring */
	float beta;
	if (randomgen(0,1)>CO_PROB){
		for (int j=0;j<VAR_SIZE;j++){
			n_gen[idx][j] = o_gen[parent1][j];
			n_gen[idx+1][j] = o_gen[parent2][j];
		}
	}
	else {
		for (int j=0;j<VAR_SIZE;j++){
			beta = a + powf(-1,rand()%2)*b*log(randomgen(0,1));
			n_gen[idx][j] = o_gen[parent1][j]+beta*fabs(o_gen[parent1][j]-o_gen[parent2][j]);
			n_gen[idx+1][j] = o_gen[parent2][j]+beta*fabs(o_gen[parent1][j]-o_gen[parent2][j]);		
		}
	}
	return;
}

void BLX_alpha(float alpha, int idx){
	/* use BLX crossover to update individuals at idx and idx+1 */
	float cmin, cmax, I;
	if (randomgen(0.0,1.0)<CO_PROB){
		for (int j=0;j<VAR_SIZE;j++){
			if (o_gen[parent1][j]>o_gen[parent2][j]){
				cmax = o_gen[parent1][j];
				cmin = o_gen[parent2][j];
			}
			else {
				cmax = o_gen[parent2][j];
				cmin = o_gen[parent1][j];			
			}
			I = cmax - cmin;
			n_gen[idx][j] = randomgen(cmin-I*alpha,cmax+I*alpha);
			n_gen[idx+1][j] = randomgen(cmin-I*alpha,cmax+I*alpha);		
		}
	}
	else {
		for (int j=0;j<VAR_SIZE;j++){
			n_gen[idx][j] = o_gen[parent1][j];
			n_gen[idx+1][j] = o_gen[parent2][j];
		}
	}
}

void powermutate(int idx){
	/* mutate indexed individual in new gen */
	for (int j=0;j<VAR_SIZE;j++){
		if (randomgen(0,1)<MU_PROB){
			float t =  (n_gen[idx][j]+5)/(5-n_gen[idx][j]);
			if (t<randomgen(0,1)){
				n_gen[idx][j] -= powf(randomgen(0,1),P)*(n_gen[idx][j]-lowlimit);
			}
			else {
				n_gen[idx][j] += powf(randomgen(0,1),P)*(uplimit-n_gen[idx][j]);
			}
		}
	}
}

float homogeneitycal(void){
	/* calculate the homogeneity of current population */
	int minimum = 0;
	float repeat = 0.0;
	minimum = min(ngen_value, PPL_SIZE, total_idx);
	for (int i=0;i<PPL_SIZE;i++){
		if (fabs(ngen_value[i]-ngen_value[minimum])<1e-8){
			repeat ++;
		}
	}
	return repeat/float(PPL_SIZE);
}

float randomgen(float minimum, float maximum){
	float ans = (double)rand()/RAND_MAX*(maximum-minimum);
	return ans+minimum;
}

void arithmcross(int idx, float alpha){
	if (randomgen(0.0,1.0)<CO_PROB){
		float r = randomgen(0.0,1.0)*(1+2.0*alpha)-alpha;
		for (int j=0;j<VAR_SIZE;j++){
			n_gen[idx][j] = r*o_gen[parent1][j]+(1-r)*o_gen[parent2][j];
		}
		for (int j=0;j<VAR_SIZE;j++){
			n_gen[idx+1][j] = r*o_gen[parent2][j]+(1-r)*o_gen[parent1][j];
		}
	}
	else {
		for (int j=0;j<VAR_SIZE;j++){
			n_gen[idx][j] = o_gen[parent1][j];
			n_gen[idx+1][j] = o_gen[parent2][j];
		}
	}
	return;
}

void uniformutate(int idx, float delta){
	/* uniform mutate */
	for (int j=0;j<VAR_SIZE;j++){
		if (randomgen(0,1)<MU_PROB){
			n_gen[idx][j] += randomgen(0.0,1.0)*delta*powf(-1,rand()%2);
		}
	}
	return;
}
