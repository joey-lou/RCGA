/* Solving Knapsack Problem using Genetic Algorithm */
/* Jiujiu Lou, 06/28/2018 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

const int PPL_SIZE = 100;					// population number
const int SIZE = 20;						// length of the sequences(number of packages)
const float CO_PROB = 0.85;					// crossover probability
const float MU_PROB = 0.05;					// mutation probability
const int ITER = 100;						// iteration(generation) number
const int SEED = 9;							// seed number to randomly generate samples
const float KEEP_PERCENT = 0.8;				// the percentange of whole population that gets updated
const int NUM_KEEP = PPL_SIZE*KEEP_PERCENT;	// number of individuals kept from generation to generation
const int NUM_UPDATE = PPL_SIZE-NUM_KEEP;	// number of individuals updated

/* global variables */
int population[PPL_SIZE][SIZE];				// 2D population representation
int ppl_scores[PPL_SIZE];					// scores of population
int ppl_weights[PPL_SIZE];					// weights of population
int parent1;								// index of parent1
int parent2;								// index of parent2
int total_index[PPL_SIZE];					// include all indices, may or may not be sorted
int keep_index[NUM_KEEP];					// indices of those we keep
int update_index[NUM_UPDATE];				// indices of those we discard/ update
int scores[SIZE];							// scores of all available packages
int weights[SIZE];							// weights of all available packages


/* list of helper functions */
int sum(int *arr, int size, int *index_arr); // take sum of given integer array entries 
int max(int *arr, int size);				// simple max finding
int weightcal(int i);						// calculate weight of ith sequence(individual)
int scorecal(int i);						// calculate score of ith sequence
void update(int i, int MAX_WEIGHT);			// update overweighted sequences by randomly deleting 1s in array
float homogeneitycal(void);					// find the percentage of max scored individuals in population
void sort_index(void);						// sort indices of total_index according to ppl_scores
int roulette(void);							// roulette selection method
void crossover(int index);					// perform crossover and replace the indexed individual with new sequence
void mutate(int index);						// mutate the indexed individual
void updateparent(void);					// update the indices to two parents using roullete
void updateindices(void);					// update the indices of all 3 index arrays
float randomgen(float min, float max);

/* visualization functions */
void printarray(int *arr,int n);			// print a 1D array, sperated by periods
void print_pplarray(int index);				// print the indexed individual in population (in binary)

int main(void){
	/* randomized sample data */
	// initialize scores 
	srand(SEED);
	for (int i=0;i<SIZE;i++){
		scores[i]=randomgen(0,100);
	}
	// initialize weights
	for (int i=0;i<SIZE;i++){
		weights[i]=randomgen(0,100);
	}
	//set the max weight
	int idx[SIZE];
	for (int i=0;i<SIZE;i++){
		idx[i] = i;
	}
	const int MAX_WEIGHT = sum(weights, SIZE, idx)/2;
	
	// update the total index to include all indices first
	for (int i=0;i<PPL_SIZE;i++){
		total_index[i]=i;
	}
	
	/*initialize population*/
	srand((unsigned)time(NULL));
	for (int i=0;i<PPL_SIZE;i++){
		for (int j=0;j<SIZE;j++){
			population[i][j] = rand()%2; 			// randomized start
		}
		/* update/calculate initial scores */
		ppl_weights[i] = weightcal(i);
		ppl_scores[i] = scorecal(i);
		while (ppl_weights[i] > MAX_WEIGHT){ 
			update(i,MAX_WEIGHT);  					// update overweighted sequences
		}
	}
	
	/* iterate */
	int iter = 0;
	float homogeneity = homogeneitycal();
	while (homogeneity<0.9 && iter<ITER){
		// update indices
		updateindices();
		for (int i=0;i<NUM_UPDATE;i++){				// populate new generation
			updateparent();
			crossover(update_index[i]);			
			mutate(update_index[i]);
			//print_pplarray(i);
			ppl_weights[update_index[i]] = weightcal(update_index[i]);
			ppl_scores[update_index[i]] = scorecal(update_index[i]);
			while (ppl_weights[update_index[i]] > MAX_WEIGHT){ 	// in case new individuals are overweighted
				update(update_index[i],MAX_WEIGHT);
				//printf("[%i] is updated\n", i);
			}
		}
		printarray(ppl_scores,PPL_SIZE);
		homogeneity = homogeneitycal();
		iter += 1;
		int total_score = sum(ppl_scores, PPL_SIZE, total_index);
		printf("Iteration# = %i. Homogeneity = %f. AveScore = %d\n\n",iter, homogeneity,total_score/PPL_SIZE);	
	}
	int index = max(ppl_scores,PPL_SIZE);
	printf("Best sequence is: ");
	print_pplarray(index);
	return 0;
}

void print_pplarray(int index){
	/* print the indexed individual in population (in binary) */
	for (int j=0;j<SIZE;j++){
		printf("%i",population[index][j]);
	}
	printf("\n");
	return;
}

void printarray(int *arr,int n){
	/* print a 1D array, sperated by periods */
	for (int i=0;i<n;i++){
		printf("%i.", arr[i]);
	}
	printf("\n");
	return;
}

int max(int *arr, int size){
	/* simple max */
	int max = 0;
	for (int i=0;i<size;i++){
		if (arr[i]>arr[max]){
			max = i;
		}
	}
	return max;
}

int sum(int *arr, int size, int *index_arr){
	/* sum up the given array */
	int sum = 0;
	for (int i=0;i<size;i++){
		sum += arr[index_arr[i]];
	}
	return sum;
}

int sequence_sum(int i){
	int sum = 0;
	for (int j=0;j<SIZE;j++){
		if (population[i][j] == 1){
			sum += 1;
		}
	}
	return sum;
}

int weightcal(int i){
	/* weight calculation for ith sequence*/
	int weight = 0;
	for (int j=0;j<SIZE;j++){
		if (population[i][j] == 1){
			weight += weights[j];
		}
	}
	return weight;
}

void update(int i, int MAX_WEIGHT){
	/* update individual i so that the weight is under constraint */
	float probility = 1.0/sequence_sum(i)*100;
	for (int j=0;j<SIZE;j++){
		if (population[i][j] == 1 && randomgen(0,1)<probility){
			population[i][j] = 0;
			if (weightcal(i)<MAX_WEIGHT){
				ppl_weights[i] = weightcal(i);
				ppl_scores[i] = scorecal(i);
				break;
			}
		}
	}
	return;
}1

int scorecal(int i){
	/* score calculation for ith sequence*/	
	int score = 0;
	for (int j=0;j<SIZE;j++){
		if (population[i][j] == 1){
			score += scores[j];
		}
	}
	return score;	
}

float homogeneitycal(void){
	/* calculate the homogeneity of current population */
	int max = 0;
	float repeat = 0;
	for (int i=0;i<PPL_SIZE;i++){
		if (ppl_scores[i]>max){
			repeat = 1;
			max = ppl_scores[i];
		}
		else if (ppl_scores[i]==max && max != 0){
			repeat += 1;
		}
	}
	return repeat/PPL_SIZE;
}

int roulette(void){
	/* find an index that corresponds to selected individual */
	float score_sum =  sum(ppl_scores,NUM_KEEP,keep_index);
	float stackp = 0.0;
	int outindex;
	for (int i=0; i<NUM_KEEP; i++){
		stackp += float(ppl_scores[keep_index[i]])/score_sum;
		if (stackp >= randomgen(0,1)){
			outindex = keep_index[i];
			//printf("[%i] is selected\n", outindex);
			return outindex;
		}
	}
	return 0;
}

void updateparent(void){
	/* update both parents using roulette method */
	parent1 = roulette();
	parent2 = roulette();
}

void sort_index(void){
	int temp;
  	for (int i=0; i<PPL_SIZE; i++){
    	for (int j=i+1; j<PPL_SIZE; j++){
      		if (ppl_scores[total_index[j]] > ppl_scores[total_index[i]]){
		        temp = total_index[i];
		        total_index[i] = total_index[j];
		        total_index[j] = temp;
      		}
    	}
  	}
  	return;
}

void updateindices(void){
	sort_index();
	for (int i=0;i<NUM_KEEP;i++){
		keep_index[i] = total_index[i];
	}
	for (int i=NUM_KEEP;i<PPL_SIZE;i++){
		update_index[i-NUM_KEEP] = total_index[i];
	}
}

void crossover(int index){
	/* replace indexed individual with crossovered sequence */
	int co_point = rand()%SIZE;
	if (randomgen(0,1)<CO_PROB){
		for (int j=0;j<co_point;j++){
			population[index][j] = population[parent1][j];
		}
		for (int j=co_point;j<SIZE;j++){
			population[index][j] = population[parent2][j];
		}
	}
	else {
		for (int j=0;j<SIZE;j++){
			population[index][j] = population[parent1][j];
		}
	}
	return;
}

void mutate(int index){
	/* mutate the indexed individual */
	float probility;
	for (int j=0;j<SIZE;j++){
			if (randomgen(0,1)<MU_PROB){
				population[index][j] = !population[index][j];
	}
	}
	return;
}

float randomgen(float min, float max){
	float ans = (float)rand()/((float)RAND_MAX/(max-min));
	return ans+min;
}
