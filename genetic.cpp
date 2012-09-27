#include "genetic.h"
#include <malloc.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;

double eval_function(double *x, int n_vars)
{
	return 0.5 + (pow(sin(sqrt(pow(x[1],2) + pow(x[2],2))),2)-0.5)/(1.0 + 0.0001*pow(pow(x[1],2) + pow(x[2],2),2));

}


void Genetic::initialize()
{
	ifstream infile("parameters.txt");
	int i,j,precision[NVARS], p[NVARS], actual_position;
	double lbound[NVARS],ubound[NVARS];;
	double random_value;

	if (!infile.is_open()){
		fprintf(galog,"\nCannot open input file!\n");
		exit(1);
	}

	srand(seed);

	for (i = 0; i < NVARS; i++){
		infile >> lbound[i];
		infile >> ubound[i];
		infile >> precision[i];
		p[i] = compute_p(lbound[i], ubound[i], precision[i]);
		GENESIZE += p[i];
	} 
	for (j = 0; j < POPSIZE ; j++){
		population[j].fitness = 0;
		population[j].rfitness = 0;
		population[j].cfitness = 0;
		population[j].gene = new char[GENESIZE];
		newpopulation[j].gene = new char[GENESIZE];
		for (i = 0;i < NVARS; i++){  
			population[j].lower[i] = lbound[i];
			population[j].upper[i] = ubound[i];
			population[j].precision[i] = precision[i];
			population[j].p[i] = p[i];
/*			random_value = randval(population[j].lower[i],population[j].upper[i]);
			char bit_string_temp[p[i]];
			bit_encode(random_value,lbound[i],ubound[i],precision[i],bit_string_temp);
			cout << bit_string_temp << endl;
			if ( i == 0 ) actual_position = 0;
			else actual_position += p[i-1];
			for ( int k = 0; k < p[i]; k++ )
				population[j].gene[k+actual_position-1] = bit_string_temp[k]; 
*/			
		}
		for (i = 0; i <= GENESIZE; i++){
			population[j].gene[i] = '0' + (int)random()%2; 
		}
	}
	population[POPSIZE].gene = new char[GENESIZE];
	newpopulation[POPSIZE].gene = new char[GENESIZE];
	infile.close();
}

int Genetic::compute_p(double lvalue, double uvalue, int precision)
{
	// TODO: Esto no funciona con intervalos menor a 1.
	int p;
	p = (int)ceil(log2(pow(10,precision)*(int)(uvalue-lvalue)));
	return p;
}

double Genetic::randval(double low, double high)
{
	double val;
	val = ((double)(rand()%1000)/1000.0) * (high - low) + low;
	return val;
}

void Genetic::evaluate()
{
	int mem,i;
	double x[NVARS+1];
	for (mem = 0; mem < POPSIZE; mem++){
		for (i = 0; i < NVARS; i++){
			if ( i == 0 )
				bit_decode(population[mem].gene,population[mem].lower[i],population[mem].upper[i],population[mem].precision[i],x[i+1]);
			else
				bit_decode(&population[mem].gene[population[mem].p[i-1]],population[mem].lower[i],population[mem].upper[i],population[mem].precision[i],x[i+1]);
		}
		/**
		* Poner aqui la funcion de evaluacion
		**/
		population[mem].fitness = eval_function(x, NVARS);

	}
}

void Genetic::keep_the_best()
{
	int mem,i,cur_best = 0;
	for (mem = 0; mem < POPSIZE; mem++){
		if (population[mem].fitness > population[POPSIZE].fitness){
		cur_best = mem;
		population[POPSIZE].fitness = population[mem].fitness;
		}
	}
	memcpy(population[POPSIZE].gene, population[cur_best].gene, sizeof(char)*GENESIZE);
}

void Genetic::elitist()
{
	int i;
	double best,worst;
	int best_mem, worst_mem;
	best = population[0].fitness; 
	worst = population[0].fitness;
	for (i = 0; i < POPSIZE - 1; i++){
		if (population[i].fitness > population[i+1].fitness){
			if (population[i].fitness >= best){
				best = population[i].fitness;
				best_mem = i;
			}
			if (population[i+1].fitness <= worst){
				worst = population[i+1].fitness;
				worst_mem = i+1;
			}
		}
		else{
			if (population[i].fitness <= worst){
				worst = population[i].fitness;
				worst_mem = i;
			}
			if (population[i+1].fitness >= best){
				best= population[i+1].fitness;
				best_mem = i+1;
			}
		}
	}
	if ( best >= population[POPSIZE].fitness){
		memcpy(population[POPSIZE].gene, population[best_mem].gene,sizeof(char)*GENESIZE);
		population[POPSIZE].fitness = population[best_mem].fitness;
	}
	else{
		memcpy(population[worst_mem].gene,population[POPSIZE].gene,sizeof(char)*GENESIZE);
		population[worst_mem].fitness = population[POPSIZE].fitness;
	}
}

void Genetic::select_proportional()
{
	int mem, i , j, k;
	double sum = 0, p, normalizer = 0.0;
	for ( mem = 0; mem < POPSIZE; mem++){
		if (population[mem].fitness < 0 && population[mem].fitness < normalizer)
			normalizer = population[mem].fitness;
	}
	for ( mem = 0; mem < POPSIZE; mem++){
		population[mem].fitness = population[mem].fitness + fabs(normalizer);
	}
	for ( mem = 0; mem < POPSIZE; mem++)
		sum += population[mem].fitness;
	for ( mem = 0; mem < POPSIZE; mem++)
		population[mem].rfitness = population[mem].fitness/sum;
	population[0].cfitness = population[0].rfitness;
	for ( mem = 1; mem < POPSIZE; mem++)
		population[mem].cfitness = population[mem-1].cfitness + population[mem].rfitness;
	for ( i = 0; i < POPSIZE; i++){
		p = rand()%1000/1000.0;
		if (p < population[0].cfitness){
			newpopulation[i] = population[0];
			memcpy(newpopulation[i].gene,population[0].gene,sizeof(char)*GENESIZE); 

		}
		else
		{
			for ( j = 0; j < POPSIZE; j++)
				if ( p >= population[j].cfitness && p < population[j+1].cfitness){
					newpopulation[i] = population[j+1];
					memcpy(newpopulation[i].gene, population[j+1].gene,sizeof(char)*GENESIZE);

			}
		}
	}
	for ( i = 0; i < POPSIZE; i++ ){
		population[i] = newpopulation[i];
		memcpy(population[i].gene, newpopulation[i].gene,sizeof(char)*GENESIZE);
	}
}

void Genetic::select_proportional_wheel()
{
	int mem, i , j, k,index = 0;
	double max_rfitness = 0.0;
	double sum = 0, beta = 0.0, normalizer = 0.0;
	for ( mem = 0; mem < POPSIZE; mem++){
		if (population[mem].fitness < 0 && population[mem].fitness < normalizer)
			normalizer = population[mem].fitness;
	}
	for ( mem = 0; mem < POPSIZE; mem++){
		population[mem].fitness = population[mem].fitness + fabs(normalizer);
	}
	for ( mem = 0; mem < POPSIZE; mem++)
		sum += population[mem].fitness;
	for ( mem = 0; mem < POPSIZE; mem++){
		population[mem].rfitness = population[mem].fitness/sum;
		if ( population[mem].rfitness > max_rfitness )
			max_rfitness = population[mem].rfitness;
	}
	index = (int) randval(0.0, (double)POPSIZE);
	for ( i = 0; i < POPSIZE; i++){
		beta = beta + randval(0.0, max_rfitness);
		while (population[index].rfitness < beta){
			beta = beta - population[index].rfitness;
			index += 1;
		}
		newpopulation[i] = population[index];
		memcpy(newpopulation[i].gene, population[index].gene,sizeof(char)*GENESIZE);

	}
	for ( i = 0; i < POPSIZE; i++ ){
		population[i] = newpopulation[i];
		memcpy(population[i].gene, newpopulation[i].gene,sizeof(char)*GENESIZE);
	}


}

void Genetic::select_tournament()
{
	int i, s1, s2;

	for (i = 0; i < POPSIZE; i++){
		s1 = (int)randval(0.0, (double)POPSIZE);
		do 
			s2 = (int)randval(0.0, (double)POPSIZE);
		while(s1 == s2);
		if (population[s1].fitness > population[s2].fitness)
		{
			newpopulation[i] = population[s1];
			memcpy(newpopulation[i].gene,population[s1].gene,sizeof(char)*GENESIZE);		
		}
		else{
			newpopulation[i] = population[s2];
			memcpy(newpopulation[i].gene,population[s2].gene,sizeof(char)*GENESIZE);		
		}
	}
	for ( i = 0; i < POPSIZE; i++ ){
		population[i] = newpopulation[i];
		memcpy(population[i].gene, newpopulation[i].gene,sizeof(char)*GENESIZE);
	}
}

int compare(const void *arg1, const void *arg2)
{
	struct genotype *ind1;
	struct genotype *ind2;
	ind1 = (struct genotype*)arg1;
	ind2 = (struct genotype*)arg2;
	if (ind1->fitness > ind2->fitness)
		return -1;
	else if(ind1->fitness > ind2->fitness)
		return 1;
	else
		return 0;
}

void Genetic::select_ranking()
{
	int mem, i , j, k;
	double sum = 0, p, normalizer = 0.0;
	for ( mem = 0; mem < POPSIZE; mem++){
		if (population[mem].fitness < 0 && population[mem].fitness < normalizer)
			normalizer = population[mem].fitness;
	}
	for ( mem = 0; mem < POPSIZE; mem++){
		population[mem].fitness = population[mem].fitness + fabs(normalizer);
	}
	qsort(population, POPSIZE, sizeof(struct genotype), compare); 
	for ( mem = 0; mem < POPSIZE; mem++)
		sum += i+1;
	for ( mem = 0; mem < POPSIZE; mem++)
		population[mem].rfitness = (double)(i+1)/sum;
	population[0].cfitness = population[0].rfitness;
	for ( mem = 1; mem < POPSIZE; mem++)
		population[mem].cfitness = population[mem-1].cfitness + population[mem].rfitness;
	for ( i = 0; i < POPSIZE; i++){
		p = rand()%1000/1000.0;
		if (p < population[0].cfitness){
			newpopulation[i] = population[0];
			memcpy(newpopulation[i].gene,population[0].gene,sizeof(char)*GENESIZE); 
		}
		else
		{
			for ( j = 0; j < POPSIZE; j++)
				if ( p >= population[j].cfitness && p < population[j+1].cfitness){
					newpopulation[i] = population[j+1];
					memcpy(newpopulation[i].gene, population[j+1].gene,sizeof(char)*GENESIZE);

			}
		}
	}
	for ( i = 0; i < POPSIZE; i++ ){
		population[i] = newpopulation[i];
		memcpy(population[i].gene, newpopulation[i].gene,sizeof(char)*GENESIZE);
	}
}

void Genetic::crossover()
{
	int i, mem, one, first = 0;
	double x;

	for ( mem = 0; mem < POPSIZE; ++mem){
		x = rand()%1000/1000.0;
		if ( x < PXOVER){
			++first;
			if (first % 2 == 0)
				Xover(one, mem);
			else 
				one = mem;
		}
	}
}

void Genetic::Xover(int one, int two)
{
	int i,point;
	if (GENESIZE > 1) {
		if (GENESIZE == 2)
			point = 1;
		else
			point = (rand() % (GENESIZE - 1)) + 1;
		for ( i = 0; i < point; i++ )
			swap(&population[one].gene[i], &population[two].gene[i]);
	
	}
}

void Genetic::swap(char *x, char *y)
{
	char temp;
	temp = *x;
	*x = *y;
	*y = temp;

}

void Genetic::mutate()
{
	int i,j;
	double lbound, hbound;
	double x;
	for (i = 0; i < POPSIZE; i++){
		for (j = 0; j < GENESIZE; j++){
			x = rand()%1000/1000.0;
			if ( x > PMUTATION){
				population[i].gene[j] = '0' + (int)random()%2;			}
		}
	}
}

void Genetic::report()
{
	int i;
	double best_val, avg, stddev, sum_square, square_sum, sum;
	sum = 0.0; sum_square = 0.0;
	for ( i = 0; i < POPSIZE; i++){
		sum += population[i].fitness;
		sum_square += population[i].fitness * population[i].fitness;
	}
	avg = sum/(double) POPSIZE;
	square_sum = avg * avg * POPSIZE;
	stddev = sqrt((sum_square - square_sum)/(POPSIZE - 1));
	best_val = population[POPSIZE].fitness;
	fprintf(galog, "%5d,	%6.3f, %6.3f, %6.3f \n", generation, best_val, avg, stddev);
}

void Genetic::bit_encode(double value, double lvalue, double uvalue, int precision, char *bit_string)
{
	int x_10,p,counter = 0;
	// TODO: Esto no funciona con intervalos menor a 1.
	p = (int)ceil(log2(pow(10,precision)*(int)(uvalue-lvalue)));
	char *bit_string_temp = new char[p];
	x_10 = (int)((value - lvalue)*(pow(2,p)/(int)(uvalue-lvalue)));
	
	for ( counter = 0; counter < p; counter++ ){
		if ( x_10 > 0){
			bit_string_temp[counter] = '0' + (int)fmod(x_10,2);
			x_10 = (int)floor(x_10/2);
		}
		else{
			bit_string_temp[counter] = '0';
		}
	}
	for ( counter  = 0; counter < p ; counter++){
		bit_string[counter] = bit_string_temp[p - counter - 1];
	}
	delete bit_string_temp;
}

void Genetic::bit_decode(char *bit_string, double lvalue, double uvalue, int precision, double &value)
{
	int x_10 = 0,p,counter = 0;
	// TODO: Esto no funciona con intervalos menor a 1.
	p = (int)ceil(log2(pow(10,precision)*(int)(uvalue-lvalue)));
	for ( int counter = p-1; counter >= 0; counter --){
		x_10 += (int)(((int)bit_string[counter]-48)*pow(2,p-counter-1));
	}
	 
	value = lvalue + (double)x_10 * ((double)(uvalue-lvalue)/(double)(pow(2,p)-1));
}

void Genetic::solve()
{
	int i;
	if ((galog = fopen(galog_file, "w")) == NULL) exit(1);
	generation = 0;
	
	fprintf(galog, "\n generation number, best value, average fitness, standard deviation \n");
//	fprintf(galog, " number value fitness deviation \n");
	
	initialize();

	evaluate();
	keep_the_best();

//	for ( i = 0; i < POPSIZE; i++)
//		cout<<population[i].gene<<endl;


	while(generation < MAXGENS){
		generation++;
		select_proportional_wheel();
		crossover();	
		mutate();
		report();
		evaluate();
		elitist();
	}	

//	printf("Simulation completed\n");
//	printf("Best member: \n");
	double x[NVARS];
	printf("%d, ", seed);
	for (i = 0; i < NVARS; i++){
		if ( i == 0 )
			bit_decode(population[POPSIZE].gene,population[POPSIZE-1].lower[i],population[POPSIZE-1].upper[i],population[POPSIZE-1].precision[i],x[i]);
		else
			bit_decode(&population[POPSIZE].gene[population[POPSIZE-1].p[i-1]],population[POPSIZE-1].lower[i],population[POPSIZE-1].upper[i],population[POPSIZE-1].precision[i],x[i]);
	}
	for (i = 0; i < NVARS; i++){
		printf("%3.3f, ", i, x[i]);
	}
	printf("%3.3f\n", population[POPSIZE].fitness);
	fclose(galog);
}
