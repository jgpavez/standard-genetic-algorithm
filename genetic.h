#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef __GENETIC__
#define __GENETIC__

#define POPSIZE 50
#define MAXGENS 500
#define NVARS 2
#define PXOVER 0.25
#define PMUTATION 0.01
#define TRUE 1
#define FALSE 0
#define PI 3.14159

/*
 * struct genotype
 * binary representation 
 *
 **/

struct genotype{
	char *gene;
	double fitness;
	double upper[NVARS];
	double lower[NVARS];
	int precision[NVARS];
	int p[NVARS];
	double rfitness;
	double cfitness;
};

/*
 * Main Class
 * Standard Genetic Algorithm
 * To optimize real functions using binary representation
 * and standard operators
 * 
**/

class Genetic{
private:
	int GENESIZE;
	int generation;
	int cur_best;
	FILE *galog;
	char *filename;
	char *galog_file;
	int seed;
	struct genotype population[POPSIZE+1];
	struct genotype newpopulation[POPSIZE+1];
public:
	Genetic(char *file, char* g_file, int i_seed): filename(file), galog_file(g_file), seed(i_seed) {};
	void initialize(void);
	inline double randval(double, double);
	void evaluate(void);
	void keep_the_best(void);
	void elitist(void);
	void select_proportional(void);
	void select_proportional_wheel(void);
	void select_tournament(void);
	void select_ranking(void);
	void crossover(void);
	void Xover(int, int);
	inline void swap(char*, char*);
	void mutate(void);
	void report(void);
	void solve(void);
	inline void bit_encode(double, double, double, int, char *);
	inline void bit_decode(char *, double, double, int, double &);
	inline int compute_p(double, double, int);
};
#endif
