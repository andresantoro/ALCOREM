#include <search.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include <math.h>
#include <zlib.h>
#include <assert.h>
#include <omp.h>

// These variables are required in the computation of the Complexity 
unsigned long long int length; 
unsigned long long int prime_layer;
unsigned int *maximum_overlap;
mpz_t *maximum_prime_product;
unsigned long long int *Nlinks_aggregate;
unsigned long long int *Nlinks_multiplex;
char **fp_string;
char **fp_string_aggregate;
unsigned long long int *dimension_string;
unsigned long long int *string_counter;
unsigned long long int *string_counter_aggregate;
double **omega;
double *complexity_multiplex_reducibility;

//Definition of the structures
struct edge{
	unsigned long long int node_i;
	unsigned long long int node_j;
	unsigned long long int key;
	unsigned int overlap_aggregate;
	mpz_t weight_link;
};

struct layer{
	struct edge *links_table;
	unsigned long long int max_node;
	unsigned long long int nlinks;
	int l1;
	int l2;
};


struct tree {
    // datum must be the first field in struct tree
    const void *datum;
    struct tree *left, *right;
};

typedef void freer(void *node);

void usage(char *argv[]){
  printf(
	"******************************************************************************\n"
	"**                                                                          **\n"
	"**         Reducibility algorithm for a multiplex system based on           **\n"
	"**                          information complexity                          **\n"
	"**                                                                          **\n"
	"**********+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*********\n"
	"**                                                                          **\n"
	"**  <unique_nodes> list containing the list of node IDs in the multiplex    **\n"    
	"**                                                                          **\n"
	"**                                                                          **\n"
	"**           <filein> should be a list of the layer composing the           **\n"
	"**                            multiplex network                             **\n"
	"**                                                                          **\n"
	"**                                                                          **\n"
	"**             ----------   Optional Variables  ---------                   **\n"
	"**                                                                          **\n"
	"**           <ITER> represents the number of times the Kolmogorov           **\n"
	"**                    complexity is averaged (default = 1000)               **\n"
    "**                                                                          **\n"
	"**    <#core> represents the number of cores used in the computation of     **\n"
	"**     the complexity (default: sequential computation, i.e. #core =1)      **\n"
	"**                                                                          **\n"
	"**               <-v>  be more verbose -- show misc extra info              **\n"
	"**                                                                          **\n"
	"**         <-s>  save optimal multiplex on files in current directory       **\n"
	"**                                                                          **\n"
	"**  <-w>  save all the multiplex aggregation on files in current directory  **\n"
	"**                                                                          **\n"
	"**      OUTPUT: by default the algorithm returns the following info:        **\n"
    "**     partition of the multiplex --> Number of layers, KC multiplex,       **\n"
    "**    # links multiplex, KC aggregate, # links aggregate, complexity, q     **\n"
	"**                                                                          **\n" 
	"******************************************************************************\n");
  printf("Usage: %s <unique_nodes> <filein>  (ITER = 100) (-p #core) (-v verbose) (-s save optimal multiplex on files) (-w save all the reduced multiplex on files) [fileout]\n\n" , argv[0]);
}

void chopN_string(char *, size_t);

unsigned long long int *load_unique_nodes(char *, unsigned long long int *, unsigned long long int *,unsigned int );

struct layer *load_data(char *listoffiles, unsigned int *, unsigned long long int *, unsigned long long int ,unsigned int );

unsigned long long int load_edgelist(char *, struct edge **, unsigned long long int *, unsigned long long int );

int cmpfunc (const void *, const void *);

void list_prime(unsigned long long int **,unsigned int );

void allocate_memory(struct layer *,unsigned int );

void reduce_multiplex(void **,struct layer *, unsigned int, unsigned int, unsigned long long int,
                      unsigned long long int *, unsigned long long int *, unsigned long long int , int ,
                      double ***, int **, unsigned int ,unsigned int ,unsigned int , unsigned int );

void evaluate_KC(struct layer *, unsigned int ,unsigned long long int, unsigned long long int *,
                 unsigned long long int *, unsigned long long int, unsigned int,unsigned int);

struct layer *multiplex_into_binarytree(void **,struct layer *, unsigned int, unsigned long long int,
                                        unsigned long long int *, unsigned long long int *, unsigned long long int,
                                        unsigned long long int *, unsigned int );

unsigned long long int *particular_primes_sequence(struct layer *, unsigned int);

unsigned long long int *shuffle_nodes(unsigned long long int, unsigned long long int *);

struct layer *update_nodelabel(struct layer *, unsigned long long int *, unsigned int,unsigned long long int, unsigned long long int *);

void free_structgmp_multiplex(struct layer *, struct layer *, unsigned int);


void tdestroy(void *root);

int delete_root(const void *node1, const void *node2);

int compare_keys(const void *, const void *);

int mpz_cmp_flag_to_01(int );

unsigned long long int print_file_bothfast_onstring(void **,unsigned long long int *, unsigned int);

void print_elements_onstring_weighted_aggregate(const void *, const VISIT, const int);

unsigned long long int  dynamic_compress(char*, unsigned long long int);

void aggregate(struct layer **,unsigned int, unsigned int, unsigned int, unsigned long long int,unsigned int);

void recursive_print (struct layer *, unsigned int,unsigned int , unsigned int, unsigned int,
                      int **,unsigned int,unsigned int);

void free_structure(unsigned long long int *,unsigned long long int *,struct layer *, unsigned int);

int * compute_sequence_layers(struct layer *, unsigned int,unsigned int );

void dump_multiplex_on_file(struct layer *multiplex,int *sequence_multiplex_toprint,
														unsigned int numberoflayers_maximum_reducibility, unsigned long long int *nodes_unique, unsigned long long int);

void dump_all_multiplex_on_file(struct layer *multiplex,int *sequence_multiplex_toprint,
														unsigned int numberoflayers_maximum_reducibility, unsigned long long int *nodes_unique, unsigned long long int);


int main(int argc, char *argv[])
{

	unsigned long long int  *primeslist,*nodes_unique,length_unique, MAX_NODE,N;
	unsigned int M,i,j,k,ITER,parallel=1,verbose=0,partition_verbose=1,save_multiplex=0,save_multiplex_total=0;
	double maximum_reducibility=0;
	unsigned int index_maximum_reducibility,numberoflayers_maximum_reducibility;
	int counter;
	int *sequence_multiplex_toprint;
  	struct layer *multiplex;
  	unsigned long int randval,esc;
  	void *root;
  	char *str_par;
	FILE *f;
	ITER = 1000;

	 // RANDOM SEED
	srand(time(NULL));
	// OR
	// f = fopen("/dev/urandom", "r");
	// esc=fread(&randval, sizeof(randval),1, f);
	// fclose(f);
	// printf("%lu\n", randval);
	// srand(randval);


	if (argc <= 2){
    usage(argv);
    exit(1);
  	}
  	if (argc > 3)
  	{
  		if (argv[3][0] !='-')
  		{
  			ITER = (unsigned int)atoi(argv[3]);
			fprintf(stderr,"Number of iterations: %u\n", ITER);
		}
  }

	//Checking for the parallel option
	for (i = 1; i < argc; i++) 
	{
 		if (argv[i][0] == '-'  && argv[i][1] == 'p') 
 		{
 			if (argc > i+1)
 				parallel = (unsigned int)atoi(argv[i+1]);
 			else
 			{
 				str_par=argv[i];
 				chopN_string(str_par, 2);
 				parallel = (unsigned int)atoi(argv[i]);
 			}
 			fprintf(stderr,"Number of core inserted: %u\n",parallel);
 		}
 		if (argv[i][0] == '-' && argv[i][1] == 'v') 
 		{
 			fprintf(stderr,"Verbose activated\n");
 			verbose = 1;
 		}

 		if (argv[i][0] == '-' && argv[i][1] == 's') 
 		{
 			fprintf(stderr,"Save multiplex activated\n");
 			save_multiplex = 1;
 		}

 		if (argv[i][0] == '-' && argv[i][1] == 'w') 
 		{
 			fprintf(stderr,"Write  all multiplex on file activated\n");
 			save_multiplex_total = 1;
 		}
 	}	

 	if(verbose)
	fprintf(stderr,"Finding the number of nodes...\n");
	nodes_unique=load_unique_nodes(argv[1],&length_unique, &MAX_NODE,verbose);
	multiplex=load_data(argv[2],&M,nodes_unique,length_unique,verbose);
	if (M==1)
	{
		fprintf(stderr,"The multiplex inserted is only made of 1 layers! Exit now!\n\n");
		exit(1);
	}
	if(verbose)
	fprintf(stderr,"Computing the prime list...\n");
	primeslist=(unsigned long long int *)malloc(M*sizeof(unsigned  long long int));
	list_prime(&primeslist,M);
	//****************************************************
	//Number of M primes used in the script
	// for (i=0;i<M;i++)
	// {
	// 	fprintf(stderr,"%u\n", primeslist[i]);
	// }
	//***************************************************
	allocate_memory(multiplex,M);
	double **distance_matrix;
	int *vectorprint;
	distance_matrix=(double **)calloc((2*M-1),sizeof(double *));
	omega=(double **)calloc((2*M-1),sizeof(double *));
		for (i = 0; i < (2*M-1); i++)
		{
			distance_matrix[i]=(double *)calloc((2*M-1),sizeof(double));
			omega[i]=(double *)malloc((2*M-1)*sizeof(double));
			multiplex[i].l1=i;
			multiplex[i].l2=i;
		}

	complexity_multiplex_reducibility=(double *) malloc(M*sizeof(double));

	reduce_multiplex(&root,multiplex,M,M,MAX_NODE,primeslist,nodes_unique,length_unique,0,&distance_matrix,&vectorprint,ITER,parallel,verbose,partition_verbose);


	
	for (i=0; i<M; i++)
	{
		if (complexity_multiplex_reducibility[i] > maximum_reducibility )
		{
			maximum_reducibility=complexity_multiplex_reducibility[i];
			index_maximum_reducibility = i;
			numberoflayers_maximum_reducibility=i+1;
		}
			//fprintf(stderr,"%d %lf\n",i+1,complexity_multiplex_reducibility[i]);
	}
	if (verbose == 1)
		fprintf(stderr,"\nM_optimal: %d\t\tMaximum quality function: %lf\n",numberoflayers_maximum_reducibility,maximum_reducibility);

	if( save_multiplex == 1)
	{
		sequence_multiplex_toprint=compute_sequence_layers(multiplex,index_maximum_reducibility,M);
		dump_multiplex_on_file(multiplex,sequence_multiplex_toprint,numberoflayers_maximum_reducibility,nodes_unique,length_unique);
	}

	if( save_multiplex_total == 1)
	{
		for (counter = M-2; counter>=0;counter--)
		{
			index_maximum_reducibility=counter;
			numberoflayers_maximum_reducibility=counter+1;
			sequence_multiplex_toprint=compute_sequence_layers(multiplex,index_maximum_reducibility,M);
			dump_all_multiplex_on_file(multiplex,sequence_multiplex_toprint,numberoflayers_maximum_reducibility,nodes_unique,length_unique);
		}
	}

	//Free process
	free_structure(primeslist,nodes_unique,multiplex,2*M-1);
	for (i = 0; i < (2*M-1); i++)
	{
		free(distance_matrix[i]);
		free(omega[i]);
	}
	free(distance_matrix);
	free(vectorprint);
	free(complexity_multiplex_reducibility);
	free(omega);

}


 void chopN_string(char *str, size_t n)
{
    assert(n != 0 && str != 0);
    size_t len = strlen(str);
    if (n > len)
        return; 
    memmove(str, str+n, len - n + 1);
}


unsigned long long int *load_unique_nodes(char *filename, unsigned long long int *length_unique, unsigned long long int *MAX_NODE,unsigned int verbose)
{
  FILE *fp;
  unsigned long long int *unique_nodes, size=10000, k=0;
  int scan=0;
  fp = fopen(filename,"r");
  unique_nodes=(unsigned long long int *)malloc(size*sizeof(unsigned long long  int));
  while(scan!=EOF)
  {
  	if(k==size)
    {
      size+=10000;
      unique_nodes=(unsigned long long int *)realloc(unique_nodes,size*sizeof(unsigned long long int));
    }
    scan=fscanf(fp,"%llu",&(unique_nodes[k]));
    k++;
  }
  k=k-1;
  unique_nodes=(unsigned long long int *)realloc(unique_nodes,k*sizeof(unsigned long long int));
  if(verbose)
  fprintf(stderr,"Number of unique nodes: %llu\n\n",k);
  fclose(fp);
  *MAX_NODE = unique_nodes[k-1];
  *length_unique=k;
  return(unique_nodes);
}

struct layer *load_data(char *listoffiles, unsigned int *M, unsigned long long int *unique_nodes, unsigned long long int length_unique,unsigned int verbose)
{
	FILE *fp,*singlefp;
	struct layer *array;
	char line[256];
	char* token;
	unsigned long long int N,i,size;
	size = 10000;
	array = (struct layer *)malloc(size * sizeof(struct layer ));

	i=0;
	fp= fopen(listoffiles,"r");
	while (1) 
	{
		if (i == size)
		{
			size+=10000;
      array= (struct layer *)realloc(array,size*sizeof(struct layer));
		}


    if (fgets(line,256, fp) == NULL) break;
		token = strtok(line,"\n");
  	N=unique_nodes[length_unique-1];
		array[i].max_node=N;
		array[i].nlinks=load_edgelist(token,&(array[i].links_table),unique_nodes,length_unique);
		if(verbose)
		fprintf(stderr,"%s --> Number of links: %llu\n\n",token,array[i].nlinks);
		i++;
	}
	fclose(fp);
	*M=i;
	array= (struct layer *)realloc(array,(2*(*M)-1)*sizeof(struct layer));
	if(verbose)
	fprintf(stderr,"Number of layers in the multiplex:%llu\n\n",i);

	return(array);
}


unsigned long long int load_edgelist(char *filename, struct edge **current, unsigned long long int *unique_nodes, unsigned long long int length_unique)
{
	unsigned long long int temp1,temp2,counter_edges=0, size =1000,Nedges,auxiliary;
	unsigned long int scan=0;
	unsigned long long int *item;
	char buff[256];
	char *ptr;
	FILE *fp;

	*(current)=(struct edge *) malloc(size * sizeof(struct edge ));
	unsigned long long int *p_a=&(unique_nodes[0]);
  	unsigned long long int *p_b;
 	unsigned long long int differenceInBytes;

	fp= fopen(filename,"r");
	while(scan!=EOF)
	{
		if (counter_edges == size)
		{
			size+=1000;
      		*(current)=(struct edge*)realloc(*(current),size*sizeof(struct edge ));
   		}
	// PARSING DATA WITH NODES_UNIQUE -> Take the corresponding index -1 of that element
	scan=fscanf(fp, "%llu", &temp1);
    item = (unsigned long long int*) bsearch (&temp1, unique_nodes, length_unique, sizeof (unsigned long long int), cmpfunc);
    p_b=item;
    differenceInBytes = (p_b - p_a);
    temp1=differenceInBytes;
		
	scan=fscanf(fp, "%llu", &temp2);	
    item = (unsigned long long int*) bsearch (&temp2, unique_nodes, length_unique, sizeof (unsigned long long int), cmpfunc);
    p_b=item;
    differenceInBytes = (p_b - p_a);
    temp2=differenceInBytes;
    //The link that we are going to add is formed by the nodes temp1 --- temp2
    //In order to maintain the first element smaller than the second do this check
		if(temp1>temp2)
		{
			auxiliary=temp2;
			temp2=temp1;
			temp1=auxiliary;
		}
		
		(*current)[counter_edges].node_i = temp1;
		(*current)[counter_edges].node_j = temp2;
		//fprintf(stderr,"%u  %u\n\n",temp1, temp2);
		counter_edges++;
	}
	Nedges=counter_edges-1;
	*(current)=(struct edge*)realloc(*(current),Nedges*sizeof(struct edge ));
	fclose(fp);
	return(Nedges);
}

int cmpfunc (const void * a, const void * b)
{
   return ( *(long long int*)a - *(long long int*)b );
}

void list_prime(unsigned long long int **listofprimes,unsigned int numberofprimes)
{

	unsigned long long int i,s;
	mpz_t prime;
 	mpz_init(prime); //Initialize prime and set its value to 0.
 	mpz_set_ui(prime, 2); //Set the value of prime to be 2
	(*listofprimes)[0]=2;
    for(i=1; i<numberofprimes; i++) 
    {
  		mpz_nextprime(prime, prime); //Set prime to be the next prime number greater than prime.
    	s=mpz_get_ui(prime);
    	(*listofprimes)[i]=s;
	}
	mpz_clear(prime);
}


void allocate_memory(struct layer *multiplex,unsigned int M)
{
	long long int i,k;
	for (k = 0; k < M; k++)
	{
		for (i = 0; i < multiplex[k].nlinks; i++)
		{
			mpz_init(multiplex[k].links_table[i].weight_link); //Initialize the link weight and set its value to 0.
			multiplex[k].links_table[i].overlap_aggregate=0;
		}
	}
}




void reduce_multiplex(void **root,struct layer *multiplex, unsigned int M, unsigned int  Mmod, unsigned long long int MAX_NODE,
 unsigned long long int *primeslist, unsigned long long int *nodes_unique, unsigned long long int length_unique, int flag,
  double ***distance_matrix, int **vectorprint, unsigned int ITER,unsigned int parallel,unsigned int verbose, unsigned int partition_verbose)
{
	unsigned long long int dimension, dimension_KC, dimension_KC_aggregate, step, i,j,k,nedges_intersection;
	unsigned int flag_intern,omegastar,nthread,*layers_order,row, column,new_row, new_column;
	double sumKC,sumKCaggregate,max_KC_step=0,fraction,min,cumulative_complexity=0.0;
	struct layer **temporarystruct;
	void **nthread_root;

	struct layer *duplex;
	struct layer *new_multiplex;
	
	*root = NULL;
	//first cycle 
	if (flag == 0)
	{
		//***********Evaluate the Kolmogorov Complexity of the whole original multiplex with M layers************
		if(partition_verbose || verbose)
		{
			for (i=0;i<M;i++)
				printf("%llu,  ",i);
			printf("-->");
		}
		evaluate_KC(multiplex,M,MAX_NODE,primeslist,nodes_unique,length_unique,ITER,parallel);
		//********************************************************************************
		/*Creation of the duplex structure for evaluating the similarity for all the (M 2) couples of layers*/
		duplex=(struct layer *)malloc(2*sizeof(struct layer));
		*vectorprint=(int *)calloc ((2*M-1),sizeof(int));
		for (i = 0; i< M-1; i++)
		{
			for (j = i+1; j< M; j++)
			{
				*root = NULL;

				duplex[0].nlinks=multiplex[i].nlinks;
				duplex[0].links_table=multiplex[i].links_table;
				duplex[0].max_node=multiplex[i].max_node;

				duplex[1].nlinks=multiplex[j].nlinks;
				duplex[1].links_table=multiplex[j].links_table;
				duplex[1].max_node=multiplex[j].max_node;

				sumKC=0;
				sumKCaggregate=0;
				nedges_intersection=0;
				cumulative_complexity=0.0;
				omp_set_num_threads(parallel);
				//Allocation for the parallel computing
				temporarystruct = (struct layer **) malloc (parallel*sizeof(struct layer *));
				nthread_root = (void **) malloc (parallel*sizeof(void *));
				Nlinks_aggregate= (unsigned long long int *)malloc(parallel*sizeof(unsigned long long int ));
				Nlinks_multiplex= (unsigned long long int *)malloc(parallel*sizeof(unsigned long long int ));
				dimension_string= (unsigned long long int *)malloc(parallel*sizeof(unsigned long long int ));
				string_counter= (unsigned long long int *)malloc(parallel*sizeof(unsigned long long int ));
				string_counter_aggregate= (unsigned long long int *)malloc(parallel*sizeof(unsigned long long int ));
				fp_string= (char **)malloc(parallel*sizeof(char *));
				fp_string_aggregate= (char **)malloc(parallel*sizeof(char *));
				maximum_overlap=(unsigned int *)malloc(parallel*sizeof(unsigned int ));
				maximum_prime_product=(mpz_t *)malloc(parallel*sizeof(mpz_t));
				for (k=0;k<parallel;k++)
					mpz_init(maximum_prime_product[k]);
				#pragma omp parallel for private (nthread)
				for (k=0;k<ITER;k++)
				{
					nthread= omp_get_thread_num();
					nthread_root[nthread]=NULL;
					length=0;
					Nlinks_aggregate[nthread]=0;
					Nlinks_multiplex[nthread]=0;
					prime_layer=primeslist[1];

					temporarystruct[nthread]=multiplex_into_binarytree(&(nthread_root[nthread]),duplex,2,MAX_NODE,primeslist,nodes_unique,length_unique,&nedges_intersection,nthread);
					if (maximum_overlap[nthread] == 0)
						maximum_overlap[nthread] = 1;
					//Printing the duplex  and the corresponding aggregate on strings
					dimension_KC=print_file_bothfast_onstring(nthread_root[nthread],&dimension_KC_aggregate,nthread);
					sumKC+=dimension_KC;
					sumKCaggregate+=dimension_KC_aggregate;
					cumulative_complexity+=(1.0*dimension_KC)/(1.0*dimension_KC_aggregate);
					tdestroy(nthread_root[nthread]);
					free_structgmp_multiplex(temporarystruct[nthread],multiplex,2);
				}
				sumKC=(1.0*sumKC)/(1.0*ITER);
				sumKCaggregate=(1.0*sumKCaggregate)/(1.0*ITER);	
				fraction=(1.0*cumulative_complexity)/(1.0*ITER);


				(*distance_matrix)[i][j]=fraction;
				(*distance_matrix)[j][i]=fraction;
				omega[i][j]=((1.0*nedges_intersection)/(1.0*ITER))/(1.0*length);
				omega[j][i]=omega[i][j];

				//Find the duplex (couple of layers) \bar{D} with the highest value of complexity, this is the reference point
				if( fraction > max_KC_step)
				{
					max_KC_step = fraction;
					row=i;
					column=j;
				}

			free(temporarystruct);
			free(nthread_root);
			free(Nlinks_aggregate);
			free(Nlinks_multiplex);
			free(dimension_string);
			free(string_counter);
			free(string_counter_aggregate);
			free(fp_string);
			free(fp_string_aggregate);
			free(maximum_overlap);
			for (k=0;k<parallel;k++)
				mpz_clear(maximum_prime_product[k]);
			free(maximum_prime_product);
			}
		}
		free(duplex);
		//*********************************************************************************************************************************
		// Finding the couple of layers D with the smallest value of complexity  (among the set of feasible candidates that have a value of 
		// overlap greater or equal than the overlap of \bar(D))
		min = 1000.0;
		for(i=0;i<M-1;i++)
			for (j = i+1; j < M; j++)
			{
				if( (*distance_matrix)[i][j] <= max_KC_step && omega[i][j] >= omega[row][column]  &&   (*distance_matrix)[i][j]  < min)
				{
					min =(*distance_matrix)[i][j] ;
					new_row=i;
					new_column=j;
				}
			}

		row=new_row;
		column=new_column;
		if (verbose)
			fprintf(stderr,"\n\nThe index with the minimum value %lf is %u %u\n",min,row, column);
		for (i=0; i<(2*M-1); i++)
  		{
  			for (j=0; j<(2*M-1); j++)
  			{
  			
  			if(i == row || j == row || i == column || j == column || i==j )
  				(*distance_matrix)[i][j]=1000.0;
  			}
		}
		// Aggregate the two layers with indexes (row,column)  and put such layer as the M+1 element in the multiplex structure
		aggregate(&multiplex,M+1,row,column,MAX_NODE,verbose);


		// l1 and l2 are the two layers that will be merged in this round
		multiplex[Mmod].l1=row;
		multiplex[Mmod].l2=column;

		// Vectorprint keeps track of the layers that are going to be merged (used for the printing procedure)
		(*vectorprint)[row]=1;
		(*vectorprint)[column]=1;
		if (verbose || partition_verbose)
		{
			for (i = 0; i < M; i++)
			{
				if ( (*vectorprint)[i] == 0)
				printf("%d,  ", multiplex[i].l1);
			}
			for (i=Mmod; i >= M; i--)
			{
				printf("( %d %d ), ", multiplex[i].l1,multiplex[i].l2);
			}
			if (partition_verbose)
				printf("-->");

		}

		layers_order=(unsigned int *)malloc((M-1)*sizeof(unsigned int ));

		//Layers_order contains the indexes of the layers composing the multiplex after the merging process
		j=0;
		for (i = 0; i <M; i++)
		{
			if ( (*vectorprint)[i] == 0)
				layers_order[j++]=multiplex[i].l1;
		}
		layers_order[j++]=M;


		//***********Evaluate the Kolmogorov Complexity of the multiplex with M-1 layers************
		new_multiplex=(struct layer *)malloc((M-1)*sizeof(struct layer ));
		for (i = 0; i < M-1; i++)
			new_multiplex[i]=multiplex[layers_order[i]];

		evaluate_KC(new_multiplex,M-1,MAX_NODE,primeslist,nodes_unique,length_unique,ITER,parallel);
		free(layers_order);
		free(new_multiplex);
		//********************************************************************************************
	reduce_multiplex(root,multiplex,M,M+1,MAX_NODE,primeslist,nodes_unique,length_unique,1,distance_matrix,vectorprint,ITER,parallel,verbose,partition_verbose);	
	}
	else
	{
		duplex=(struct layer *)malloc(2*sizeof(struct layer));
		flag_intern=0;
		max_KC_step=0.0;
		duplex[0].nlinks=multiplex[Mmod-1].nlinks;
		duplex[0].links_table=multiplex[Mmod-1].links_table;
		duplex[0].max_node=multiplex[Mmod-1].max_node;
		for (i = 0; i< Mmod-1; i++)
		{
			if((*distance_matrix)[Mmod-1][i]!=1000.0)
			{
	
				duplex[1].nlinks=multiplex[i].nlinks;
				duplex[1].links_table=multiplex[i].links_table;
				duplex[1].max_node=multiplex[i].max_node;
				sumKC=0;
				sumKCaggregate=0;
				nedges_intersection=0;
				cumulative_complexity =0.0;
				omp_set_num_threads(parallel);
				//Allocation for the parallel computing
				temporarystruct = (struct layer **) malloc (parallel*sizeof(struct layer *));
				nthread_root = (void **) malloc (parallel*sizeof(void *));
				Nlinks_aggregate= (unsigned long long int *)malloc(parallel*sizeof(unsigned long long int ));
				Nlinks_multiplex= (unsigned long long int *)malloc(parallel*sizeof(unsigned long long int ));
				dimension_string= (unsigned long long int *)malloc(parallel*sizeof(unsigned long long int ));
				string_counter= (unsigned long long int *)malloc(parallel*sizeof(unsigned long long int ));
				string_counter_aggregate= (unsigned long long int *)malloc(parallel*sizeof(unsigned long long int ));
				fp_string= (char **)malloc(parallel*sizeof(char *));
				fp_string_aggregate= (char **)malloc(parallel*sizeof(char *));
				maximum_overlap=(unsigned int *)malloc(parallel*sizeof(unsigned int ));
				maximum_prime_product=(mpz_t *)malloc(parallel*sizeof(mpz_t));
				for (k=0;k<parallel;k++)
					mpz_init(maximum_prime_product[k]);
				#pragma omp parallel for private (nthread)
				for (k=0;k<ITER;k++)
				{
					nthread= omp_get_thread_num();
					nthread_root[nthread]=NULL;
					length=0;
					prime_layer=primeslist[1];
					temporarystruct[nthread]=multiplex_into_binarytree(&(nthread_root[nthread]),duplex,2,MAX_NODE,primeslist,nodes_unique,length_unique,&nedges_intersection,nthread);
					if (maximum_overlap[nthread] == 0)
						maximum_overlap[nthread]=1;

					//Printing the duplex and the corresponding aggregate
					dimension_KC=print_file_bothfast_onstring(nthread_root[nthread],&dimension_KC_aggregate,nthread);
					sumKC+=dimension_KC;
					sumKCaggregate+=dimension_KC_aggregate;
					cumulative_complexity+=(1.0*dimension_KC)/(1.0*dimension_KC_aggregate);
					tdestroy(nthread_root[nthread]);
					free_structgmp_multiplex(temporarystruct[nthread],multiplex,2);
				}
				sumKC=(1.0*sumKC)/(1.0*ITER);
				sumKCaggregate=(1.0*sumKCaggregate)/(1.0*ITER);
				fraction = (1.0*cumulative_complexity)/(1.0*ITER);


				(*distance_matrix)[Mmod-1][i]=fraction;
				(*distance_matrix)[i][Mmod-1]=fraction;
				omega[Mmod-1][i]=((1.0*nedges_intersection)/(1.0*ITER))/(1.0*length);
				omega[i][Mmod-1]=omega[Mmod-1][i];
				free(temporarystruct);
				free(nthread_root);
				free(Nlinks_aggregate);
				free(Nlinks_multiplex);
				free(dimension_string);
				free(string_counter);
				free(string_counter_aggregate);
				free(fp_string);
				free(fp_string_aggregate);
				free(maximum_overlap);
				for (k=0;k<parallel;k++)
					mpz_clear(maximum_prime_product[k]);
				free(maximum_prime_product);
				}

		}
		(*distance_matrix)[Mmod-1][Mmod-1]=1000.0;
	
		// Find the duplex (couple of layers) \bar{D} with the highest value of complexity, this is the reference point
		for (i = 0; i < Mmod-1; i++)
		{
			for (j = i+1; j < Mmod; j++)
			{
				if( (*distance_matrix)[i][j] >= max_KC_step && (*distance_matrix)[i][j]!= 1000.0)
					{
						max_KC_step = (*distance_matrix)[i][j];
						row=i;
						column=j;
					}
			}
		}
		//*********************************************************************************************************************************
		// Finding the couple of layers D with the smallest value of complexity  (among the set of feasible candidates that have a value of 
		// overlap greater or equal than the overlap of \bar(D))
		min=1000.0;
		for(i=0; i < Mmod-1;i++)
			for (j = i+1; j <  Mmod; j++)
			{
				if( (*distance_matrix)[i][j] <= max_KC_step && omega[i][j] >= omega[row][column] && (*distance_matrix)[i][j] < min)
				{
					min=(*distance_matrix)[i][j];
					//flag_intern=1;
					new_row=i;
					new_column=j;
				}


			}
		row=new_row;
		column=new_column;


		if (verbose)
		fprintf(stderr,"\n\nThe index with the minimum value %lf is %u %u\n",min,row, column);
		for (i=0; i<(2*M-1); i++)
  		{
  			for (j=0; j<(2*M-1); j++)
  			{
  			
  			if(i == row || j == row || i == column || j == column )
  				(*distance_matrix)[i][j]=1000.0;
  			}
		}	
		multiplex[Mmod].l1=row;
		multiplex[Mmod].l2=column;
		
		if (row < M)
		(*vectorprint)[row]=1;
		if (column < M)
		(*vectorprint)[column]=1;
		int *copy=(int *)calloc((2*M-1),sizeof(int));
		for ( i = 0; i < 2*M-1; i++)
		{
			copy[i]=(*vectorprint)[i];
		}

		//If the number of layers are greater than 2 continue, otherwise stop. You have done all the steps
		if ( M > 2)
		{
			aggregate(&multiplex,Mmod+1,row,column,MAX_NODE,verbose);	
			//Layers order contains the indexes of the layers composing the multiplex after the merging process

			step=Mmod-M+1;
			layers_order=(unsigned int *)malloc((M-step)*sizeof(unsigned int ));
			j=0;
			for (i = 0; i < M; i++)
				{
					if ( copy[i] == 0)
					{
					if (verbose || partition_verbose)
					printf("%d,  ", multiplex[i].l1);
					layers_order[j++]=multiplex[i].l1;
					copy[i]=1;
					}
				}
			for (i=Mmod; i >= M; i--)
			{
				if(copy[i]==0)
				{
					layers_order[j++]=i;
					
					if (verbose || partition_verbose)
					printf("( ");
					recursive_print(multiplex,i,i,M,Mmod,&copy,verbose,partition_verbose);
					if (verbose|| partition_verbose)
					printf("), ");
				}
			}
			free(copy);
			if (partition_verbose)
				printf("-->");
			//Layers order contains the indexes of the layers composing the multiplex after the merging process
			//***********Evaluate the Kolmogorov Complexity of the multiplex with M-1 layers************
			new_multiplex=(struct layer *)malloc((M - step)*sizeof(struct layer ));
			for (i = 0; i < M-step; i++)
				new_multiplex[i]=multiplex[layers_order[i]];
			evaluate_KC(new_multiplex,M-step,MAX_NODE,primeslist,nodes_unique,length_unique,ITER,parallel);
			free(layers_order);
			free(new_multiplex);
			free(duplex);
		}

		if (Mmod+1 < 2*M-1)
		{
			reduce_multiplex(root,multiplex,M,Mmod+1,MAX_NODE,primeslist,nodes_unique,length_unique,1,distance_matrix,vectorprint,ITER,parallel,verbose,partition_verbose);
		}
		else
		{
			if(verbose)
			fprintf(stderr,"\nFinish!\n");
		}
	}
}






void evaluate_KC(struct layer *multiplex, unsigned int M,unsigned long long int MAX_NODE, unsigned long long int *primeslist, unsigned long long int *nodes_unique, unsigned long long int length_unique, unsigned int ITER,unsigned int parallel)
{
	
	unsigned long long int dimension_KC,dimension_KC_aggregate,k,nedges_intersection;
	unsigned int nthread;
	struct layer **temporarystruct;
	void **nthread_root;
	double sumKC,sumKCaggregate,complexity,numerator,denominator,cumulative_complexity;

	//***********Evaluate the Kolmogorov Complexity of the whole multiplex************
	sumKC=0;
	sumKCaggregate=0;
	nedges_intersection=0;
	cumulative_complexity=0.0;
	omp_set_num_threads(parallel);
	//Allocation for the parallel computing
	temporarystruct = (struct layer **) malloc (parallel*sizeof(struct layer *));
	nthread_root = (void **) malloc (parallel*sizeof(void *));
	Nlinks_aggregate= (unsigned long long int *)malloc(parallel*sizeof(unsigned long long int ));
	Nlinks_multiplex= (unsigned long long int *)malloc(parallel*sizeof(unsigned long long int ));
	dimension_string= (unsigned long long int *)malloc(parallel*sizeof(unsigned long long int ));
	string_counter= (unsigned long long int *)malloc(parallel*sizeof(unsigned long long int ));
	string_counter_aggregate= (unsigned long long int *)malloc(parallel*sizeof(unsigned long long int ));
	fp_string= (char **)malloc(parallel*sizeof(char *));
	fp_string_aggregate= (char **)malloc(parallel*sizeof(char *));
	maximum_overlap=(unsigned int *)malloc(parallel*sizeof(unsigned int ));
	maximum_prime_product=(mpz_t *)malloc(parallel*sizeof(mpz_t));	
	for (k=0;k<parallel;k++)
		mpz_init(maximum_prime_product[k]);
	//Parallel computation
	#pragma omp parallel for private (nthread)
	for (k=0;k<ITER;k++)
		{
			nthread= omp_get_thread_num();
			nthread_root[nthread]=NULL;
			length=0;
			Nlinks_aggregate[nthread]=0;
			Nlinks_multiplex[nthread]=0;
			temporarystruct[nthread]=multiplex_into_binarytree(&(nthread_root[nthread]),multiplex,M,MAX_NODE,primeslist,nodes_unique,length_unique,&nedges_intersection,nthread);
			//Printing the whole multiplex  & the corresponding aggregate
			prime_layer=primeslist[M-1];
			dimension_KC=print_file_bothfast_onstring(nthread_root[nthread],&dimension_KC_aggregate,nthread);
			sumKC+=dimension_KC;
			sumKCaggregate+=dimension_KC_aggregate;
			cumulative_complexity+=(1.0*dimension_KC)/(1.0*dimension_KC_aggregate);
			tdestroy(nthread_root[nthread]);
			free_structgmp_multiplex(temporarystruct[nthread],multiplex,M);
		}
	numerator=(1.0*sumKC)/(1.0*ITER);
	denominator= (1.0*sumKCaggregate)/(1.0*ITER);
	
	complexity=cumulative_complexity/(1.0*ITER);
	printf("%d\t%lf\t%llu\t%lf\t%llu\t%lf\t%lf\n", M,(1.0*sumKC)/(1.0*ITER),Nlinks_multiplex[0], (1.0*sumKCaggregate)/(1.0*ITER),Nlinks_aggregate[0],complexity,complexity/log(Nlinks_multiplex[0]));
	complexity_multiplex_reducibility[M-1]=complexity/log(Nlinks_multiplex[0]);
	free(temporarystruct);
	free(nthread_root);
	free(Nlinks_aggregate);
	free(Nlinks_multiplex);
	free(dimension_string);
	free(string_counter);
	free(string_counter_aggregate);
	free(fp_string);
	free(fp_string_aggregate);
	free(maximum_overlap);
	for (k=0;k<parallel;k++)
		mpz_clear(maximum_prime_product[k]);
	free(maximum_prime_product);
}


struct layer *multiplex_into_binarytree(void **root,struct layer *multiplex, unsigned int M, unsigned long long int MAX_NODE,
 unsigned long long int *primeslist, unsigned long long int *nodes_unique, unsigned long long int length_unique,
 unsigned long long int  *length_nedges_intersection, unsigned int nthread)
{
 
	unsigned long long int  i,j,k, nodei,nodej,*shuffle;
	struct edge *update;
	unsigned long long int *distr1;
	struct layer *temporarystruct;
	distr1=particular_primes_sequence(multiplex, M);
	shuffle=shuffle_nodes(length_unique,nodes_unique);
	temporarystruct = update_nodelabel(multiplex,shuffle,M,length_unique,nodes_unique);
	maximum_overlap[nthread]=0;
	mpz_set_ui(maximum_prime_product[nthread],primeslist[M-1]);
	
	for (k=0;k<M;k++)
	{
		for (i = 0; i < temporarystruct[distr1[k]].nlinks; i++)
		{
			nodei=temporarystruct[distr1[k]].links_table[i].node_i;
			nodej=temporarystruct[distr1[k]].links_table[i].node_j;
			if(nodei > nodej)
			{
				temporarystruct[distr1[k]].links_table[i].node_i = temporarystruct[distr1[k]].links_table[i].node_j;
				temporarystruct[distr1[k]].links_table[i].node_j = nodei;
				nodei=temporarystruct[distr1[k]].links_table[i].node_i;
				nodej=temporarystruct[distr1[k]].links_table[i].node_j;
			}
			Nlinks_multiplex[nthread]++;
			temporarystruct[distr1[k]].links_table[i].key= MAX_NODE *  (nodei) + (nodej);
			update=tfind (& temporarystruct[distr1[k]].links_table[i], root, compare_keys);
			// If the element does not belong to the tree, then update is egual to NULL and the element will be
			// inserted with the function tsearch
			if( update== NULL)
			{
				Nlinks_aggregate[nthread]++;
				mpz_set_ui(temporarystruct[distr1[k]].links_table[i].weight_link, primeslist[k]);
				temporarystruct[distr1[k]].links_table[i].overlap_aggregate++;
				(void) tsearch(& temporarystruct[distr1[k]].links_table[i], root, compare_keys);
			}
			else
			{
				//If the element already exists in the tree then//
  				//update the prime weight of the link considered//
  				(*length_nedges_intersection)++;
  				update=*(struct edge **)update;	
  				(*update).overlap_aggregate++;
				mpz_mul_ui((*update).weight_link,(*update).weight_link,primeslist[k]);
				if (mpz_cmp_flag_to_01(mpz_cmp((*update).weight_link,maximum_prime_product[nthread])))
        {
          //gmp_printf("changing %Zd in %Zd\n",maximum_prime_product[nthread],(*update).weight_link);
          mpz_set(maximum_prime_product[nthread],(*update).weight_link);  
        }
			}
		}
	}
	free(distr1);
	free(shuffle);
	return(temporarystruct);
}



unsigned long long int *particular_primes_sequence(struct layer *multiplex, unsigned int M)
{
	unsigned long long int *distr1,*distr2,*temp,i,j;

	distr1=(unsigned long long int*)malloc(M*sizeof(unsigned long long int));
	distr2=(unsigned long long int*)malloc(M*sizeof(unsigned long long int));
	temp=(unsigned long long  int*)malloc(M*sizeof(unsigned long long int));

	for(i=0;i<M;i++)
	{
		distr1[i]=multiplex[i].nlinks;
		distr2[i]=0;
	}
	//increasing sequence
	qsort(distr1, M, sizeof(unsigned long long int), cmpfunc);
	for(i=0;i<M;i++)
		for(j=0;j<M;j++)
		{
			if(distr1[i]==multiplex[j].nlinks && distr2[j]==0)
			{
				temp[i]=j;
				distr2[j]+=1;
        break;  
			}
		}

	free(distr1);
	free(distr2);
	return(temp);
}


unsigned long long int *shuffle_nodes(unsigned long long int length_unique, unsigned long long int *nodes_unique)
{
	unsigned long long int *distr,i,j,temp;
	distr=(unsigned long long int*)malloc(length_unique*sizeof(unsigned long long int));
	for(i=0;i<length_unique;i++)
		distr[i]=i;

	//shuffle the label of the nodes
	for (i=length_unique-1;i>0;i--)
	{
		j=rand()%(i+1);
		temp = distr[i];
		distr[i]=distr[j];
		distr[j]=temp;
	}
		//for(i=0;i<length_unique;i++)
		//printf("%Ld  ", distr[i]);
	return(distr);
}



struct layer *update_nodelabel(struct layer *multiplex, unsigned long long int *shufflenodes, unsigned int M,unsigned long long int length_unique, unsigned long long int *nodes_unique)
{
	unsigned int i;
	unsigned long long int nodei,nodej,j;
	struct layer *new_duplex= (struct layer *) malloc(M * sizeof(struct layer));
	for (i=0; i<M;i++)
	{
		new_duplex[i].nlinks=multiplex[i].nlinks;
		new_duplex[i].max_node=multiplex[i].max_node;
		new_duplex[i].links_table = (struct edge *) malloc (multiplex[i].nlinks * sizeof(struct edge));
		for(j=0;j<multiplex[i].nlinks;j++)
		{
			nodei=multiplex[i].links_table[j].node_i;
			new_duplex[i].links_table[j].node_i=shufflenodes[nodei];

			nodej=multiplex[i].links_table[j].node_j;
			new_duplex[i].links_table[j].node_j=shufflenodes[nodej];
			mpz_init(new_duplex[i].links_table[j].weight_link);
			new_duplex[i].links_table[j].overlap_aggregate=0;
		}
	}
	return(new_duplex);
}



void free_structgmp_multiplex(struct layer *temporarystruct, struct layer *multiplex, unsigned int M)
{
	unsigned int i;
	unsigned long long int j;
	for (i = 0; i <M; i++)
	{
		
		for (j = 0; j < temporarystruct[i].nlinks; j++)
		{
			mpz_clear(temporarystruct[i].links_table[j].weight_link);
		}
		free(temporarystruct[i].links_table);	
	}

	free(temporarystruct);	
}


void tdestroy(void *root)
{
	 struct edge * elementptr;
    while (root != NULL) {
        elementptr = *(struct edge **)root;
        tdelete((void *)elementptr, &(root), delete_root); 
    }
}

int delete_root(const void *node1, const void *node2)
{
    return 0;
}


int compare_keys(const void *e1p, const void *e2p)
{
	const struct edge *e1, *e2;
	long long int last, first;

	e1 = (const struct edge *) e1p;
	e2 = (const struct edge *) e2p;
	/* check key values */
	if (e1->key < e2->key)
		return -1;
	else if (e1->key == e2->key)
		return 0;
	else
		return 1;
}


int mpz_cmp_flag_to_01(int mpz_cmp_flag){
  if(mpz_cmp_flag <= 0)
    return 0;
  else
    return 1;
}


unsigned long long int print_file_bothfast_onstring(void **root,unsigned long long int *dimension_aggregate,  unsigned int thread_number)
{
	unsigned long long int dimension=0;
	dimension_string[thread_number]=100000;
	fp_string[thread_number]= (char *) malloc (dimension_string[thread_number] * sizeof(char));
	fp_string_aggregate[thread_number]= (char *)malloc (dimension_string[thread_number] * sizeof(char));
	string_counter[thread_number]=0;
	string_counter_aggregate[thread_number]=0;

	twalk(root,print_elements_onstring_weighted_aggregate);
	dimension = dynamic_compress(fp_string[thread_number],string_counter[thread_number]+1);
	*dimension_aggregate = dynamic_compress(fp_string_aggregate[thread_number],string_counter_aggregate[thread_number]+1);
	free(fp_string[thread_number]);
	free(fp_string_aggregate[thread_number]);
	return(dimension);
}




void print_elements_onstring_weighted_aggregate(const void *nodep, const VISIT which, const int depth)
{
	
	struct edge *e = *((struct edge **) nodep);

	switch (which) {
	case leaf:
	case postorder:
		length++;
		unsigned int thread_number= omp_get_thread_num();
		mpz_t temporary_mpz;
		if( (string_counter[thread_number]  > (dimension_string[thread_number] - 5000) && string_counter[thread_number] < dimension_string[thread_number]) ||  
			(string_counter_aggregate[thread_number]  > (dimension_string[thread_number] - 5000) && string_counter_aggregate[thread_number] < dimension_string[thread_number] ))
		{
			dimension_string[thread_number]+=100000;	
			fp_string[thread_number]= (char *)realloc (fp_string[thread_number],dimension_string[thread_number] * sizeof(char));
			fp_string_aggregate[thread_number]= (char *)realloc (fp_string_aggregate[thread_number],dimension_string[thread_number] * sizeof(char));
		}
		//Printing on the string associated to the multiplex, the edge list (i,j,Omega_ij) in lexicographic order
		string_counter[thread_number] += gmp_sprintf(fp_string[thread_number]+string_counter[thread_number],"%llu %llu %Zd\n", e->node_i, e->node_j,e->weight_link);

		mpz_init(temporary_mpz);
    	mpz_set(temporary_mpz,maximum_prime_product[thread_number]);
    	//Printing on the string associated to the aggregate, the edge list (i,j,max(Omega_ij)) in lexicographic order
		string_counter_aggregate[thread_number] += gmp_sprintf(fp_string_aggregate[thread_number]+string_counter_aggregate[thread_number],"%llu %llu %Zd\n", e->node_i, e->node_j,temporary_mpz);
		mpz_clear(temporary_mpz);
		
		break;
	default:
		break;
	}
}

unsigned long long int  dynamic_compress(char* str, unsigned long long int size)
{
  char *temporary_dyn;
  unsigned long long int dimension;
  temporary_dyn= (char *)malloc(size*sizeof(char));
	// deflate
	// zlib struct
	z_stream defstream;
	defstream.zalloc = Z_NULL;
	defstream.zfree = Z_NULL;
	defstream.opaque = Z_NULL;
	defstream.avail_in = (uInt)strlen(str)+1; // size of input, string + terminator
	defstream.next_in = (Bytef *)str; // input char array
	//printf("Size of the output:%d\n", size);
	defstream.avail_out = (uInt)size; // size of output
	defstream.next_out = (Bytef *)temporary_dyn; // output char array

	deflateInit(&defstream, Z_DEFAULT_COMPRESSION);
	deflate(&defstream, Z_FINISH);
	deflateEnd(&defstream);

	// This is one way of getting the size of the output
	//printf("Deflated size is: %lu\n", (char *) defstream.next_out - temporary_dyn);
	dimension= (char *) defstream.next_out - temporary_dyn;
	free(temporary_dyn);
	return(dimension);
}


void aggregate(struct layer **multiplex,unsigned int Nlayers, unsigned int row, unsigned int column, unsigned long long int MAX_NODE,unsigned int verbose)
{
	void *root=NULL;
	unsigned int *layer_to_collapse,k;
	unsigned long long int num_links=0,max_nodelayer=0,size=1000,a,b,i;
	struct edge *update,**new_link;
	layer_to_collapse = (unsigned int *)malloc(2*sizeof(unsigned int));
	layer_to_collapse[0]=row;
	layer_to_collapse[1]=column;
	new_link = &(*multiplex)[Nlayers-1].links_table;
	*new_link = (struct edge *) malloc(size * sizeof(struct edge ));
	if (verbose)
	{
	fprintf(stderr,"%u---Link: %llu\n",layer_to_collapse[0],(*multiplex)[layer_to_collapse[0]].nlinks);
	fprintf(stderr,"%u---Link: %llu\n",layer_to_collapse[1],(*multiplex)[layer_to_collapse[1]].nlinks);
	}

	for (k=0;k<2;k++)
	{
		for (i = 0; i < (*multiplex)[layer_to_collapse[k]].nlinks; i++)
		{
			a=(*multiplex)[layer_to_collapse[k]].links_table[i].node_i;
			b=(*multiplex)[layer_to_collapse[k]].links_table[i].node_j;
			(*multiplex)[layer_to_collapse[k]].links_table[i].key= a*MAX_NODE+b;
			update=tfind (& (*multiplex)[layer_to_collapse[k]].links_table[i], &root, compare_keys);
		
			// If the element does not belong to the tree, then update is egual to NULL and the element will be
			// inserted with the function tsearch
			if( update == NULL)
			{
				
				//Find the max node between the two layers to aggregate
				if(max_nodelayer < b )
					max_nodelayer = b;
				//Unweighted union of the two layers
				(void) tsearch(& (*multiplex)[layer_to_collapse[k]].links_table[i], &root, compare_keys);
				(*new_link)[num_links].node_i= a;
				(*new_link)[num_links].node_j= b;
				(*new_link)[num_links].key= MAX_NODE *  (a) + (b);
				(*new_link)[num_links].overlap_aggregate=0;
				mpz_init((*new_link)[num_links].weight_link);
				//mpz_set_ui((*new_link)[num_links].weight_link, size);
				//printf("unione: %Ld  %Ld\n",(*new_link)[num_links].node_i,(*new_link)[num_links].node_i);
				num_links++;
				if(num_links == size)
				{
					size+=1000;
					(*new_link)=(struct edge *)realloc((*new_link),size*sizeof(struct edge));
				}
			}
		}

	}
	(*multiplex)[Nlayers-1].links_table=*new_link;
	(*multiplex)[Nlayers-1].nlinks=num_links;
	if (verbose)
	fprintf(stderr,"Number of links aggregate:%llu\n\n", num_links);
	(*multiplex)[Nlayers-1].max_node=max_nodelayer;
	free(layer_to_collapse);

	tdestroy(root);

}	



void recursive_print (struct layer *multiplex, unsigned int i,unsigned int iprev, unsigned int M, unsigned int Mmod, int **copy,unsigned int verbose,unsigned int partition_verbose)
{
	unsigned int j,k,flag=0;
	if(i!= iprev)
	{
		(*copy)[iprev]=i;
	}

	if ((*copy)[i]==0)
	{

		if(multiplex[i].l1 < M)
		{
			if (verbose|| partition_verbose)
			printf("%d ",multiplex[i].l1 );
		}
		else
		{
			recursive_print(multiplex, multiplex[i].l1,i, M,Mmod, copy,verbose,partition_verbose);
		}

		if(multiplex[i].l2 < M)
		{
			if (verbose|| partition_verbose)
			printf("%d ",multiplex[i].l2 );
		}
		else
		{
			recursive_print(multiplex, multiplex[i].l2,i, M,Mmod, copy,verbose,partition_verbose);
		}
	}

	if(multiplex[i].l1 < M && multiplex[i].l2 <M )
	{
		(*copy)[i]=1;
	}

}



void free_structure(unsigned long long int *list_prime,unsigned long long int *nodes_unique,struct layer *multiplex, unsigned int M)
{
	unsigned long long int i,k;
	for (k=0;k<M;k++)
	{
		for (i = 0; i < multiplex[k].nlinks; i++)
			mpz_clear(multiplex[k].links_table[i].weight_link); 
		free(multiplex[k].links_table);
	}
	free(list_prime);
	free(nodes_unique);
	free(multiplex);
}


int *compute_sequence_layers(struct layer *multiplex,unsigned int index_maximum_reducibility,unsigned int M)
{
	int *vector_layers_ID,*indexes_not_to_print,*indexes_to_print,*indexes_to_evaluate,i,j,r,s,t,flag=0;
	unsigned  int maximum_layer_index=2*(M-1)-index_maximum_reducibility;
	unsigned int number_layer_optimal_multiplex= index_maximum_reducibility+1;
	vector_layers_ID = (int *)malloc ((number_layer_optimal_multiplex)*sizeof(int));
	indexes_not_to_print = (int *) malloc (M*sizeof(int));
	indexes_to_print = (int *)malloc (M*sizeof(int));
	indexes_to_evaluate = (int *)malloc (M*sizeof(int));

	for (i=0;i<M; i++)
	{
		indexes_not_to_print[i]=-10;
		indexes_to_print[i]=-10;
		indexes_to_evaluate[i]=-10;
	}
	r=0;
	s=0;
	t=0;
	//printf("maximum_layer_index:%d\n",maximum_layer_index);
	for (i=maximum_layer_index; i>=M ;i --)
	{
		flag =0;
		// if i belongs to indexes_to_evaluate -> flag =1
		for (j =0;j<s;j++)
			if (i == indexes_to_evaluate[j])
				flag =1;

		if (flag == 0)
		{
			indexes_to_print[t]=i;
			t++;
		}

		if (multiplex[i].l1 <  M)
		{
			indexes_not_to_print[r]=multiplex[i].l1;
			r++;
		}
		else
		{
			indexes_to_evaluate[s]=multiplex[i].l1;
			s++;
		}
		if (multiplex[i].l2 <  M)
		{
			indexes_not_to_print[r]=multiplex[i].l2;
			r++;
		}
		else
		{
			indexes_to_evaluate[s]=multiplex[i].l2;
			s++;
		}

	}
	

	qsort(indexes_not_to_print,r,sizeof(int),cmpfunc);
	qsort(indexes_to_print,t,sizeof(int),cmpfunc);
	j=0;
	s=0;
	for (i =0;i<M;i++)
	{
		if (indexes_not_to_print[s]!=i)
			vector_layers_ID[j++]=i;
		else
			s++;
	}
	for (i=0;i<t;i++)
		vector_layers_ID[j+i]=indexes_to_print[i];

	// for (i=0;i<number_layer_optimal_multiplex;i++)
	// 	printf("%d ",vector_layers_ID[i]);
	// printf("\n");
	return(vector_layers_ID);
}


void dump_multiplex_on_file(struct layer *multiplex,int *sequence_multiplex_toprint,
														unsigned int numberoflayers_maximum_reducibility, unsigned long long int *nodes_unique, unsigned long long int length_unique)
{

	unsigned int i,j,k;
	unsigned long long int *p_a=&(nodes_unique[length_unique]);
	unsigned long long int *p_b;
 	unsigned long long int differenceInBytes,temp1,temp2,index1,index2;
	char buff[256];
	FILE *fp;
	for(k=0;k<numberoflayers_maximum_reducibility;k++)
	{
		sprintf(buff, "layer%d_optimal.txt",k);
		fp = fopen(buff,"w+");
		for (i=0;i<multiplex[sequence_multiplex_toprint[k]].nlinks;i++)
		{
			index1 = multiplex[sequence_multiplex_toprint[k]].links_table[i].node_i;
			index2 = multiplex[sequence_multiplex_toprint[k]].links_table[i].node_j;
			temp1=nodes_unique[index1];
			temp2=nodes_unique[index2];
			fprintf(fp,"%llu %llu\n",temp1,temp2);	
		}
		fclose(fp);
	}

}





void dump_all_multiplex_on_file(struct layer *multiplex,int *sequence_multiplex_toprint,
														unsigned int numberoflayers_maximum_reducibility, unsigned long long int *nodes_unique, unsigned long long int length_unique)
{

	unsigned int i,j,k;
	unsigned long long int *p_a=&(nodes_unique[length_unique]);
	unsigned long long int *p_b;
 	unsigned long long int differenceInBytes,temp1,temp2,index1,index2;
	char buff[256];
	FILE *fp;
	for(k=0;k<numberoflayers_maximum_reducibility;k++)
	{
		sprintf(buff, "multiplexM=%d_layer%d.txt",numberoflayers_maximum_reducibility,k);
		fp = fopen(buff,"w+");
		for (i=0;i<multiplex[sequence_multiplex_toprint[k]].nlinks;i++)
		{
			index1 = multiplex[sequence_multiplex_toprint[k]].links_table[i].node_i;
			index2 = multiplex[sequence_multiplex_toprint[k]].links_table[i].node_j;
			temp1=nodes_unique[index1];
			temp2=nodes_unique[index2];
			fprintf(fp,"%llu %llu\n",temp1,temp2);	
		}
		fclose(fp);
	}

}