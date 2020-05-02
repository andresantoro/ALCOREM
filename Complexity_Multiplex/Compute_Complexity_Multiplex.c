#include <search.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <zlib.h>
#include <omp.h>
#include <gmp.h>


// These variables are required for the computation of the Complexity 
unsigned long long int length; 
unsigned long long int prime_layer;
unsigned int *maximum_overlap;
mpz_t *maximum_prime_product;
unsigned long long int *Nlinks_aggregate;
unsigned long long int *Nlinks_multiplex;
char **fp_string;
char **fp_string_aggregate;
FILE *fp_save_prime_weight;
unsigned long long int *dimension_string;
unsigned long long int *string_counter;
unsigned long long int *string_counter_aggregate;

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
	"**              Computation of the Kolmogorov complexity                    **\n"
	"**                         of a multiplex system                            **\n"
	"**                                                                          **\n"
	"**                                                                          **\n"
	"**      <unique_nodes> list containing the ordered list of node IDs         **\n"    
	"**                              in the multiplex                            **\n"
	"**                                                                          **\n"
	"**           <filein> should be a list of the layers composing the          **\n"
	"**                            multiplex network                             **\n"
	"**                                                                          **\n"
	"**                     ----   Optional Variables  ----                      **\n"
	"**                                                                          **\n"
	"**           <ITER> represents the number of times the Kolmogorov           **\n"
	"**                    complexity is averaged (default = 1000)               **\n"
    "**                                                                          **\n"
	"**    <#core> represents the number of cores used in the computation of     **\n"
	"**        the complexity   (default: sequential computation)                **\n"
	"**                                                                          **\n"
	"**    <-s 'filename'> save the prime-weight matrix in the file 'filename'   **\n"
	"**                         as an ordered weighted edge list                 **\n"
	"**                                                                          **\n"
	"**               <-v>  be more verbose -- show misc extra info              **\n"
    "**      OUTPUT: by default the algorithm returns the following info:        **\n"
    "**  KC_multiplex, NLinks_multiplex,KC_aggregate,Nlinks_aggregate,Complexity **\n"
    "**                                                                          **\n"	
	"******************************************************************************\n");
  printf("Usage: %s <unique_nodes> <filein>  [ITER = 100] [-p #core] [-s <filename>]\n\n" , argv[0]);
}

unsigned long long int *load_unique_nodes(char *, unsigned long long int *, unsigned long long int *,unsigned int);

struct layer *load_data(char *, unsigned int *, unsigned long long int *, unsigned long long int , unsigned int);

unsigned long long int load_edgelist(char *, struct edge **, unsigned long long int *, unsigned long long int );

int cmpfunc (const void *, const void *);

void list_prime(unsigned long long int **,unsigned int );

void allocate_memory(struct layer *,unsigned int);

void evaluate_KC(struct layer *, unsigned int,unsigned long long int, unsigned long long int *,
									unsigned long long int *, unsigned long long int, unsigned int,unsigned int ,unsigned int, char *,unsigned int );
struct layer *multiplex_into_binarytree(void **,struct layer *, unsigned int, unsigned long long int ,
																				unsigned long long int *, unsigned long long int *, unsigned long long int , 
																				unsigned long long int  *, unsigned int);

unsigned long long int *particular_primes_sequence(struct layer *, unsigned int );

unsigned long long int *shuffle_nodes(unsigned long long int , unsigned long long int *);

struct layer *update_nodelabel(struct layer *, unsigned long long int *, unsigned int,unsigned long long int , unsigned long long int *);

int mpz_cmp_flag_to_01(int);

void tdestroy(void *);

int delete_root(const void *, const void *);

int compare_keys(const void *, const void *);

void free_structgmp_multiplex(struct layer *, struct layer *, unsigned int );

unsigned long long int print_file_bothfast_onstring(void **,unsigned long long int *,  unsigned int );

void print_elements_onstring_weighted_aggregate(const void *, const VISIT, const int );

unsigned long long int  dynamic_compress(char* , unsigned long long int );

void free_structure(unsigned long long int *,unsigned long long int *,struct layer *, unsigned int );

void chopN_string(char *, size_t );



struct layer *multiplex_into_binarytree_saving(void **,struct layer *, unsigned int , unsigned long long int,
 			unsigned long long int *, unsigned long long int *, unsigned long long int , unsigned long long int  *, unsigned int);

struct layer *update_nodelabel_saving(struct layer *,unsigned int,unsigned long long int , unsigned long long int *);

void print_elements_onfile_saving(const void *, const VISIT , const int );




int main(int argc, char *argv[])
{
	unsigned int M,i,j,k,ITER,parallel=1,save_prime_weight=0,verbose=0;
	unsigned long long int *nodes_unique,length_unique, MAX_NODE,*primeslist;
	struct layer *multiplex;
  	unsigned int randval,esc;
  	char *str_par,*prime_weight_filename;
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
	if (argc <= 2)
	{
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

	
	for (i = 1; i < argc; i++) 
	{
	  //Checking for the parallel option
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
 		//Checking for the verbose option
 		if (argv[i][0] == '-' && argv[i][1] == 'v') 
 		{
 			fprintf(stderr,"Verbose activated\n");
 			verbose = 1;
 		}
 		//Checking for saving the prime-weight edge list option
 	 	if (argv[i][0] == '-' && argv[i][1] == 's') 
 		{
 			fprintf(stderr,"Save prime-weight edge list activated\n");
 			save_prime_weight = 1;
 			prime_weight_filename = argv[i+1];
 		}
 	}

 	if(verbose)
 	fprintf(stderr,"Finding the number of nodes...\n");
	nodes_unique=load_unique_nodes(argv[1],&length_unique, &MAX_NODE,verbose);
	multiplex=load_data(argv[2],&M,nodes_unique,length_unique,verbose);
	if(verbose)
	fprintf(stderr,"Computing the prime list...\n");
	primeslist=(unsigned long long int *)malloc(M*sizeof(unsigned long long  int));
	list_prime(&primeslist,M);
	allocate_memory(multiplex,M);
	// //****************************************************
	//Check if the data has been loaded in the correct way
	// for (i=0; i<M;i++)
	// {
	// 	printf("Layer %u\n\n",i+1);
	// 	for(j=0;j<multiplex[i].nlinks;j++)
	// 	{
	// 		printf("%llu %llu\n", multiplex[i].links_table[j].node_i,multiplex[i].links_table[j].node_j);
	// 	}
	// 	printf("\n");
	// }
	//***************************************************
	evaluate_KC(multiplex,M,MAX_NODE,primeslist,nodes_unique,length_unique,ITER,parallel,save_prime_weight,prime_weight_filename,verbose);
	free_structure(primeslist,nodes_unique,multiplex,M);

}



unsigned long long int *load_unique_nodes(char *filename, unsigned long long int *length_unique, unsigned long long int *MAX_NODE,unsigned int verbose)
{
  FILE *fp;
  unsigned long long int *unique_nodes, size=10000, k=0;
  int scan=0;
  fp = fopen(filename,"r");
  unique_nodes=(unsigned long long int *)malloc(size*sizeof(unsigned long long int));
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
	char line[512];
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


    if (fgets(line,512, fp) == NULL) break;
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
	array= (struct layer *)realloc(array,(*M)*sizeof(struct layer));
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
  		mpz_nextprime(prime, prime); //Set prime to the next prime greater than prime.
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


void evaluate_KC(struct layer *multiplex, unsigned int M,unsigned long long int MAX_NODE, unsigned long long int *primeslist, unsigned long long int *nodes_unique, unsigned long long int length_unique, unsigned int ITER,unsigned int parallel, unsigned int save_prime_weight, char *prime_weight_filename, unsigned int verbose)
{
	
	unsigned long long int dimension_KC,dimension_KC_aggregate,k,nedges_intersection;
	unsigned int nthread;
	struct layer **temporarystruct;
	void **nthread_root;
	struct edge *elementptr;
	double sumKC,sumKCaggregate,complexity,numerator,denominator,complexity_cumulative=0.0;

	
	//***********Evaluate the Kolmogorov Complexity of the whole multiplex************
	sumKC=0;
	sumKCaggregate=0;
	nedges_intersection=0;
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
    	//mpz_set_ui(maximum_prime_product[nthread],primeslist[M-1]);
			temporarystruct[nthread]=multiplex_into_binarytree(&(nthread_root[nthread]),multiplex,M,MAX_NODE,primeslist,nodes_unique,length_unique,&nedges_intersection,nthread);
			// //Printing the whole multiplex  & the corresponding weighted aggregate
			// if (maximum_overlap[nthread] == 0)
				// maximum_overlap[nthread]=1;
			//gmp_fprintf(stderr,"prime in the aggregate: %Zd\n",maximum_prime_product[nthread]);
			prime_layer=primeslist[M-1];
			dimension_KC=print_file_bothfast_onstring(nthread_root[nthread],&dimension_KC_aggregate,nthread);
			sumKC+=dimension_KC;
			sumKCaggregate+=dimension_KC_aggregate;
			complexity_cumulative+=(1.0*dimension_KC)/(1.0*dimension_KC_aggregate);
			tdestroy(nthread_root[nthread]);
			free_structgmp_multiplex(temporarystruct[nthread],multiplex,M);
		}
	//numerator=(1.0*sumKC)/(1.0*ITER);
	//denominator=  (1.0*sumKCaggregate)/(1.0*ITER);
	//complexity = numerator/denominator;
	complexity = complexity_cumulative/ITER;
	if(verbose)
		fprintf(stderr,"#KC_multiplex\tNLinks_multiplex\tKC_aggregate\tNlinks_aggregate\tComplexity\n");
	printf("%lf\t%llu\t%lf\t%llu\t%lf\n", (1.0*sumKC)/(1.0*ITER),Nlinks_multiplex[0], (1.0*sumKCaggregate)/(1.0*ITER),Nlinks_aggregate[0],complexity);
	
	if (save_prime_weight == 1)
	{
		nthread= 0;
		nthread_root[nthread]=NULL;
		temporarystruct[nthread]=multiplex_into_binarytree_saving(&(nthread_root[nthread]),multiplex,M,MAX_NODE,primeslist,nodes_unique,length_unique,&nedges_intersection,nthread);
		fp_save_prime_weight = fopen(prime_weight_filename,"w+");
		twalk(nthread_root[nthread],print_elements_onfile_saving);
		fclose(fp_save_prime_weight);
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

struct layer *multiplex_into_binarytree(void **root,struct layer *multiplex, unsigned int M, unsigned long long int MAX_NODE,
 unsigned long long int *primeslist, unsigned long long int *nodes_unique, unsigned long long int length_unique, unsigned long long int  *length_nedges_intersection, unsigned int nthread)
{
 
	unsigned long long int i,j,k, nodei,nodej,*shuffle;
	struct edge *update;
	unsigned long long int *distr1;
	struct layer *temporarystruct;
	distr1=particular_primes_sequence(multiplex, M);
	shuffle=shuffle_nodes(length_unique,nodes_unique);
	temporarystruct = update_nodelabel(multiplex,shuffle,M,length_unique,nodes_unique);
	//mpz_init(maximum_prime_product[nthread]);
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
	//fprintf(stderr,"%d ", length_nedges_intersection );
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


int mpz_cmp_flag_to_01(int mpz_cmp_flag){
  if(mpz_cmp_flag <= 0)
    return 0;
  else
    return 1;
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

void free_structgmp_multiplex(struct layer *temporarystruct, struct layer *multiplex, unsigned int M)
{
	unsigned int i;
	unsigned long long int j;
	for (i = 0; i <M; i++)
	{
		
		for (j = 0; j < multiplex[i].nlinks; j++)
		{
			mpz_clear(temporarystruct[i].links_table[j].weight_link);
		}
		free(temporarystruct[i].links_table);	
	}

	free(temporarystruct);	
}


unsigned long long int print_file_bothfast_onstring(void **root,unsigned long long int *dimension_aggregate,  unsigned int thread_number)
{
	unsigned long long int dimension=0;
	dimension_string[thread_number]=100000;
	fp_string[thread_number]= (char *) malloc (dimension_string[thread_number] * sizeof(char));
	fp_string_aggregate[thread_number]= (char *)malloc (dimension_string[thread_number] * sizeof(char));
	string_counter[thread_number]=0;
	string_counter_aggregate[thread_number]=0;

	//Printing the weighted multiplex on string 
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
		string_counter[thread_number] += gmp_sprintf(fp_string[thread_number]+string_counter[thread_number],"%llu %llu %Zd\n", e->node_i, e->node_j,e->weight_link);
		mpz_init(temporary_mpz);

		//mpz_ui_pow_ui(temporary_mpz,prime_layer,maximum_overlap[thread_number]);
		mpz_set(temporary_mpz,maximum_prime_product[thread_number]);
		// gmp_printf("%u %u %Zd --- prime_layer: %u -- maximum_overlap: %u\n", e->node_i, e->node_j,temporary_mpz,prime_layer,maximum_overlap[thread_number]);
		string_counter_aggregate[thread_number] += gmp_sprintf(fp_string_aggregate[thread_number]+string_counter_aggregate[thread_number],"%llu %llu %Zd\n", e->node_i, e->node_j,temporary_mpz);
		//string_counter_aggregate[thread_number] += gmp_sprintf(fp_string_aggregate[thread_number]+string_counter_aggregate[thread_number],"%u %u %u\n", e->node_i, e->node_j,e->overlap_aggregate);
		//string_counter_aggregate[thread_number] += gmp_sprintf(fp_string_aggregate[thread_number]+string_counter_aggregate[thread_number],"%u %u %u\n", e->node_i, e->node_j,2);
		//gmp_printf("%Ld %Ld %.0lf\n", e->node_i, e->node_j,pow(prime_layer,(double)e->overlap_aggregate));
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

 void chopN_string(char *str, size_t n)
{
    assert(n != 0 && str != 0);
    size_t len = strlen(str);
    if (n > len)
        return; 
    memmove(str, str+n, len - n + 1);
}




struct layer *multiplex_into_binarytree_saving(void **root,struct layer *multiplex, unsigned int M, unsigned long long int MAX_NODE,
 unsigned long long int *primeslist, unsigned long long int *nodes_unique, unsigned long long int length_unique, unsigned long long int  *length_nedges_intersection, unsigned int nthread)
{
 
	unsigned long long int i,j,k, nodei,nodej;
	struct edge *update;
	unsigned long long int *distr1;
	struct layer *temporarystruct;
	distr1=particular_primes_sequence(multiplex, M);
	temporarystruct = update_nodelabel_saving(multiplex,M,length_unique,nodes_unique);
	//mpz_init(maximum_prime_product[nthread]);
	//mpz_set_ui(maximum_prime_product[nthread],primeslist[M-1]);
	
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
				//if (mpz_cmp_flag_to_01(mpz_cmp((*update).weight_link,maximum_prime_product[nthread])))
        {
          //gmp_printf("changing %Zd in %Zd\n",maximum_prime_product[nthread],(*update).weight_link);
          //mpz_set(maximum_prime_product[nthread],(*update).weight_link);  
        }
			}
		}
	}
	//fprintf(stderr,"%d ", length_nedges_intersection );
	free(distr1);
	return(temporarystruct);
}



struct layer *update_nodelabel_saving(struct layer *multiplex, unsigned int M,unsigned long long int length_unique, unsigned long long int *nodes_unique)
{
	unsigned int i;
	unsigned long long int nodei,nodej,j;
	struct layer *new_duplex= (struct layer *) malloc(M * sizeof(struct layer));
	for (i=0; i<M;i++)
	{
		new_duplex[i].nlinks=multiplex[i].nlinks;
		new_duplex[i].max_node=nodes_unique[multiplex[i].max_node];
		new_duplex[i].links_table = (struct edge *) malloc (multiplex[i].nlinks * sizeof(struct edge));
		for(j=0;j<multiplex[i].nlinks;j++)
		{
			nodei=multiplex[i].links_table[j].node_i;
			new_duplex[i].links_table[j].node_i=nodes_unique[nodei];

			nodej=multiplex[i].links_table[j].node_j;
			new_duplex[i].links_table[j].node_j=nodes_unique[nodej];
			mpz_init(new_duplex[i].links_table[j].weight_link);
			new_duplex[i].links_table[j].overlap_aggregate=0;
			//printf("%lld %lld\n", nodei,nodej);

		}
	}

	return(new_duplex);
}


void print_elements_onfile_saving(const void *nodep, const VISIT which, const int depth)
{
	
	struct edge *e = *((struct edge **) nodep);

	switch (which) {
	case leaf:
	case postorder:
		length++;
		unsigned int thread_number= omp_get_thread_num();
		mpz_t temporary_mpz;
		mpz_init(temporary_mpz);
		mpz_set(temporary_mpz,maximum_prime_product[thread_number]);
		gmp_fprintf(fp_save_prime_weight,"%llu %llu %Zd\n", e->node_i, e->node_j,e->weight_link);
		mpz_clear(temporary_mpz);
		break;
	default:
		break;
	}
}