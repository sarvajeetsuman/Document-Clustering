//Cure Algorithm Used for comparison 
#include<iostream>
#include<limits>
#include<stxxl/vector>
#include<stxxl/stream>
#include<stxxl/bits/containers/matrix.h>
#include<fstream>
#include<string>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include<unistd.h>


using namespace stxxl;
using namespace std;

typedef double value_type ;


stxxl::unsigned_type block =64;//32
stxxl::unsigned_type threshold=7,minPts=10;


typedef stxxl::VECTOR_GENERATOR<unsigned int, 1, 1, 64*1024*1024, stxxl::RC, stxxl::lru>::result vector_type;
vector_type my_vector;



bool isInCorePoints(int element)
	{
		vector_type::const_iterator iter = my_vector.begin();
		for (stxxl::uint64 j = 0; j < my_vector.size(); j++)
			{
				if(*(iter)==element)
					return true;
				//for(int z=0;z<k;z++)
				//{
					++iter;
					//j++;
				//}
				
			}
		return false;
	}
int similarity(double v1[],double v2[],int size )
	{ 
		
		int flag=0,x=0,y=0,n=0;
		for( x=0;x<size && flag!=1;x++)
			{
				if(v1[0]==v2[x])
					flag=1;
			}
		for( x=0;x<size && flag!=2;x++)
			{
				if(v2[0]==v1[x])
					flag=flag+1;
			}
		if(flag==2)
			{
				for(x=0;x<size;x++)
					for(y=0;y<size;y++)
					{
						if(v1[x]==v2[y])
							n++;
						
					}
			}
			return n;
	}
		
double distance2(double a[],double b[],int width)
	{
		double sum=0;
		for(int i=0;i<width;i++)
		{
			sum+=(a[i]-b[i])*(a[i]-b[i]);
		}
		return sqrt(sum);
	}
		
 		

		char file[]="/media/livingroom/stxxl/stxxl-1.4.1/local/files.txt";
		ifstream doc;
		double idf;
		
		
		int size=21578; //number of document 
		
		int dimension=64;//number of term in doucment
		
		int k=7;
		
		int init=size*0.1;
		
		
		const int large_block_order=64;
		stxxl::unsigned_type no_of_vector=size;
		int rows=floor(large_block_order*large_block_order/dimension);
	
		const int small_block_order=64;
		int no_of_blocks=ceil(size/rows);

		
		stats_data stats_begin, stats_end, stats_cumulative;
	
		int internal_memory1 = 2 * sizeof(value_type)*large_block_order*large_block_order;
		int internal_memory2 = 2* sizeof(value_type)*small_block_order*small_block_order;
	
		typedef block_scheduler< matrix_swappable_block<value_type, small_block_order> > bst;
        typedef matrix<value_type, small_block_order> mt;
        
        typedef block_scheduler< matrix_swappable_block<value_type, large_block_order> > bst2;
        typedef matrix<value_type, large_block_order> mt2;
	 	
		bst * b_s=new bst(internal_memory2);
		bst2 * b_s2=new bst2(internal_memory1);
		bst &bs=*b_s;
		bst2 &bs2=*b_s2;
	
	// Create 2 matrices with given dimensions
		mt2 *a=new mt2(bs2, no_of_vector, dimension);
		mt *b=new mt(bs,no_of_vector,k+1);
		
	
		typedef mt2::row_major_iterator row_iterator;
		typedef mt::row_major_iterator row_iterator2;
		
		
#ifndef FLOAT
#define FLOAT	float			/* float type used throughout */
#endif

#define FLOAT_FORMAT	"%.15lg"	/* printf format use for floats */
#define FLOAT2_FORMAT	"%.15lg %.15lg"	/* printf format use for floats */


#define new(type) \
	(type *)calloc(1, sizeof(type))

#define new_array_of(n, type) \
	(type *)calloc((unsigned)(n), sizeof(type))

#define change_array_size(array, n, type) \
	(type *)realloc(array, (unsigned)(n) * sizeof(type))

#define new_2d_array_of(n1, n2, type) \
	(type **)calloc_2d((unsigned)(n1), (unsigned)(n2), sizeof(type))

#define new_3d_array_of(n1, n2, n3, type) \
	(type ***)calloc_3d((unsigned)(n1), (unsigned)(n2), (unsigned)(n3), \
								sizeof(type))

#define change_2d_array_size(array, n1, n2, m1, m2, type) \
	(type **)realloc_2d(array, (unsigned)(n1),(unsigned)(n2), \
				(unsigned)(m1), (unsigned)(m2), sizeof(type)) 


#define FirstPruneRatio 0.3333
#define SecondPruneMulti 3
#define FirstCut 3
#define SecondCut 10

#define MY_OK	1		/* avoid conflict with curses result code */
#define MY_ERR	0

/* error handling macros */
#define Erreturn(msg) { \
	(void)strcpy(ERR_MSG, msg); \
	return MY_ERR; \
    }

#define Erreturn1(msg, x) { \
	(void)sprintf(ERR_MSG, msg, x); \
	return MY_ERR; \
    }

#define Erreturn2(msg, x, y) { \
	(void)sprintf(ERR_MSG, msg, x, y); \
	return MY_ERR; \
    }

#define Erreturn3(msg, x, y, z) { \
	(void)sprintf(ERR_MSG, msg, x, y, z); \
	return MY_ERR; \
    }

#define IfErr(x) 	if ((x) == MY_ERR)
#define IfEOF(x)	if ((x) == EOF)

/*
 * macros handling D/C values
 */
#ifndef NO_DONTCARES
/*
 * we use IEEE NaN to represent don't care values -- ugly, but it works
 */
static long _nan = 0x7fffffff;
#define DC_VAL		((double)*(float *)(&_nan))
#define IS_DC(x)	isnan(x)

#ifdef sgi
#define isnan(x)	((x) == DC_VAL)
#endif

#endif /* !NO_DONTCARES */


/*
 * error.c -- Error message handling for cluster
 */
 


#define Nmsg 200		/* error msg can be fairly long */
char ERR_MSG[Nmsg];
int  ERR_FLAG;
/*
 * cluster.h -- declarations for clustering routines
 */
 

typedef struct _tree {
    FLOAT  *pat;
    FLOAT *mean;
    FLOAT **rep;
    int     pruned;
    int     size;
    int     root;
    int     leaf;
    FLOAT   y;
    FLOAT   distance;
    struct _tree *r_tree, *l_tree;
}       BiTree;

extern BiTree *new_tree();

#define LEAF -1

#ifndef TRUE
#define TRUE	1
#endif
#ifndef FALSE
#define FALSE	0
#endif


#define NONE (-2)

#define BUFSIZE 256

#ifndef SCALE
#define SCALE "_SCALE_"
#endif

#ifndef DONTCARE
#define DONTCARE "D/C"
#endif

#ifndef MAXFLOAT
#define     MAXFLOAT        ((float)3.40282346638528860e+38)
#endif
/*
static FLOAT   distance();
static FLOAT   root();
static FLOAT   cure_distance();
static void    merge();
static void    repInTree();
static void print_rep_points(); */

static char buffer[BUFSIZE];	/* temporary string read buffer */

#define PCA	"pca"
char   *program;		/* program name */

int     vflag = 0;		/* print explanatory messages */
static int pflag = 0;		/* do PCA instead of HCA */
static int dflag = 0;		/* print distances */
static int tflag = 0;		/* print tree in ASCII */
static int Tflag = 0;		/* display tree using curses */
static int gflag = 0;		/* print tree in graph(1) format */
static int bflag = 0;		/* print tree in graph(1) format, with breaks */
static int Bflag = 0;		/* print patterns as bit vectors */
static int sflag = 0;		/* suppress scaling  */
int	Eflag = 0;		/* print eigenvalues */

static int width = 0;		/* width of ASCII tree representation */

static int normp = 2;		/* which l-norm to use for distance metric */

static int csize = 0; 		/* number of clusters */
static int rsize = 10; 		/* number of rep points in a cluster */
static float alpha = 0.3;	/* shrinking factor alpha */
static FLOAT root(FLOAT dist);
char ***calloc_3d(unsigned N1, unsigned N2, unsigned N3, unsigned size);
void repInTree(BiTree *item, int rcount,FLOAT *mean, FLOAT **rep, int lpat, FLOAT *maxDist, FLOAT **repPoint);
void merge(BiTree *item, int pair1, int pair2, int lpat, float min_dist);
void print_rep_points( BiTree *item, int cno, int lpat);
void cluster(FLOAT **pattern, char  **name,int lpat,int  npat,int  *pinfo);
int nstrings(FILE *fp, FLOAT **patternP, char  **nameP);
char  **calloc_2d(unsigned N1, unsigned N2, unsigned size);
void free_3d_array(char ***array, int N1, int N2);
char  **realloc_2d(char **array, unsigned N1, unsigned N2,unsigned M1,unsigned M2,unsigned size);
void free_2d_array(char  **array,int N);
int read_pattern(FILE   *Pfile, FLOAT *** patternP, int *lpatternP, int *npatternP, char ***nameP);
int read_names(FILE   *fp, char  **name,  int npat);
char *new_string(char *string);
void print_names(BiTree *tree, char **name);
static FLOAT distance(FLOAT *pat1, FLOAT *pat2, int lpat);
static FLOAT cure_distance(BiTree *x, BiTree *y,int lpat);
void find_nnb(BiTree items, int which,  int lpat,int *index,FLOAT *ndist);
static void usage()
{
	cout<<"\n In usage Method";
    if (!pflag)
	fprintf(stderr,"usage: %s [-dtTgbBsv] [-w width] [-n norm] [-k no-cluster] [-r rep-size] [-a alpha] [vectorfile [namesfile]]\n", program);
    else
	fprintf(stderr, "usage: %s [-Esv] [-e eigenbase] [-c pcs] [vectorfile [namesfile]]\n", PCA);
    exit(2);
}

int main(int argc, char *argv[])
{

	STXXL_MSG("Restitution\n");
	stats_data stats_before, stats_after;
	matrix_operation_statistic_data matrix_stats_before, matrix_stats_after;
	matrix_stats_before.set();
	stats_before = *stats::get_instance();
	stxxl::stats_data begin_1(*stats::get_instance());
	cout<<"Start the operation\n\n";
	
	
    FILE   *fp, *pfile;
    char   *efile = NULL;
    char   *comps = NULL;
    FLOAT **pattern = NULL;
    int     lpat, npat;
    char  **name = NULL;

    int     opt;
    extern char *optarg;
    extern int optind;

    int *pinfo, i;
    char fname[80];

    /* who are we ? */
    if ((program = strrchr(argv[0], '/')) == NULL)
	program = argv[0];
    else
	program += 1;

    if (strcmp(program, PCA) == 0)
	pflag = 1;

    while ((opt = getopt(argc, argv, "pdtTgbBvsn:a:r:k:e:c:w:E")) != -1)
	switch (opt) {
	cout<<"\n in main method and switch() opt is: "<<opt<<endl;
	case 'p':
	    pflag = 1;
	    break;
	case 'd':
	    if (pflag)
		usage();
	    dflag = 1;
	    break;
	case 't':
	    if (pflag)
		usage();
	    tflag = 1;
	    break;
	case 'T':
	    if (pflag)
		usage();
	    Tflag = 1;
	    break;
	case 'g':
	    if (pflag)
		usage();
	    gflag = 1;
	    break;
	case 'b':
	    if (pflag)
		usage();
	    bflag = 1;
	    gflag = 1;
	    break;
	case 'B':
	    if (pflag)
		usage();
	    Bflag = 1;
	    break;
	case 'v':
	    vflag = 1;
	    break;
	case 'n':
	    if (pflag)
		usage();
	    if (sscanf(optarg, "%d", &normp) != 1 || normp < 0)
		usage();
	    break;
	case 'k':
	    if (pflag)
		usage();
	    if (sscanf(optarg, "%d", &csize) != 1 || csize < 0)
		usage();
	    break;
	case 'a':
	    if (pflag)
		usage();
	    if (sscanf(optarg, "%f", &alpha) != 1 || alpha < 0 || alpha > 1.0)
		usage();
	    break;
	case 'r':
	    if (pflag)
		usage();
	    if (sscanf(optarg, "%d", &rsize) != 1 || rsize < 0)
		usage();
	    break;
	case 'e':
	    if (!pflag)
		usage();
	    efile = optarg;
	    break;
	case 'c':
	    if (!pflag)
		usage();
	    comps = optarg;
	    break;
	case 's':
	    sflag = 1;
	    break;
	case 'w':
	    if (pflag)
		usage();
	    if (sscanf(optarg, "%d", &width) != 1 || width <= 0)
		usage();
	    break;
	case 'E':
	    if (!pflag)
		usage();
	    Eflag = 1;
	    break;
	case '?':
	    usage();
	    break;
	
	}
	cout<<"\nout side of switch case"<<endl;
    if (!(pflag || dflag || tflag || Tflag || gflag || Bflag))
	dflag = tflag = vflag = 1;	/* default behavior */

    if (optind + 2 < argc)
	usage();

    if (!(optind < argc) || !strcmp(argv[optind], "-"))
	fp = stdin;
    else
	IfErr(fp = fopen(argv[optind], "r")) {
	    fprintf(stderr, "%s: cannot open file %s\n", argv[0], argv[optind]);
	    exit(1);
	}
	cout<<"\n argv[optiond] "<<argv[optind]<<endl;//
    IfErr(read_pattern(fp, &pattern, &lpat, &npat, &name)) {
	fprintf(stderr, "%s: %s: cannot read pattern\n", program, ERR_MSG);
	exit(1);
    }
    if (vflag)
	fprintf(stderr, "read %d patterns:  size = %d\n", npat, lpat);

    if (optind + 1 < argc) {
	IfErr (name = new_array_of(npat, char *)) {;
	    fprintf(stderr, "%s: not enough core for name array\n", program);
	    exit(1);
	}

	if (!strcmp(argv[optind + 1], "-"))
	    fp = stdin;
	else
	    IfErr(fp = fopen(argv[optind + 1], "r")) {
		fprintf(stderr, "%s: cannot open file %s\n",
			program, argv[optind + 1]);
		exit(1);
	    }

	IfErr(read_names(fp, name, npat)) {
	    fprintf(stderr, "%s: %s: cannot read names\n", program, ERR_MSG);
	    exit(1);
	}
    }

    IfErr(pinfo = new_array_of(npat, int)) {
            fprintf(stderr, "%s: not enough core\n", program);
            exit(1);
    }

    sprintf(fname, "%s-partition", argv[optind]);
    IfErr(pfile = fopen(fname, "w")) {
            fprintf(stderr, "%s: cannot open file %s\n", argv[0], argv[optind]);
            exit(1);
    }
	cout<<"\n above cluster";
    cluster(pattern, name, lpat, npat, pinfo);

    for (i=0; i<npat; ++i)
        fprintf(pfile, "%d\n", pinfo[i]);
    free(pinfo);
    fclose(pfile);
    
    
	cout<<"operation is over\n\n";
	stxxl::stats_data begin_2(*stats::get_instance());
	STXXL_MSG(begin_2-begin_1);
	stats_after = *stats::get_instance();
	matrix_stats_after.set();
	cout<<"Now History\n\n";
	STXXL_MSG(stats_after - stats_before);
	STXXL_MSG("END OF THE PROGRAM\n");
	
	return 1;

   
}

/* skip blanks and next end of line */
char skip_blanks(FILE *fp)
{
	
    char    c;

    while ((c = getc(fp)) == ' ' || c == '\t');
    if (c != '\n')
	ungetc(c, fp);
    return c;
}

int read_names(FILE *fp, char **name, int npat)
{
	cout<<"\n In read_names method\n";
    register int i;

    for (i = 0; i < npat; i++) {
	char    c = skip_blanks(fp);

	if (c == '\"') {
	    getc(fp);
	    fgets(buffer, sizeof(buffer), fp);
	    buffer[strlen(buffer) - 1] = '\0';
	}
	else {
	    IfEOF(fscanf(fp, "%s", buffer))
		Erreturn("not enough names");
	    skip_blanks(fp);
	}
	//cout<<"\n In read_names method buffer="<<buffer<<endl;
	IfErr(name[i] = new_string(buffer))
	    Erreturn("not enough core");
    }
    return MY_OK;
}

void print_names(BiTree *tree, char **name)
{
    if (tree->leaf != LEAF) {
	if (name)
	    printf(" %s", name[tree->leaf]);
	else
	    printf(" %d", tree->leaf);
    }
    else {
	print_names(tree->r_tree, name);
	print_names(tree->l_tree, name);
    }
}

void print_k_cluster(BiTree *tree, char **name,int *pinfo, int cno)
{
    if (tree->leaf != LEAF) {
/*
        if (name)
            printf("%d : %s\n", tree->leaf, name[tree->leaf]);
        else
            printf("%d : %d\n", tree->leaf, tree->leaf);
*/

        pinfo[tree->leaf] = cno;
    }
    else {
        print_k_cluster(tree->r_tree, name, pinfo, cno);
        print_k_cluster(tree->l_tree, name, pinfo, cno);
    }
}


void find_nnb(BiTree *items, int which,  int lpat,int *index,FLOAT *ndist)
{
	//cout<<"\n In find_nnb method"<<endl;
    int i;
    FLOAT dist, min_dist;
    int min_index;

    if (items[which].root == FALSE) {
	*index = NONE;
	return;
    }

    min_index = NONE;
    min_dist = 0.0;
    /*
     * find minimum distance neighbor -- to avoid duplication
     * only pairs with 1st index < 2nd index are considered
     */
    for (i = 0; i < which; i++) {
	if (items[i].root == FALSE || items[i].pruned == TRUE)
	    continue;
	
	dist = cure_distance(&items[which], &items[i], lpat);
	if (min_index == NONE || dist < min_dist) {
	    min_index = i;
	    min_dist = dist;
	}
    }

    *index = min_index;
    if (min_index >= 0)
	*ndist = min_dist;
    return;
}



void cluster(FLOAT **pattern, char **name, int lpat, int npat, int *pinfo)
{
	cout<<"\n In Cluster method \n";
    register int i, j, k, nodes = npat, cno;
    BiTree *item = new_array_of(npat, BiTree);
    /*
     * for each data point or cluster center, we keep the index of the nearest
     * neighbor, as well as the distance to it.
     */
    int *nnb_index = new_array_of(npat, int);
    FLOAT *nnb_dist = new_array_of(npat, FLOAT);

    if (item == NULL || nnb_index == NULL || nnb_dist == NULL) {
	    fprintf(stderr, "%s: not enough core\n", program);
	    exit(1);
    }

    /*
     * initialize leaf nodes
     */
    for (i = 0; i < npat; i++) {
	item[i].pat = pattern[i];
	IfErr (item[i].mean = new_array_of(lpat, FLOAT)) {
	    fprintf(stderr, "%s: not enough core\n", program);
	    exit(1);
	}
	memcpy(item[i].mean, item[i].pat, lpat*sizeof(FLOAT));
	IfErr (item[i].rep = new_array_of(rsize, FLOAT *)) {
	    fprintf(stderr, "%s: not enough core\n", program);
	    exit(1);
	}
	for (j=0; j < rsize; ++j) {
		IfErr (item[i].rep[j] = new_array_of(lpat, FLOAT)) {
	    		fprintf(stderr, "%s: not enough core\n", program);
	    		exit(1);
		}
	}
	memcpy(item[i].rep[0], item[i].pat, lpat*sizeof(FLOAT));
	item[i].root = TRUE;
	item[i].size = 1;
	item[i].leaf = i;
	item[i].distance = 0.0;
	item[i].l_tree = item[i].r_tree = NULL;
    }

    /*
     * initialize nearest neighbors
     */
    for (i = 0; i < npat; i++) {
	find_nnb(item, i, lpat, &nnb_index[i], &nnb_dist[i]);
/*
printf("%d : %d %f %f %f\n", i, nnb_index[i], nnb_dist[i], item[i].rep[0][0], item[i].rep[0][1]);
*/
    }

    /*
     * cluster until done
     */
  
    while (nodes > csize) {
	BiTree  *newitem;
	FLOAT   dist, min_dist;
	int     pair1, pair2;
	
	/*
	 * find minimum distance pair
	 */
	
	pair1 = NONE;
	min_dist = 0.0;
	for (i = 0; i < npat; i++) {
	    if (item[i].root == FALSE || item[i].pruned == TRUE)
		continue;
	    if (nnb_index[i] != NONE &&
		(pair1 == NONE || nnb_dist[i] < min_dist)) {
		pair1 = i;
		pair2 = nnb_index[i];
		min_dist = nnb_dist[i];
	    }
	}
	if (pair1 == NONE)
	    break;		/* analysis finished */

	min_dist = root(min_dist);

	if (dflag) {
	    //printf("minimum distance = %f\t(", (float)min_dist);
	    //print_names (&item[pair1], name);
	    //printf(" )\t(");
	    //print_names (&item[pair2], name);
	    //printf(" )\n");
	}

        merge(item, pair1, pair2, lpat, min_dist);

	/*
	 * update nearest neighbors
	 */
	for (i = 0; i < npat; i++) {
	    if (nnb_index[i] == NONE || item[i].root == FALSE || 
	        item[i].pruned == TRUE) {
		continue;
	    } else if (nnb_index[i] == pair1 || nnb_index[i] == pair2) {
		/*
		 * worst case: the old nnb is the node that will disappear.
		 * recompute nnb from scratch
		 */
		find_nnb(item, i, lpat, &nnb_index[i], &nnb_dist[i]);
	    } else if (pair1 < i &&
		       (dist = cure_distance(&item[pair1], &item[i], lpat))
			    < nnb_dist[i]) {
		/*
		 * distance to new node is smaller than previous nnb,
		 * so make it the new nnb.
		 */
		nnb_index[i] = pair1;
		nnb_dist[i] = dist;
	    }
	}

	/* number of nodes have been reduced */
	--nodes;

	/* prune clusters based on their size */
	
	if (nodes == (int)(npat * FirstPruneRatio)) {
		printf("==== First phase of pruning at %d nodes remaining ====\n", nodes);
	    for (i = 0; i < npat; i++) {
	    //cout<<"\n npat="<<npat<<"\t i="<<i<<endl;
		if (item[i].root == TRUE && item[i].size < FirstCut) {
		//printf("====      %d of size %d removed\n", i, item[i].size);
			item[i].pruned = TRUE;
			--nodes;
		}
	    }

	    /* recalc min dist */
	   
	    for (i = 0; i < npat; i++) {
	    
	    if(nnb_index[i]>=0){
	    
		if (item[i].root == TRUE && item[nnb_index[i]].pruned == TRUE)
			find_nnb(item, i, lpat, &nnb_index[i], &nnb_dist[i]);//cout<<"\ncome line 682 i="<<i;
	    }
	    }
	}
	else if (nodes == csize * SecondPruneMulti) {
		printf("==== Second phase of pruning at %d nodes remaining ====\n", nodes);
	    for (i = 0; i < npat; i++) {
		if (item[i].root == TRUE && item[i].pruned == FALSE &&
		    item[i].size < SecondCut) {
		//printf("====      %d of size %d removed\n", i, item[i].size);
			item[i].pruned = TRUE;
			--nodes;
		}
	    }

	    /* recalc min dist */
	    cout<<"\n recalc min dist \n";
	    for (i = 0; i < npat; i++) {
	    if(nnb_index[i]>=0){
		if (item[i].root == TRUE && item[nnb_index[i]].pruned == TRUE)
			find_nnb(item, i, lpat, &nnb_index[i], &nnb_dist[i]);
	    }
	    }
	}
    }
	
#if 0
    /*
     * search for root
     */
    j = 0;
    for (i = 0; i < npat; i++) {
      if (item[i].root == TRUE) {

        if (tflag) {
	  if (vflag)
	    printf("Cluster %d Tree = \n", j);
	   IfErr (print_tree(&item[i], name, npat, width)) {
	    fprintf(stderr, "%s: %s: error printing tree\n", program, ERR_MSG);
	    exit(1);
	  }
        }

        if (gflag) {
	  if (vflag)
	    printf("Cluster %d Tree Graph= \n", j);
	  graph_tree(&item[i], name, npat, bflag);
        }

        if (Tflag) {
#ifdef HAVE_CURSES
	  if (vflag)
	    printf("Cluster %d Displaying tree with curses...\n", j);
	  curses_tree(&item[i], name, npat, width);
#else
	  fprintf(stderr, "Sorry, no curses support available.\n");
#endif
        }

        if (Bflag) {
	    print_bits(&item[i], name, npat);
        }

	printf("\n");
	++j;
      }
    }
#endif
	cout<<"\nexplicit size is given\n";
    if (csize > 0) {  /* explicit size is given */
	cno = 0;
	for (i=0; i < npat; ++i)
		if (item[i].root == TRUE) {
/*
			printf("%d **** Cluster %d from Hierarchical Clustering ****\n", 
                                item[i].size, cno);
*/
			if (item[i].pruned == FALSE) {
				print_k_cluster(&item[i], name, pinfo, cno);
				print_rep_points(&item[i], cno, lpat);
				++cno;
			}
			else {
				print_k_cluster(&item[i], name, pinfo, nodes);
			}
		}
/*
	printf("-1 ****\n");
*/
    }
}

static FLOAT cure_distance(BiTree *x, BiTree *y,int lpat)
{
	//cout<<"\n In cure_distance method\n";
    register int i, j;
    FLOAT min_dist = MAXFLOAT, dist;
    int xmax, ymax;

    if (x->size > rsize)
	xmax = rsize;
    else
	xmax = x->size;

    if (y->size > rsize)
	ymax = rsize;
    else
	ymax = y->size;

    for (i=0; i<xmax; ++i) 
	for (j=0; j<ymax; ++j)
		if ((dist = distance(x->rep[i], y->rep[j], lpat)) < min_dist)
			min_dist = dist;

    return min_dist;
}

static FLOAT distance(FLOAT *pat1, FLOAT *pat2, int lpat)
{
	
    register int i;
    FLOAT   dist = 0.0;

    for (i = 0; i < lpat; i++) {
	FLOAT   diff = 0.0;

#ifndef NO_DONTCARES
	if (!IS_DC(pat1[i]) && !IS_DC(pat2[i]))
#endif
	    diff = pat1[i] - pat2[i];
	
	switch (normp) {
	    FLOAT   adiff;

	case 2:
	    dist += diff * diff;
	    break;
	case 1:
	    dist += fabs(diff);
	    break;
	case 0:
	    if ((adiff = fabs(diff)) > dist)
		dist = adiff;
	    break;
	default:
	    dist += pow(fabs(diff), (double) normp);
	    break;
	}
    }
    return dist;
}

static FLOAT root(FLOAT dist)
{
	//cout<<"\n In root method\n";
    switch (normp) {
    case 2:
	return sqrt(dist);
    case 1:
    case 0:
	return dist;
    default:
	return pow(dist, 1 / (double) normp);
    }
}

BiTree *new_tree(BiTree *item)
{
    BiTree *tree;

    IfErr (tree = new(BiTree))
	return NULL;

    tree->r_tree = item->r_tree;
    tree->l_tree = item->l_tree;
    tree->leaf = item->leaf;
    tree->root = item->root;
    tree->size = item->size;
    tree->distance = item->distance;
    tree->pat = item->pat;

    return tree;
}

#define Blksize 128

int read_pattern(FILE *Pfile, FLOAT ***patternP, int *lpatternP, int *npatternP, char ***nameP)
{
	cout<<"\n In read_pattern\n";
    register int i;
    int     status;
    int     Asize = Blksize;	/* current array size */

    /***** these local variables are temporary storage and are ****
     ***** copied to the real variables if there is no error ******/
    int     lpattern, npattern;	/* size and # of patterns */
    FLOAT **pattern;
    FLOAT  *first_pattern;
    char  **name = NULL;
    char   *first_name = NULL;
    FLOAT  *scales = NULL;

    IfErr(lpattern = nstrings(Pfile, &first_pattern, &first_name))
	return MY_ERR;

    /* check for scaling info */
    if (first_name != NULL && strcmp(first_name, SCALE) == 0) {
	scales = first_pattern;
	free(first_name);

	/* read first vector again */
	IfErr(status = nstrings(Pfile, &first_pattern, &first_name))
	    return MY_ERR;
	if (status != lpattern)
	    Erreturn("scaling vector not matching data");
    }

    /* allocate space for input/target arrays */
    IfErr(pattern = new_2d_array_of(Asize, lpattern, FLOAT))
	Erreturn("cannot allocate memory for patterns");

    if (first_name != NULL)
	IfErr(name = new_array_of(Asize, char *))
	    Erreturn("cannot allocate memory for names");

    /**** this loop reads in one line from pattern file,**
     **** stores each pattern into pattern buffer    *****/
	//cout<<"\n pattern="<<pattern<<endl;
    for (npattern = 0;; npattern++) {
		//cout<<"\n npattern="<<npattern<<endl;
	if (npattern >= Asize) {/* need to allocate more space for arrays */
	    IfErr(pattern = change_2d_array_size(
		(char **)pattern, Asize, lpattern, Asize + Blksize, lpattern, FLOAT))
		Erreturn("cannot allocate memory for pattern ");

	    if (name != NULL)
		IfErr(name = change_array_size(name, Asize + Blksize, char *))
		    Erreturn("cannot allocate memory for pattern ");

	    Asize += Blksize;	/* array size is now Blksize bigger */
	    //cout<<"\n array size is now Blksize bigger Asize="<<Asize<<endl;
	}

	if (first_pattern != NULL) {	/* copy data from first line */
	    for (i = 0; i < lpattern; i++)
		pattern[npattern][i] = first_pattern[i];
	    free(first_pattern);
	    first_pattern = NULL;

	    if (first_name != NULL) {
		name[npattern] = first_name;
		first_name = NULL;
	    }
	}
	else {			/* read data from file */
		//cout<<"\n read data from file \n";
		
	    for (i = 0; i < lpattern; i++) {
		IfEOF(status = fscanf(Pfile, "%s", buffer)) {
		    if (i == 0)
			break;
		    Erreturn1("cannot read pattern # %d", npattern);
		}
		
		//cout<<"\n buffer="<<buffer<<endl;
#ifndef NO_DONTCARES
		if (strcmp(buffer, DONTCARE) == 0)
		    pattern[npattern][i] = DC_VAL;
		else
#endif
		{
		    double f;
		    IfErr (status = sscanf(buffer, "%lf", &f))
			Erreturn1("cannot read pattern # %d", npattern);
		    pattern[npattern][i] = f;
		    //cout<<" f="<<f<<"\t i="<<i<<endl;
		   
		
	    
	    
		}
	    }
	    IfEOF(status) break;

	    if (name != NULL) {
		char    c = skip_blanks(Pfile);

		if (c == '\"') {
		    getc(Pfile);
		    fgets(buffer, sizeof(buffer), Pfile);
		    buffer[strlen(buffer) - 1] = '\0';
		}
		else {
		    IfEOF(status = fscanf(Pfile, "%s", buffer))
			Erreturn("not enough names");
		    skip_blanks(Pfile);
		}

		IfErr(name[npattern] = new_string(buffer))
		    Erreturn("not enough core");
	    }
	    IfEOF(status) break;
	}

	if (!sflag && scales != NULL)	/* scale data if requested */
	    for (i = 0; i < lpattern; i++)
		pattern[npattern][i] *= scales[i];

    }
    cout<<"after read data";
    int ip = 0;
    
    for (row_iterator it_a = (*a).begin();npattern<size && it_a != (*a).end(); npattern++,++it_a, ++ip) {
    
    for (i = 0; i < lpattern; i++) {
    	
    	 
		
		*it_a =pattern[npattern][i];
		
		
    }
    
    }
    
    //doc.open("/media/livingroom/stxxl/stxxl-1.4.1/local/myfile_finalS.txt");
     //i = 0;
		
		/*for (row_iterator it_a = (*a).begin(); it_a != (*a).end(); ++it_a, ++i)
			{
				
				doc>>idf; 
				*it_a =idf;
				cout<<idf<<" ";
				
			}*/
			
			
		i=0;
		//k-nn matrix initialization
		for(row_iterator2 it_b=(*b).begin();it_b!=(*b).end();++it_b,++i)
			{
				if(i%(k+1)==0)
					*it_b=i/(k+1);
				else
					*it_b=99999999;
			}
		
	    
	    cout<<"\n after k-nn matrix initialization";
		int j,p,d;
		double tempV1,tempV2;
		double vect1[dimension],vect2[dimension];
		row_iterator it_c=(*a).begin();
		row_iterator2 it_d=(*b).begin();
		
		double *distanceVect=(double*)malloc(k*sizeof(double));
		double *neighbors=(double*)malloc(k*sizeof(double));
		for(i=0;i < no_of_vector;i++)
			{	
				++it_d;
				for(int ppp=0;ppp<k;ppp++)
				{
					distanceVect[ppp]=999999;
					neighbors[ppp]=99999999;
				}
				for(j=0;j<dimension;j++)
					{
						vect1[j]=*it_c;
						++it_c;
					}
				row_iterator it_a=(*a).begin();
					
				for(p=0;p<no_of_vector;p++)
					{
						
						for(j=0;j<dimension;j++)
							{
								vect2[j]=*it_a;
								++it_a;
							}
						/*double dist=distance2(vect1,vect2,dimension);
					
						if(distanceVect[k-1]>dist)
							{	
								distanceVect[k-1]=dist;
								neighbors[k-1]=p;
								d=k-1;
								while(d>0 && distanceVect[d-1]>distanceVect[d])
									{
										tempV2=distanceVect[d-1];
										tempV1=neighbors[d-1];
										distanceVect[d-1]=distanceVect[d];
										distanceVect[d]=tempV2;
										neighbors[d-1]=neighbors[d];
										neighbors[d]=tempV1;
										d--;	
									}	
					
							}*/
			
						}
						
				/*for(p=0;p<k;p++)
					{
						*it_d=neighbors[p];
						++it_d;
					}
				*/
		
			}
			
	
    /* if there is any error, these pointers below aren't changed */

    /* free any space already allocated for patterns */
    if (*patternP != NULL)
	free_2d_array((char **)*patternP, *npatternP);
    if (*nameP != NULL)
	free(*nameP);

    *npatternP = npattern;	/* # of patterns read in from file */
    *lpatternP = lpattern;	/* # of elements in pattern */
    *patternP = pattern;	/* input array */
    *nameP = name;		/* name array */

    return MY_OK;		/* patterns were read in without error */
}

int nstrings(FILE *fp, FLOAT **patternP, char  **nameP)
{		
	cout<<"\nIn nstring method\n";
    register int i;
    int     c;
    FLOAT  *pattern;
    char   *name = NULL;

    IfErr (pattern = new_array_of(1, FLOAT))
	Erreturn("not enough core");

    for (i = 0;; i++) {
	double f;
		
	if ((c = skip_blanks(fp)) == '\n')
	    break;

	/* read field (number of name) */
	if (c == '\"') {
	    getc(fp);		/* discard quote */
	    fgets(buffer, sizeof(buffer), fp);
	    buffer[strlen(buffer) - 1] = '\0';
	}
	else if (fscanf(fp, "%s", buffer) != 1)
	    break;

#ifndef NO_DONTCARES
	if (c != '\"' && strcmp(buffer, DONTCARE) == 0)
	    pattern[i] = DC_VAL;
	else
#endif
	if (c == '\"' || sscanf(buffer, "%lf", &f) != 1) {
	    IfErr (name = new_string(buffer))
		Erreturn("not enough core");

	    if (c != '\"')
		skip_blanks(fp);
	    break;
	}
	pattern[i] = f;

	IfErr(pattern = change_array_size(pattern, (i + 2), FLOAT))
	    Erreturn("not enough core");
    }
	
	
	
    if (i == 0)
	Erreturn("empty pattern");

    *patternP = pattern;
    *nameP = name;
    return i;
}

void merge(BiTree *item, int pair1, int pair2, int lpat, float min_dist)
{
	//cout<<"\n In Merge Method";
	register int i, j, k;
	BiTree *newitem;
	int newsize = item[pair1].size+item[pair2].size;
	FLOAT r1 = (float)item[pair1].size/newsize;
	FLOAT r2 = (float)item[pair2].size/newsize;
	FLOAT maxDist, *repPoint;

	IfErr (newitem = new_tree(&item[pair1])) { /* copy */
	    fprintf(stderr, "%s: not enough core\n", program);
	    exit(1);
	}

	/*
	 * replace left child node with new tree node
	 * link right child node into parent
	 */
	item[pair1].l_tree = newitem;
	item[pair1].r_tree = &item[pair2];
	item[pair1].leaf = LEAF;	/* ith item cannot be a leaf it has
					 * subtrees */

	/* set the size and clear up some memory */
	item[pair1].size = newsize;
	item[pair1].distance = min_dist;

	newitem->root = FALSE;
	item[pair2].root = FALSE;	/* jth item is no longer a root.its a
					 * subtree */

	for (i=0; (i<rsize && i<item[pair2].size); ++i)
		free(item[pair2].rep[i]);
	free(item[pair2].rep);

	/* find mean of two clusters as weighted average*/
	for (i=0; i<lpat; ++i)
		item[pair1].mean[i] = r1*item[pair1].mean[i] + r2*item[pair2].mean[i];
	free(item[pair2].mean);

/*
printf("merge mean : %f %f\n", item[pair1].mean[0], item[pair1].mean[1]);
*/
	/* find new representative points */
	for (i=0; (i < rsize && i < newsize); ++i) {
		maxDist = 0;
		repInTree(&item[pair1], i, item[pair1].mean, item[pair1].rep, 
			  lpat, &maxDist, &repPoint);
		memcpy(item[pair1].rep[i], repPoint, lpat*sizeof(FLOAT));
/*
printf("rep %d before : %f %f\n", i, item[pair1].rep[i][0], item[pair1].rep[i][1]);
*/
	}

	/* shrink based on alpha */
	for (j=0; j<i; ++j)
		for (k=0; k<lpat; ++k)
			item[pair1].rep[j][k] += (alpha*(item[pair1].mean[k]-
					                 item[pair1].rep[j][k]));
}

void repInTree(BiTree *item, int rcount,FLOAT *mean, FLOAT **rep, int lpat, FLOAT *maxDist, FLOAT **repPoint)
{
	//cout<<"\n In repInTree Method";
    register int i;
    FLOAT minDist, dist;

    if (item->leaf != LEAF) {
	if (rcount == 0)
		minDist = distance(item->pat, mean, lpat);
	else {
		minDist = MAXFLOAT;
		for (i=0; i<rcount; ++i) 
			if ((dist = distance(item->pat, rep[i], lpat)) < minDist)
				minDist = dist;
	}

	if (minDist >= *maxDist) {
		*maxDist = minDist;
		*repPoint = item->pat;
	}
    }
    else {
		repInTree(item->l_tree, rcount, mean, rep, lpat, maxDist, 
			  repPoint);
		repInTree(item->r_tree, rcount, mean, rep, lpat, maxDist, 
			  repPoint);
    }
}

void print_rep_points( BiTree *item, int cno, int lpat)
{
	cout<<"\n print_rep_points Method\n";
    register int i, j;

    printf("**** Cluster %d from CURE Rep Points ****\n", cno);
    for (i=0; i<item->size && i<rsize; ++i) {
	for (j=0; j<lpat; ++j)
		printf("%f ", item->rep[i][j]);
	printf("\n");
    }
}

/* allocates space for a N1 * N2 * N3 array of specified element size.
	returns the pointer (char*) - type must be casted by the caller */

char ***calloc_3d(unsigned N1, unsigned N2, unsigned N3, unsigned size)
{
	cout<<"\n In calloc MEthod";
    register int i, j;
    char ***ptr = (char ***) calloc(N1, sizeof(char *));
    if (ptr == NULL)
	return NULL;

    for (i = 0; i < N1; i++) {
	if (NULL == (ptr[i] = (char **) calloc(N2, sizeof(char *))))
	    return NULL;
	for (j = 0; j < N2; j++)
	    if (NULL == (ptr[i][j] = (char *)calloc(N3, size)))
		return NULL;
    }
    return ptr;
}

/* change size of 2d array pointed by "array" from N1xN2 to M1xM2. */
/* contents are undisturbed up to the lesser of N1/M1 and N2/M2 */

char  **realloc_2d(char **array, unsigned N1, unsigned N2,unsigned M1,unsigned M2,unsigned size)
{
	cout<<"\n In realloc Method";
    char  **ptr;
    register int i, j;

    if (NULL == (ptr = (char **) calloc_2d(M1, M2, size)))
	return NULL;
    N2 *= size;
    M2 *= size;
    for (i = 0; i < N1 && i < M1; i++)
	for (j = 0; j < N2 && j < M2; j++)
	    ptr[i][j] = array[i][j];
    free_2d_array(array, N1);

    return ptr;
}

/* frees a 2d array.  N is range of 1st index */

void free_2d_array(char  **array,int N)
{
	cout<<"\n In free_2d_array MEthod ";
    register int i;
    for (i = 0; i < N; i++)
	free(array[i]);
    free((char *) array);
}

/* frees a 3d array.  N1/N2 are ranges of 1st/2nd index */

void free_3d_array(char ***array, int N1, int N2)
{
	cout<<"\n In free_3d_array MEthod ";
    register int i;
    for (i = 0; i < N1; i++)
	free_2d_array(array[i], N2);
    free((char *) array);
}

/* allocates space for a N1 * N2 array of specified element size.
	returns the pointer (char*) - type must be casted by the caller */

char  **calloc_2d(unsigned N1, unsigned N2, unsigned size)
{
	cout<<"\n In calloc_2d MEthod ";
    register int i;
    char  **ptr = (char **) calloc(N1, sizeof(char *));
    if (ptr == NULL) {
	return NULL;
    }
    for (i = 0; i < N1; i++)
	if (NULL == (ptr[i] = (char *)calloc(N2, size)))
	    return (NULL);
    return ptr;
}

char  *new_string(char *string)
{
	cout<<"\n In new_string Method";
    char   *buf;
    IfErr(buf = new_array_of(strlen(string) + 1, char))
	return MY_ERR;
    strcpy(buf, string);
    return buf;
}
