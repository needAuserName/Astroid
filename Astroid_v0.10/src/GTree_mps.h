/*    
 *    GTree_mps.h		
 *
 *    Copyright (C) 2014 University of Kentucky and
 *                       Yan Huang
 *
 *    Authors: Yan Huang
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef GTREE_MPS
#define GTREE_MPS

#include "Configuration.h"
#include "splice_graph.h"
/************************************************************************/
/* HEADER FILE                                                          */
/************************************************************************/

bool DOTRIMMING = true;
bool DORETAIN = true;


bool CLEAN_SPLICE_DIRECTION = false; //need correcting transcription strand for splice junctions (too splice junctions have same location but different strands)


long junctionNum = 0;

GenomeTree *gTree;

const long MAX_STACK_SIZE = MAX_JUNCTION_NUM;
GTvertex* stack_vertex[MAX_STACK_SIZE];
long stackTop; 


void constructGTree(RangeJunctionList* origList);



JuncGraph *junctionGraph;

//for sorting junctions
rangeJunction* sortList_Junction[MAX_JUNCTION_NUM * 3]; //fragment list to be sorted
long sortKey_Junction[MAX_JUNCTION_NUM * 3];
long sortList_Junction_Num = 0;

long mergeSort_Larray[MAX_JUNCTION_NUM * 2 + 1];
long mergeSort_Rarray[MAX_JUNCTION_NUM * 2 + 1];
rangeJunction* mergeSort_LorderedList[MAX_JUNCTION_NUM * 2 + 1];
rangeJunction* mergeSort_RorderedList[MAX_JUNCTION_NUM * 2 + 1];

//for sorting splice sites
spliceSite* sortList_spliceSite[2 * MAX_JUNCTION_NUM]; //fragment list to be sorted
long sortKey_spliceSite[2 * MAX_JUNCTION_NUM];
long sortList_spliceSite_Num;

spliceSite* mergeSort_LorderedList_spliceSite[MAX_JUNCTION_NUM + 1];
spliceSite* mergeSort_RorderedList_spliceSite[MAX_JUNCTION_NUM + 1];

long quicksortStack[10000];

//junction graph path queue for computing all paths
const int MAX_PATH_NUM = 1000000; 
JuncGraphPath* juncPathQueue[MAX_PATH_NUM + 1]; 
long juncPathQueueTail;
long juncPathQueueHead;
long juncPathQueueMaxLength;




void mergeSort_JunctionSort(long sortList_size);
RangeJunctionList* separateIndepRegion(RangeJunctionList* origList, bool &separable);
RangeJunctionList* separateIndepPath(RangeJunctionList* origList, bool &separable);
double abundance_estimation(GTvertex *rootVertex, int index_tissue, double &MSE);
bool countGTree(GTvertex *rootVertex);
double JSD_significance_level(double JSDvalue, long dimension, double N);
void calculate_ASM_group_meanExpression(GTvertex *rootVertex);
bool alterSpliceReliability(GTvertex *curVertex);
void output_ASMpath_gtf(GTvertex *targetVertex, RangeJunctionList *pathList);
void filterFalseFragments(bool	);

//main output
ofstream gtreefile;
ofstream gTreePlotTree;
ofstream leafsizefile;
ofstream diffRegionFile;
ofstream unresolvableFile;
ofstream annotationASM_file;

//other output
ofstream targetfile;
ofstream abnormalfile;
ofstream exon_skipping_file;
ofstream mutual_exclusive_file;
ofstream retainedIntron;
ofstream alterStartEnd;
long example_exon_skipping_id = 1;
ofstream example_exon_skipping_file;
ofstream test_D_statistic_dist_file;



///////////////////////////////////////////////////////////
//General parameters


const int unknown = 0;
const int exon_skipping = 1;
const int mutual_exclusive = 2;
const int intron_retention = 3;
const int alter_splice_site = 4;
const int alter_start = 5;
const int alter_end = 6;

const int diff_expression = 7;

long exon_skipping_cnt = 0;
long mutual_exclusive_cnt = 0;
long intron_retention_cnt = 0;

long CHROMOSOME_START = 0; // 43737953; // 763064;
long CHROMOSOME_END   = MAX_CHR_LENGTH; // 43754221; // 789740;
// const long CHROMOSOME_START = 763064;
// const long CHROMOSOME_END   = 789740;

unsigned long fragmentID_Cnt = 0;

const long targetRegion_start = 16500118;
const long targetRegion_end = 16516729;



//class for filtering
class deletedSites
{
public:
	long sites;
	deletedSites *next;

	deletedSites();
};

//distribution pool for the sampling of JSD space
int JSDpool_dimension_list[1000];
int JSDpool_dimension_list_cnt = 0;


long DatasetReadCount_total[SUPPORT_VECTOR_SIZE];
const long DatasetReadCount_normalization = 30000000;
double DatasetNormalizationRatio[SUPPORT_VECTOR_SIZE];


//For collecting d statistics
const long MAX_VERVEX_NUM_FOR_STATISTICS = 1000000;
long vertexListForStatistics_cnt = 0;
GTvertex* vertexListForStatistics[MAX_VERVEX_NUM_FOR_STATISTICS];


ofstream statistics_expression_file;
ofstream statistics_asm_file;
ofstream ASMpath_gtf_file;
ofstream NonASM_gtf_file;

ofstream splice_all_file;
ofstream splice_filtered_file;

ofstream STAT_ASM_composition; //output composition of ASMs

string inputPath;
string inputPath_SAM;
string resultPath;
string GTFFile;
string GapFile;
string resultNamePrefix;

long GENEcount = 0;
long ASMcount = 0;



/************************************************************************/
/* For Assembly                                                         */
/************************************************************************/

class graphpath
{
public:
	trans_direction direction;
	long pathlength;
	long gaplength;
	JuncGraphVertex *vertex_start;
	JuncGraphVertex *vertex_end;
	JuncGraphEdge *edgelist;
	JuncGraphEdge *edgelist_tail;
	int pathindex;

	graphpath *next;

	graphpath();
	graphpath* clone(); 
	void add_edge(JuncGraphEdge *new_edge); //add new edge to edge list tail 
	~graphpath();
};

class pathCluster
{
public:
	rangeJunction *fragstart;
	rangeJunction *fragend;
	graphpath *paths;

	pathCluster();
	~pathCluster();

};

class pathlength_record
{
public:
	long lowerbound;
	long  higherbound;
	pathlength_record *next;

	pathlength_record();
};


class connectedRead
{
public:
	long readlowIndex;
	long readhighIndex;
	int pathselectedIndex;
	pathlength_record *lengthlist;
	long minLength_inbetween;

	connectedRead();
	~connectedRead();
};


class graphpathCluster
{
public:
	long readLowIndex;
	long readHighIndex;
	graphpath *paths;

	graphpathCluster();
	graphpathCluster* clone();
	~graphpathCluster();
};

class ReadVertex
{
public:
	string readname;
	long start;
	long end; // after clustering, keep min end / max start
	long start_min;
	long end_max;
	long readLength;
	long Id; // index in the flow network
	int includedReadNum; // for clustering
	vector<int> junctionsIndex; // a number illustrate the junction set contained by each read
	trans_direction transDirection;
	bool spliced;
	rangeJunction *list;
	rangeJunction *listTail;
	ReadVertex *next;

	ReadVertex* clone();
	ReadVertex();
	~ReadVertex();
};


class ReadPath
{
public:
	ReadVertex *pathVertices;
	ReadVertex *firstVertex;
	ReadVertex *lastVertex;
	trans_direction transDirection;
	bool spliced;
	ReadPath *next;

	ReadPath* clone();
	ReadPath();
	~ReadPath();
};

// For aggregate transcript fragments
class ReadIndex
{
public:
	int readIndex;
	ReadIndex *next;

	ReadIndex(){readIndex = MAX_NUMBER; next = NULL;}
};


class ReadIndicesPath
{
public:
	ReadIndex *readIndices;
	ReadIndex *lastreadIndex;
	int copyNum;
	trans_direction transDirection;
	ReadIndicesPath *next;

	ReadIndicesPath();
	~ReadIndicesPath();
};


class Transcript
{
public:
	string transName;
	long copyNum;
	ReadPath *readpaths;
	RangeJunctionList *transcript;
	Transcript *nextList;

	Transcript *clone();
	Transcript();
	~Transcript();
};

class geneExonBoundary
{
public:
	long position;
	geneExonBoundary *next;

	geneExonBoundary();
};


/************************************************************************/
/* Date: 09/19/2013 */
/* Structure saves ASM and deterministic read information */
/************************************************************************/
class TransInfo
{
public:
	vector<int> FragIndexCovered;
	long transLength;
	TransInfo *nextTransinfo;

	TransInfo();
	~TransInfo();
};

class ReadInfo
{
public:
	vector<int> FragIndexCovered;
	ReadInfo *nextreadInfo;

	ReadInfo();
	~ReadInfo();
};


class AdjASMsInfo
{
	// save two adjacent ASMs
public:
	long firstASM_start;
	long secondASM_end;
	long Hubstart;
	long Hubend;
	long Hublength; // the exonic length of the hub
	TransInfo *ASM_Transinfo;
	ReadInfo *deterministicReadinfo;

	AdjASMsInfo *nextinfo;
	AdjASMsInfo();
	~AdjASMsInfo();
};


JuncGraph *spliceGraph;
vector<graphpathCluster *> ReadpathCluster; // for any two gapped reads, build all possible paths between them
vector<pathCluster *> FragpathCluster;
vector<connectedRead *> ConnectedReadCluster;

// Input read files
vector<ReadVertex *> ReadList;
vector<ReadVertex *> ReadList_cluster;
long ReadNum;
long gapConstraint;
vector<long> gapList, gapList_partial;


ifstream inputfile_read;
bool input_read = true;
ofstream outputfile_gtf, outputfile_ASMgtf, outputfile_readpath, outputfile_stat, outputfile_gaps;

// Gene locus and adjacent ASMs locus
vector<long> ReadVertexIndex; 


const long MAX_ITER = 3;
const double MAX = 1000000000; 
const double epsilon = 1e-3;

// heuristic
const long MAX_READNUM_ALLOWED = 300;
const long MAX_DOABLE_READNUM = 500000;

// heuristic for removing readpaths and transcript
// (disabled) 1. Calculate the probability of each read path, remove the non significant 5% of them 
// 2. If the copy number of one isoform is less than 5% of the total isoform copy number, remove the isoform
const double threshold = 0.05;

long GeneStart, GeneEnd;

// likelihood ratio test parameters
long TranscriptNum;
Transcript **TransVec;

// global parameters for matching
long vertexNum;

// For gradient descent method to find the transcript copy number with minimum joint likelihood
const int blocksize = 10;
const int maxgapConstraint = 1000;

// For read fragment clustering
long ReadNumConstraint = 30;
const long MAXJUNTOLERENT = 20;
const long MAXCLUSTERSIZE = 1000;

// For Deterministic read fragment part
vector<double> TransProb_ASM;
double constant_Penalty = 3;
long TransNum; // roughly estimate the transcript number

// For removing transcript copy paths which cover very few of a transcript
double threshold_coverage = 0.1;

#endif