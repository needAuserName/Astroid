/*    
 *    splice_graph.h		
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

#ifndef SPLICE_GRAPH
#define SPLICE_GRAPH

#include "Configuration.h"



class GTvertex;
class fragment;

enum fragment_type {frag_exon, frag_intron, frag_retained_intron, frag_junction, frag_insertion, frag_deletion, frag_SNP};
enum trans_direction {undetermined, sense, antisense};

/************************************************************************/
/* fragment information                                                 */
/************************************************************************/
class alter_junction
{
public:
	fragment *juncInfo;
	//	int category; //0 for the merged target junction, -1 for 5', 1 for 3', 2 for exon
	//	double proportion;
	alter_junction *next;

	alter_junction();
};

class fragment
{
public:
	char name[30];
	unsigned long ID; //unique ID throughout the pipeline
	char chromosome_start[100];
	char chromosome_end[100];
	trans_direction transDirection;
	long start;
	long end;
	fragment_type type; //fragment type
	int altersite; //if the fragment is an exon, then this field indicates whether this fragment is an alternative splice site: 0 for exon, 1 for donor site, 2 for acceptor site

	double support[SUPPORT_VECTOR_SIZE + 1]; //support from the fragment file

	alter_junction* alter;
	long alterFragCnt; //fragment count of alternative splice sites, note it's the number of fragments
	//	alter_junction* fivePalterTail; //tail pointer of 5' alternative junction
	long *coverage;
	long start_real; //real start position in case the position is changed due to alternative splice sites
	long end_real; //real end position

	fragment();
	fragment* clone(); //make a clone, used for exon, so alter_junction not included
	~fragment();
};

class spliceSite
{
public:
	char chromosome[100];
	long position;
	bool directionOut; //true for out-going junction, false for in-coming junction

	spliceSite();
};


/************************************************************************/
/* GenomeTree                                                           */
/************************************************************************/
class rangeJunction
{
	//fragment within a range
public:
	fragment *junc;
	rangeJunction *next;

	rangeJunction();
};

class RangeJunctionList
{
	//fragment list within a range
public:
	long rangeLow;
	long rangeHigh;
	unsigned long cnt_rangeJunction;
	trans_direction transDirection;
	rangeJunction *list;
	rangeJunction *listtail;
	RangeJunctionList *nextList;

	RangeJunctionList();
	RangeJunctionList* clone();
	void insert_feature(rangeJunction *new_rangeJunction);
	void count_featurelist();
	~RangeJunctionList();
};


class GTedge
{
	//edge of Genome Tree
public:
	GTvertex *linkedVertex;
	GTedge *next;

	GTedge();
};

class alternative_path
{
public:
	long path_start; //site that starts diverging
	long path_end; //site that ends diverging
	long whole_path_start; //the whole path including the starting exon
	long whole_path_end; //the whole path including the ending exon
	trans_direction transDirection;
	double support[SUPPORT_VECTOR_SIZE]; //support at this vertex
	double proportion[SUPPORT_VECTOR_SIZE]; //proportion at this vertex
	long junctionNum;
	long exonNum;
	GTvertex *pathVertex;
	alternative_path *next;

	alternative_path();
};

class GTvertex
{
	//vertex of Genome Tree
public:
	long ID;
	long level;
	long rangeLow;
	long rangeHigh;
	GTedge *child;
	int childType; //0 for no child (i.e. leaf node), 1 for independent region, 2 for independent path, 3 for dependent path
	long childNum;
	RangeJunctionList *junctionInRange; //junctions within the vertex's range
	long junctionNum;
	long exonNum;
	long IndepASMNum;
	long exonicLength;
	GTvertex *prevSibling; //the sibling on the left (5') side
	GTvertex *nextSibling; //the sibling on the right (3') side

	GTedge *alterSpliceSite; //alternative splice sites involved in this vertex

	//for difference analysis
	double support[SUPPORT_VECTOR_SIZE]; //support at this vertex
	double proportion[SUPPORT_VECTOR_SIZE]; //proportion at this vertex (for genes, this array stores mean 25%-75% coverage)
	double MSE_estimation[SUPPORT_VECTOR_SIZE]; //mean squared error of estimated transcript abundance
	double min_path_support[SUPPORT_VECTOR_SIZE];
	double obs_support[SUPPORT_VECTOR_SIZE]; //observed support directly calculated from raw support (with no estimation)

	double total_inflow;
	double total_outflow;

	bool estimated; //if true, the vertex has been estimated and will be treated as a whole
	fragment *representative;
	long estimate_exonNum; //number of blocks under estimation, equivalent to number of exons previously

	//statistics
	//these tests are conducted about the children, and are summarized at this vertex
	double expression_F_score;
	double expression_pvalue;
	double expression_winthinGroupRank;
	double JSDsqrt_mean;
	double JSDsqrt_combination;
	double JSD_Pvalue_mean;
	double JSD_Pvalue_combination;
	double JSD_withinGroupRank_mean;
	double JSD_withinGroupRank_combination;
	bool JSD_reliable;
	double pearsonCorr_support;
	double pearsonCorr_proportion;
	double divergence_stat;
	double divergence_stat_rank;

	int ASMcategory;
	double ASMsupport_group1;
	double ASMsupport_group2;

	alternative_path* major_alter_paths;
	long major_alter_paths_num;

	GTvertex();
};

class GenomeTree
{
public:
	GTvertex* root;

	GenomeTree();
	~GenomeTree();
};

//global variables
const long MAX_CHR_LENGTH = 1000000000; //maximal chromosome length
const long MAX_JUNCTION_NUM = 1000000;




/************************************************************************/
/* JUNCTION GRAPH                                                       */
/************************************************************************/


class JuncGraphVertex;
enum JuncGraphVertexType {virStart, virEnd, normal};

class JuncGraphEdge
{
	//edge of fragment graph
public:
	JuncGraphVertex *linkedVertex;
	JuncGraphEdge *next;

	JuncGraphEdge();
};

class JuncGraphVertex
{
	//vertex of fragment graph
public:
	rangeJunction* corresJunc;
	JuncGraphEdge* edges;
	JuncGraphVertex* next;
	bool traversed;
	bool hasInEdge; //if hasInEdge is false, consider it as a start vertex when enumerating transcripts
	bool hasOutEdge;

	JuncGraphVertexType vertexType;

	JuncGraphVertex();
	~JuncGraphVertex();
};

class JuncGraphPath
{
	//path of fragment graph
public:
	RangeJunctionList *pathJuncList;
	JuncGraphVertex *arrivedVertex;

	JuncGraphPath();
	~JuncGraphPath();
};

class JuncGraph
{
	//fragment graph
public:
	JuncGraphVertex *vertices;

	JuncGraph();
	~JuncGraph();
};


const double coverageThreshold_exon = 5; //throw-away := <= coverageThreshold * SUPPORT_VECTOR_SIZE
const double coverageThreshold_intron = 10; //throw-away := <= coverageThreshold * SUPPORT_VECTOR_SIZE
const long MIN_EXON_LENGTH = 1;
const long MIN_ALTER_SPLICE_SITE_LENGTH = 1; //minimum exonic length for alternative splice sites. need a maximum??
const long MAX_NOCOVERAGE_LENGTH = 1; //maximal length for an exon to have no coverage
const double MAX_NOCOVERAGE_THRESH_LOW_COV_REGION = 5;
const long MAX_NOCOVERAGE_LENGTH_LOW_COV_REGION = 50;
//const long MAX_NOCOVERAGE_LENGTH = 100; //for fake gene examples
const double MAX_JUNSUPPORT = 5; // for junction filtering
const double MEAN_JUNSUPPORT = 5; // for junction filtering

//filter exon in gene if its coverage is less than 0.05 of mean gene coverage in more than 0.95 samples
const double coverageThreshold_GeneExon = 0; //percentage of mean gene coverage

const double MIN_COVERAGE_FILTER_FRAGMENT = 0; //minimum coverage for fragments when filtering standalone fragments.

const int MAXJUNCNUMINDEPPATH = 50; //maximum number of junctions when separating dependent paths. if #junctions > this number, do not enumerate paths
const int MAX_NUM_ENUMERATED_PATH = 50; //maximum number of dependent paths when enumerating, this is set in the purpose to replace the above one. try to enumerate first, stop and keep NOTHING if exceeding this limit

const int MINJUNCTIONLENGTH = 5; // minimal length of the junction that is considered valid (For real data: original setting = 10)

const long MAX_ITERATION = 1000; //for calculating p value

const int COVERAGE_CHANGE_WINDOW = 10;
const double COVERAGE_CHANGE_THRESH = 5.0;

#endif

