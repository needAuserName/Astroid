/*    
 *    Assembly.cpp		
 *    Assembly
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

#include "GTree_mps.h"
#include "mathfunction.h"
#include "general_functions.h"

void merge_JunctionSort(long p, long q, long r)
{
	long n1, n2, i, j, k;

	n1 = q - p + 1;
	n2 = r - q;

	for (i = 1; i <= n1; i++)
	{
		mergeSort_Larray[i] = sortKey_Junction[p + i - 1];
		mergeSort_LorderedList[i] = sortList_Junction[p + i - 1];
	}
	for (j = 1; j <= n2; j++)
	{
		mergeSort_Rarray[j] = sortKey_Junction[q + j];
		mergeSort_RorderedList[j] = sortList_Junction[q + j];
	}

	mergeSort_Larray[n1 + 1] = MAX_CHR_LENGTH * 2;
	mergeSort_Rarray[n2 + 1] = MAX_CHR_LENGTH * 2;

	i = 1;
	j = 1;

	for (k = p; k <= r; k++)
	{
		if (mergeSort_Larray[i] <= mergeSort_Rarray[j])
		{
			sortKey_Junction[k] = mergeSort_Larray[i];
			sortList_Junction[k] = mergeSort_LorderedList[i];

			i++;
		} 
		else
		{
			sortKey_Junction[k] = mergeSort_Rarray[j];
			sortList_Junction[k] = mergeSort_RorderedList[j];

			j++;
		}
	}

	return;
}


void mergeSort_JunctionSort(long sortList_size)
{
	//non-recursive merge sort for sorting junctions
	long m, n, i, r;
	m = 1;
	n = sortList_size;

	while (m <= n)
	{
		i = 1;
		while (i <= n - m)
		{
			r = (i + 2 * m - 1) < n ? (i + 2 * m - 1) : n;
			merge_JunctionSort(i, i + m - 1, r);
			i = i + 2 * m;
		}

		m = m * 2;
	}

	return;
}

long partition(long p,long r, double *sortArray)
{
	long i, j;
	double x, tmp;

	//randomized partition
	i = p + (double)rand()/ (RAND_MAX) * (r - p);

	if (sortArray[r] != sortArray[i])
	{
		tmp = sortArray[r];
		sortArray[r] = sortArray[i];
		sortArray[i] = tmp;
	}

	x = sortArray[r];
	i = p - 1;

	for(j = p; j < r; j++)
	{
		if (sortArray[j] <= x)
		{
			i++;

			if (sortArray[i] != sortArray[j])
			{
				tmp = sortArray[j];
				sortArray[j] = sortArray[i];
				sortArray[i] = tmp;
			}
		}
	}

	if (sortArray[r] != sortArray[i+1])
	{
		tmp = sortArray[r];
		sortArray[r] = sortArray[i+1];
		sortArray[i+1]=tmp;
	}

	return i+1;
}

void quicksort(double *sortArray, long length)
{
	long top = 0, p, r, q;

	quicksortStack[top++] = 1;
	quicksortStack[top++] = length;

	while (top != 0)
	{
		r = quicksortStack[--top];
		p = quicksortStack[--top];

		if(p>=r)
			continue;

		q = partition(p, r, sortArray);

		quicksortStack[top++] = p;
		quicksortStack[top++] = q - 1;

		quicksortStack[top++] = q + 1;
		quicksortStack[top++] = r;
	}

	return;
}

//RangeJunctionList Queue
//for separating dependent paths
void juncPathQueue_initialization()
{
	juncPathQueueHead = 0;
	juncPathQueueTail = 0;
	juncPathQueueMaxLength = MAX_PATH_NUM;

	return;
}

bool juncPathQueue_enqueue(JuncGraphPath *x)
{
	juncPathQueue[juncPathQueueTail] = x;

	if (juncPathQueueTail == juncPathQueueMaxLength)
	{
		juncPathQueueTail = 0;
	}
	else 
	{
		juncPathQueueTail++;
	}

	return true;
}

JuncGraphPath* juncPathQueue_dequeue()
{
	JuncGraphPath *x;

	if (juncPathQueueHead < juncPathQueueTail)
	{
		x = juncPathQueue[juncPathQueueHead];

		if (juncPathQueueHead == juncPathQueueMaxLength)
		{
			juncPathQueueHead = 0;
		}
		else 
		{
			juncPathQueueHead++;
		}

		return x;
	} 
	else
	{
		return NULL;
	}
}

/************************************************************************/
/* Vertex Stack for DFS                                                 */
/************************************************************************/
void stack_initial()
{
	stackTop = 0;

	return;
}

bool stack_empty()
{
	if (stackTop == 0)
	{
		return true;
	} 
	else
	{
		return false;
	}
}

void stack_push(GTvertex *x)
{
	if (stackTop >= MAX_STACK_SIZE)
	{
		cout << "stack error" << endl;
		exit(1);
	}

	stackTop++;
	stack_vertex[stackTop] = x;

	return;
}

GTvertex* stack_pop()
{
	if (stack_empty() == true)
	{
		return NULL;
	} 
	else
	{
		stackTop--;
		return stack_vertex[stackTop+1];
	}
}


/************************************************************************/
/* JUNCTION GRAPH                                                       */
/************************************************************************/


JuncGraphEdge::JuncGraphEdge()
{
	linkedVertex = NULL;
	next = NULL;
}

JuncGraphVertex::JuncGraphVertex()
{
	corresJunc = NULL;
	edges = NULL;
	next = NULL;
	traversed = false;
	hasInEdge = false;
	hasOutEdge = false;
	vertexType = normal;
}

JuncGraphVertex::~JuncGraphVertex()
{
	if (corresJunc != NULL)
	{
		delete corresJunc;
	}

	JuncGraphEdge *delEdge;

	while (edges != NULL)
	{
		delEdge = edges;
		edges = delEdge->next;
		delete delEdge;
	}
}

JuncGraphPath::JuncGraphPath()
{
	pathJuncList = NULL;
	arrivedVertex = NULL;
}

JuncGraphPath::~JuncGraphPath()
{
	arrivedVertex = NULL;
	delete pathJuncList;
}

JuncGraph::JuncGraph()
{
	vertices = NULL;
}

JuncGraph::~JuncGraph()
{
	JuncGraphVertex *delVertex;

	while (vertices != NULL)
	{
		delVertex = vertices;
		vertices = delVertex->next;
		delete delVertex;
	}
}



/************************************************************************/
/* GENOME TREE	                                                        */
/************************************************************************/


alter_junction::alter_junction()
{
	juncInfo = NULL;
	//	category = 0;
	//	proportion = 0.0;
	next = NULL;
}


fragment::fragment()
{
	for (int i = 0; i < 30; i++)
	{
		name[i] = '\0';
	}
	ID = 0;
	transDirection = undetermined;
	for (int i = 0; i < 100; i++)
	{
		chromosome_start[i] = '\0';
		chromosome_end[i] = '\0';
	}

	start = 0;
	end = 0;

	type = frag_junction;
	altersite = 0;

	for (int i = 0; i <= SUPPORT_VECTOR_SIZE; i++)
	{
		support[i] = 0.0;
	}

	alter = NULL;
	//	fivePalterTail = NULL;
	coverage = NULL;

	start_real = 0;
	end_real = 0;
}

fragment* fragment::clone()
{
	fragment *newFrag = new fragment;

	strcpy(newFrag->name, name);
	newFrag->ID = ID;
	strcpy(newFrag->chromosome_start, chromosome_start);
	strcpy(newFrag->chromosome_end, chromosome_end);
	newFrag->start = start;
	newFrag->end = end;
	newFrag->type = type;
	newFrag->transDirection = transDirection;

	for (int iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
	{
		newFrag->support[iLoop] = support[iLoop];
	}

	newFrag->alter = NULL;
	//	newFrag->fivePalterTail = NULL;
	newFrag->coverage = NULL;
	newFrag->start_real = start_real;
	newFrag->end_real = end_real;

	return newFrag;
}

fragment::~fragment()
{
	alter_junction *delJunc;
	delJunc = alter;
	while (delJunc != NULL)
	{
		alter = delJunc->next;
		delete delJunc;
		delJunc = alter;
	}

	if (coverage != NULL)
	{
		delete [] coverage;
	}
}

spliceSite::spliceSite()
{
	position = 0;
	directionOut = true;
}

rangeJunction::rangeJunction()
{
	junc = NULL;
	next = NULL;
}

RangeJunctionList::RangeJunctionList()
{
	rangeLow = MAX;
	rangeHigh = 0;
	cnt_rangeJunction = 0;
	transDirection = undetermined;
	list = NULL;
	listtail = NULL;
	nextList = NULL;
}

RangeJunctionList* RangeJunctionList::clone()
{
	//clone a same RangeJunctionList
	RangeJunctionList *resultList;
	resultList = new RangeJunctionList;

	resultList->rangeLow = rangeLow;
	resultList->rangeHigh = rangeHigh;
	resultList->transDirection = transDirection;
	resultList->list = NULL;
	resultList->listtail = NULL;
	resultList->nextList = NULL;

	rangeJunction *curList, *newList;
	curList = list;
	while (curList != NULL)
	{
		newList = new rangeJunction;
		newList->junc = curList->junc;
		newList->next = NULL;

		if (resultList->list == NULL)
		{
			resultList->list = newList;
			resultList->listtail = newList;
		}
		else
		{
			resultList->listtail->next = newList;
			resultList->listtail = newList;
		}

		curList = curList->next;
	}

	return resultList;
}

void RangeJunctionList::insert_feature(rangeJunction *new_rangeJunction)
{
	//insert the new feature to tail of the feature list
	if (listtail == NULL)
	{
		assert(list == NULL);

		list = new_rangeJunction;
		listtail = new_rangeJunction;
	}
	else
	{
		assert(listtail->next == NULL);

		listtail->next = new_rangeJunction;
		listtail = new_rangeJunction;
	}

	++cnt_rangeJunction;

	return;
}

void RangeJunctionList::count_featurelist()
{
	//count number of features in the feature list
	cnt_rangeJunction = 0;
	rangeJunction *cur_feature = list;
	while (cur_feature != NULL)
	{
		++cnt_rangeJunction;
		cur_feature = cur_feature->next;
	}

	return;
}

RangeJunctionList::~RangeJunctionList()
{
	nextList = NULL;

	rangeJunction *delList;

	while (list != NULL)
	{
		delList = list;
		list = delList->next;
		delete delList;
	}
}

GTedge::GTedge()
{
	linkedVertex = NULL;
	next = NULL;
}

alternative_path::alternative_path()
{
	path_start = 0;
	path_end = 0;
	whole_path_start = 0;
	whole_path_end = 0;
	transDirection = undetermined;
	for (int i = 0; i < SUPPORT_VECTOR_SIZE; i++)
	{
		support[i] = 0.0;
		proportion[i] = 1.0;
	}
	junctionNum = 0;
	exonNum = 0;
	pathVertex = NULL;
	next = NULL;
}

GTvertex::GTvertex()
{
	ID = 0;
	level = 0;
	rangeLow = 0;
	rangeHigh = 0;
	child = NULL;
	childType = 0;
	childNum = 0;
	junctionInRange = NULL;
	junctionNum = 0;
	exonNum = 0;
	IndepASMNum = 0;
	exonicLength = 0;
	prevSibling = NULL;
	nextSibling = NULL;

	alterSpliceSite = NULL;

	for (int i = 0; i < SUPPORT_VECTOR_SIZE; i++)
	{
		support[i] = 0.0;
		proportion[i] = 1.0;
		MSE_estimation[i] = 0.0;
		min_path_support[i] = 0.0;
		obs_support[i] = 0.0;
	}

	total_inflow = 0.0;
	total_outflow = 0.0;

	estimated = false;
	representative = NULL;
	estimate_exonNum = 0;

	expression_F_score = 0.0;
	expression_pvalue = 1.0;
	expression_winthinGroupRank = SUPPORT_VECTOR_SIZE;
	JSDsqrt_mean = 0.0;
	JSDsqrt_combination = 0.0;
	JSD_Pvalue_mean = 1.0;
	JSD_Pvalue_combination = 1.0;
	JSD_withinGroupRank_mean = SUPPORT_VECTOR_SIZE;
	JSD_withinGroupRank_combination = SUPPORT_VECTOR_SIZE;
	JSD_reliable = true;
	pearsonCorr_support = 0.0;
	pearsonCorr_proportion = 0.0;
	divergence_stat = 0.0;
	divergence_stat_rank = SUPPORT_VECTOR_SIZE;

	ASMcategory = -1;
	ASMsupport_group1 = 0.0;
	ASMsupport_group2 = 0.0;

	major_alter_paths = NULL;
	major_alter_paths_num = 0;
}


GenomeTree::GenomeTree()
{
	root = NULL;
}

GenomeTree::~GenomeTree()
{
	//
}

deletedSites::deletedSites()
{
	sites = 0;
	next = NULL;
}


graphpath::graphpath()
{
	direction = undetermined;
	pathlength = 0;
	vertex_start = NULL;
	vertex_end = NULL;
	edgelist = NULL;
	edgelist_tail = NULL;
	pathindex = 0;

	next = NULL;
}

graphpath* graphpath::clone()
{
	//clone a same graph path
	graphpath *result_path = new graphpath;

	result_path->direction = direction;
	result_path->pathlength = pathlength;
	result_path->vertex_start = vertex_start;
	result_path->vertex_end = vertex_end;
	result_path->pathindex = pathindex;

	JuncGraphEdge *cur_edge, *new_edge;
	cur_edge = edgelist;
	while (cur_edge != NULL)
	{
		new_edge = new JuncGraphEdge;
		new_edge->linkedVertex = cur_edge->linkedVertex;

		if (result_path->edgelist == NULL)
		{
			result_path->edgelist = new_edge;
			result_path->edgelist_tail = new_edge;
		}
		else
		{
			result_path->edgelist_tail->next = new_edge;
			result_path->edgelist_tail = new_edge;
		}

		cur_edge = cur_edge->next;
	}

	return result_path;
}

void graphpath::add_edge(JuncGraphEdge *new_edge)
{
	//insert the newfeature to tail of the feature list
	if (edgelist_tail == NULL)
	{
		edgelist = new_edge;
		edgelist_tail = new_edge;
	}
	else
	{
		edgelist_tail->next = new_edge;
		edgelist_tail = new_edge;
	}

	return;
}

graphpath::~graphpath()
{
	//	next = NULL;

	JuncGraphEdge *del_edge;
	while (edgelist != NULL)
	{
		del_edge = edgelist;
		edgelist = del_edge->next;
		delete del_edge;
	}
}

pathCluster::pathCluster()
{
	fragstart = NULL;
	fragend = NULL;
	paths = NULL;
}

pathCluster::~pathCluster()
{
	delete fragstart;
	delete fragend;
	graphpath *del_path;
	while (paths != NULL)
	{
		del_path = paths;
		paths = del_path->next;
		delete del_path;
	}
}

pathlength_record::pathlength_record()
{
	lowerbound = MAX;
	higherbound = 0;
	next = NULL;
}

connectedRead::connectedRead()
{
	readlowIndex = MAX;
	readhighIndex = 0;
	pathselectedIndex = -1;
	lengthlist = NULL;
	minLength_inbetween = MAX;
}

connectedRead::~connectedRead()
{
	pathlength_record *dellist;
	while(lengthlist != NULL)
	{
		dellist = lengthlist;
		lengthlist = dellist->next;
		delete dellist;
	}
}

graphpathCluster::graphpathCluster()
{
	readLowIndex = MAX;
	readHighIndex = 0;
	paths = NULL;
}

graphpathCluster* graphpathCluster::clone()
{
	graphpathCluster *resultcluster;
	graphpath *curpath, *newpath, *pathtail;

	resultcluster = new graphpathCluster;
	resultcluster->readLowIndex = readLowIndex;
	resultcluster->readHighIndex = readHighIndex;
	curpath = paths;
	while(curpath != NULL)
	{
		newpath = curpath->clone();
		if (resultcluster->paths == NULL)
		{
			resultcluster->paths = newpath;
			pathtail = newpath;
		}
		else
		{
			pathtail->next = newpath;
			pathtail = newpath;
		}
		curpath = curpath->next;
	}

	return resultcluster;
}


graphpathCluster::~graphpathCluster()
{
	graphpath *del_path;
	while (paths != NULL)
	{
		del_path = paths;
		paths = del_path->next;
		delete del_path;
	}
}

// Read Graph Structure 

ReadVertex::ReadVertex()
{
	readname = "";
	start = MAX;
	end = 0;
	start_min = MAX;
	end_max = 0;
	readLength = 0;
	Id = -1;
	includedReadNum = 0;
	transDirection = undetermined;
	spliced = false;
	list = NULL;
	listTail = NULL;
	next = NULL;
}

ReadVertex* ReadVertex::clone()
{
	ReadVertex *newVertex;
	rangeJunction *newfrag, *curfrag;

	newVertex = new ReadVertex;
	newVertex->readname = readname;
	newVertex->start = start;
	newVertex->end = end;
	newVertex->start_min = start_min;
	newVertex->end_max = end_max;
	newVertex->readLength = readLength;
	newVertex->Id = Id;
	newVertex->includedReadNum = includedReadNum;
	newVertex->transDirection = transDirection;
	newVertex->spliced = spliced;
	newVertex->junctionsIndex = junctionsIndex;

	curfrag = list;
	while(curfrag != NULL)
	{
		newfrag = new rangeJunction;
		newfrag->junc = curfrag->junc->clone();
		newfrag->next = NULL;

		if (newVertex->list == NULL)
		{
			newVertex->list = newfrag;
			newVertex->listTail = newfrag;
		}
		else
		{
			newVertex->listTail->next = newfrag;
			newVertex->listTail = newfrag;
		}
		curfrag = curfrag->next;
	}
	newVertex->listTail->next = NULL;
	newVertex->next = NULL;

	return newVertex;
}

ReadVertex::~ReadVertex()
{
	next = NULL;

	rangeJunction *del_list;
	while(list != NULL)
	{
		del_list = list;
		list = del_list->next;
		delete del_list;
	}
	junctionsIndex.clear();
}

ReadPath::ReadPath()
{
	pathVertices = NULL;
	firstVertex = NULL;
	lastVertex = NULL;
	transDirection = undetermined;
	spliced = false;
	next = NULL;
}

ReadPath* ReadPath::clone()
{
	ReadPath *resultPath;
	resultPath = new ReadPath;

	ReadVertex *curVertex, *newVertex;
	curVertex = pathVertices;
	while(curVertex != NULL)
	{
		newVertex = curVertex->clone();
		if (resultPath->pathVertices == NULL)
		{
			resultPath->firstVertex = newVertex;
			resultPath->pathVertices = newVertex;
			resultPath->lastVertex = newVertex;
		}
		else
		{
			resultPath->lastVertex->next = newVertex;
			resultPath->lastVertex = newVertex;
		}
		curVertex = curVertex->next;
	}
	resultPath->lastVertex->next = NULL;
	resultPath->transDirection = transDirection;
	resultPath->spliced = spliced;

	return resultPath;
}

ReadPath::~ReadPath()
{
	ReadVertex *del_vertex;
	while(pathVertices != NULL)
	{
		del_vertex = pathVertices;
		pathVertices = del_vertex->next;
		delete del_vertex;
	}
}

ReadIndicesPath :: ReadIndicesPath()
{
	readIndices = NULL; 
	lastreadIndex = NULL; 
	copyNum = 0; 
	transDirection = undetermined;
	next = NULL;
}

ReadIndicesPath :: ~ReadIndicesPath()
{
	ReadIndex *del_Index;
	while(readIndices != NULL)
	{
		del_Index = readIndices;
		readIndices = del_Index->next;
		delete del_Index;
	}
}

Transcript::Transcript()
{
	copyNum = 0;
	readpaths = NULL;
	transcript = NULL;
	nextList = NULL;
}

Transcript* Transcript::clone()
{
	Transcript *newTrans;
	newTrans = new Transcript;

	newTrans->nextList = NULL;
	newTrans->copyNum = copyNum;
	newTrans->transcript = transcript->clone();
	ReadPath *curpath, *newpath;
	curpath = readpaths;
	while(curpath != NULL)
	{
		newpath = curpath->clone();

		newpath->next = newTrans->readpaths;
		newTrans->readpaths = newpath;

		curpath = curpath->next;
	}

	return newTrans;
}

Transcript::~Transcript()
{
	ReadPath *del_path;
	while(readpaths != NULL)
	{
		del_path = readpaths;
		readpaths = del_path->next;
		delete del_path;
	}
	if (transcript != NULL)
	{
		delete transcript;
	}

}

geneExonBoundary::geneExonBoundary()
{
	position = MAX;
	next = NULL;
}


TransInfo::TransInfo()
{
	transLength = 0;
	nextTransinfo = NULL;
}

TransInfo::~TransInfo()
{
	FragIndexCovered.clear();
}

ReadInfo::ReadInfo()
{
	nextreadInfo = NULL;
}

ReadInfo::~ReadInfo()
{
	FragIndexCovered.clear();
}

AdjASMsInfo::AdjASMsInfo()
{
	firstASM_start = MAX;
	secondASM_end = 0;
	Hubstart = MAX;
	Hubend = 0;
	Hublength = 0;
	ASM_Transinfo = NULL;
	deterministicReadinfo = NULL;
	nextinfo = NULL;
}

AdjASMsInfo::~AdjASMsInfo()
{
	ReadInfo *delinfo;
	while(deterministicReadinfo != NULL)
	{
		delinfo = deterministicReadinfo;
		deterministicReadinfo = delinfo->nextreadInfo;
		delete delinfo;
	}
	TransInfo *delinfo_trans;
	while(ASM_Transinfo != NULL)
	{
		delinfo_trans = ASM_Transinfo;
		ASM_Transinfo = delinfo_trans->nextTransinfo;
		delete delinfo_trans;
	}
}

int IdenticalTrans(RangeJunctionList *ref, RangeJunctionList *assembly)
{
	int relation = -1; // 0 indicates trans 2 is different from trans 1; 1 indicates trans1 contains trans2; 2 indicates trans1 is contained within trans1
	// 3 indicates overlapping (trans1 is downstream); 4 indicates overlapping (trans2 is downstream); 5 indicates identical.

	bool overlap;
	rangeJunction *ref_frag, *assembly_frag, *ref_tmpfrag, *assembly_tmpfrag, *ref_prevfrag, *assembly_prevfrag;

	overlap = false;
	ref_prevfrag = assembly_prevfrag = NULL;
	ref_frag = ref->list;
	assembly_frag = assembly->list;
	if (ref_frag != NULL && assembly_frag != NULL)
	{
		if (ref_frag->junc->start < assembly_frag->junc->start)
		{
			while(ref_frag != NULL)
			{
				if (ref_frag == ref->list)
				{
					if (assembly_frag->junc->start > ref_frag->junc->start && abs(assembly_frag->junc->end - ref_frag->junc->end) <= 1)
					{
						overlap = true;
						ref_prevfrag = ref_frag;
						assembly_prevfrag = assembly_frag;
						ref_frag = ref_frag->next;
						assembly_frag = assembly_frag->next;
						break;
					}
				}
				else
				{
					if ((ref_frag->junc->start == assembly_frag->junc->start) && (ref_frag->junc->end == assembly_frag->junc->end))
					{
						// overlap one exonic segment
						if ((ref_frag->next == NULL && ref_frag != ref->list && assembly_frag != assembly->list) || (assembly_frag->next == NULL && ref_frag != ref->list && assembly_frag != assembly->list))
						{
							overlap = false;
						}
						else
						{
							ref_prevfrag = ref_frag;
							assembly_prevfrag = assembly_frag;
							ref_frag = ref_frag->next;
							assembly_frag = assembly_frag->next;
							overlap = true;
						}
						break;
					}
					else if ((ref_frag->junc->start <= assembly_frag->junc->start) && (ref_frag->junc->end >= assembly_frag->junc->end) && assembly_frag == assembly->list && assembly_frag->next == NULL)
					{
						if (ref_frag->junc->type == frag_exon || ref_frag->junc->type == frag_retained_intron)
						{
							overlap = true;
							ref_prevfrag = ref_frag;
							assembly_prevfrag = assembly_frag;
							ref_frag = ref_frag->next;
							assembly_frag = assembly_frag->next;
						}
						else
						{
							overlap = false;
						}
						break;
					}
				}

				ref_frag = ref_frag->next;
			}
		}
		else if (ref_frag->junc->start > assembly_frag->junc->start)
		{
			while(assembly_frag != NULL)
			{
				if (assembly_frag == assembly->list && ref_frag == ref->list)
				{
					if (ref_frag->junc->start > assembly_frag->junc->start && ref_frag->junc->end == assembly_frag->junc->end)
					{
						overlap = true;
						ref_prevfrag = ref_frag;
						assembly_prevfrag = assembly_frag;
						ref_frag = ref_frag->next;
						assembly_frag = assembly_frag->next;
						break;
					}
				}
				else
				{
					if ((ref_frag->junc->start == assembly_frag->junc->start) && (ref_frag->junc->end == assembly_frag->junc->end))
					{
						// overlap one exonic segment
						if ((ref_frag->next == NULL && ref_frag != ref->list && assembly_frag != assembly->list) || (assembly_frag->next == NULL && ref_frag != ref->list && assembly_frag != assembly->list))
						{
							overlap = false;
						}
						else
						{
							ref_prevfrag = ref_frag;
							assembly_prevfrag = assembly_frag;
							ref_frag = ref_frag->next;
							assembly_frag = assembly_frag->next;
							overlap = true;
						}
						break;
					}
				}
				assembly_frag = assembly_frag->next;
			}
		}
		else if ((ref_frag->junc->start == assembly_frag->junc->start) && (ref_frag->junc->end == assembly_frag->junc->end))
		{
			// overlap one exonic segment
			if ((ref_frag->next == NULL && ref_frag != ref->list && assembly_frag != assembly->list) || (assembly_frag->next == NULL && ref_frag != ref->list && assembly_frag != assembly->list))
			{
				overlap = false;
			}
			else
			{
				ref_prevfrag = ref_frag;
				assembly_prevfrag = assembly_frag;
				ref_frag = ref_frag->next;
				assembly_frag = assembly_frag->next;
				overlap = true;
			}
		}
		else if (abs(ref_frag->junc->start - assembly_frag->junc->start) <= 1 && abs(ref_frag->junc->end - assembly_frag->junc->end) <= 1 && ref_frag->next == NULL && assembly_frag->next == NULL)
		{
			overlap = true;
			relation = 5;
			return relation;
		}
		else
		{
			overlap = false;
		}
	}

	if (overlap == true)
	{
		while(ref_frag != NULL && assembly_frag != NULL)
		{
			if (assembly_frag->next == NULL)
			{
				if (ref_frag->junc->start == assembly_frag->junc->start)
				{
					// do nothing
				}
				else
				{
					relation = 0;
					break;
				}
			}
			else
			{
				if ((ref_frag->junc->start != assembly_frag->junc->start) || (ref_frag->junc->end != assembly_frag->junc->end))
				{
					relation = 0;
					break;
				}
			}

			ref_prevfrag = ref_frag;
			assembly_prevfrag = assembly_frag;
			ref_frag = ref_frag->next;
			assembly_frag = assembly_frag->next;
		}

		if (relation != 0)
		{
			ref_tmpfrag = ref->list;
			assembly_tmpfrag = assembly->list;
			// 0 indicates trans 2 is different from trans 1; 1 indicates trans1 contains trans2; 2 indicates trans1 is contained within trans1
			// 3 indicates overlapping (trans1 is downstream); 4 indicates overlapping (trans2 is downstream); 5 indicates identical.
			if (ref_frag == NULL && assembly_frag == NULL)
			{
				if (ref_tmpfrag->junc->start == assembly_tmpfrag->junc->start && ref_tmpfrag->junc->end == assembly_tmpfrag->junc->end)
				{
					relation = 5;
				}
				else if (ref_tmpfrag->junc->start < assembly_tmpfrag->junc->start)
				{
					if (abs(ref_tmpfrag->junc->end - assembly_tmpfrag->junc->end) <= 1)
					{
						relation = 5;
					}
					else
					{
						relation = 1;
					}

				}
				else if (ref_tmpfrag->junc->start > assembly_tmpfrag->junc->start)
				{
					if (ref_tmpfrag->junc->end == assembly_tmpfrag->junc->end)
					{
						relation = 5;
					}
					else
					{
						relation = 2;
					}
				}
			}
			else if (ref_frag == NULL)
			{
				if (abs(assembly_prevfrag->junc->end - ref_prevfrag->junc->end) <= 1)
				{
					if (ref_tmpfrag->junc->start < assembly_tmpfrag->junc->start)
					{
						relation = 3;
					}
					else
					{
						relation = 2;
					}
				}
				else
				{
					relation = 0;
				}
			}
			else if (assembly_frag == NULL)
			{
				if (abs(assembly_prevfrag->junc->end - ref_prevfrag->junc->end) <= 1)
				{
					if (ref_tmpfrag->junc->start > assembly_tmpfrag->junc->start)
					{
						relation = 4;
					}
					else
					{
						relation = 1;
					}
				}
				else
				{
					relation = 0;
				}
			}
		}
	}
	else
	{
		relation = 0;
	}

	return relation;
}

// Transform long to string
string itostr(long t)
{
	ostringstream oss;
	oss << t;
	return oss.str();
}


// Merge adjacent exons of each assembled transcript
void MergeAdjExons(RangeJunctionList *prematurelist)
{
	rangeJunction *curfrag, *prevfrag, *newfrag;

	prevfrag = NULL;
	curfrag = prematurelist->list;
	while(curfrag->next != NULL)
	{
		if (curfrag->junc->type != frag_junction && curfrag->next->junc->type != frag_junction && (curfrag->next->junc->start - curfrag->junc->end < 2) )
		{
			// merge two adjacent exons
			newfrag = new rangeJunction;
			newfrag->junc = new fragment;
			newfrag->junc->start = curfrag->junc->start;
			newfrag->junc->end = curfrag->next->junc->end;
			newfrag->junc->type = frag_exon;
			newfrag->junc->transDirection = curfrag->junc->transDirection;
			if (prevfrag == NULL)
			{
				newfrag->next = curfrag->next->next;
				prematurelist->list = newfrag;
			}
			else
			{
				prevfrag->next = newfrag;
				newfrag->next = curfrag->next->next;
			}
			delete curfrag->next;
			delete curfrag;
			curfrag = newfrag;
		}
		else
		{
			prevfrag = curfrag;
			curfrag = curfrag->next;
		}
	}
	return;
}

/************************************************************************/
/* Build GTree */
/************************************************************************/

//clear a vector and release its memory
template <class T>
void free_vector(T &t)
{
	T tmp;
	tmp.swap(t);
}

void initialization()
{
	string filename, info;

	//read statFile
	ifstream statfile;
	filename = inputPath + "stat.txt";
	statfile.open(filename.c_str());


	long tmp = -1, tmp_small = MAX_CHR_LENGTH, tmp_large = -1;
	statfile >> tmp;
	while (tmp >= 0)
	{
//		statfile >> tmp;

		if (tmp < tmp_small)
		{
			tmp_small = tmp;
		}

		statfile >> tmp;
		if (tmp > tmp_large)
		{
			tmp_large = tmp;
		}

		getline(statfile, info);

		tmp = -1;
		statfile >> tmp;
	}
	CHROMOSOME_START = tmp_small;
	CHROMOSOME_END = tmp_large;
	statfile.close();

	srand ( time(NULL) );


	//read dataset stat file
	for (int tmpCnt = 0; tmpCnt < SUPPORT_VECTOR_SIZE; tmpCnt++)
	{
		DatasetReadCount_total[tmpCnt] = 0;
		DatasetNormalizationRatio[tmpCnt] = 1.0;
	}
	ifstream datasetStatFile;
	filename = inputPath + "DatasetStat.txt";
	datasetStatFile.open(filename.c_str());
	for (int tmpCnt = 0; tmpCnt < SUPPORT_VECTOR_SIZE; tmpCnt++)
	{
		datasetStatFile >> DatasetReadCount_total[tmpCnt];
	}
	datasetStatFile.close();

	return;
}

void cleanAll()
{
	// For assembly
	inputfile_read.close();
	// 	outputfile_gtf.close();
	// 	outputfile_stat.close();
	// 	outputfile_readpath.close();
	// 	outputfile_gaps.close();

	return;
}


void input_junction(string junctionFilename)
{
	ifstream junctionFile;
	junctionFile.open(junctionFilename.c_str());

	fragment *newJunction;
	alter_junction *newAlterJunc;
	rangeJunction *newRangeJunc;
	char name[30];
	bool junctionExist, thresholdValid, strand;
	long temp, support1, support2;
	double tmpSupport[8];
	string info;
	for (int i = 0; i < 30; i++)
	{
		name[i] = '\0';
	}

	junctionFile >> name;
	while (name[1] != '\0')
	{
		newJunction = new fragment;
		newJunction->type = frag_junction;

		///* for fragment db
		strcpy(newJunction->chromosome_start, name);
		strcpy(newJunction->chromosome_end, name);
		junctionFile >> newJunction->start;
		junctionFile >> newJunction->end;
		junctionFile >> strand; newJunction->transDirection = strand == true ? sense : antisense;
		getline(junctionFile, info);

#ifndef TRANSCRIPTION_DIRECTION
		newJunction->transDirection = sense;
#endif

		thresholdValid = true;

		if (newJunction->end - newJunction->start - 1 < MINJUNCTIONLENGTH)
			thresholdValid = false;

		if (thresholdValid == true)
		{

			junctionExist = false;
			for (temp = sortList_Junction_Num>0?sortList_Junction_Num:1; temp <= sortList_Junction_Num; temp++)
			{
				if (strcmp(sortList_Junction[temp]->junc->chromosome_start, newJunction->chromosome_start) == 0 && strcmp(sortList_Junction[temp]->junc->chromosome_end, newJunction->chromosome_end) == 0 && sortList_Junction[temp]->junc->start == newJunction->start && sortList_Junction[temp]->junc->end == newJunction->end)
				{
					if (sortList_Junction[temp]->junc->transDirection == newJunction->transDirection)
					{
						junctionExist = true;
						delete newJunction;
						break;
					}
					else
					{
						CLEAN_SPLICE_DIRECTION = true;
						cout << "Warning: two junctions at the same location but of different transcription direction: " << newJunction->start << " - " << newJunction->end << "\n";
					}
				}
			}

			if (junctionExist == false)
			{
				newJunction->ID = ++fragmentID_Cnt;
				newRangeJunc = new rangeJunction;
				newRangeJunc->junc = newJunction;

				sortList_Junction_Num++;
				junctionNum++;
				sortList_Junction[sortList_Junction_Num] = newRangeJunc;
			}
		}

		name[1] = '\0';
		junctionFile >> name;	
	}

	junctionFile.close();

	return;
}


void getJunctionSupport(int index)
{
	//get junction support from fragent file
	ifstream junctionSupportFile;
	string inputfilename;

	inputfilename = inputPath + "junction_" + itostr(index+1) + ".txt";
	junctionSupportFile.open(inputfilename.c_str());

	long i;
	for (i = 1; i <= sortList_Junction_Num; i++)
	{
		sortKey_Junction[i] = sortList_Junction[i]->junc->end;
	}
	mergeSort_JunctionSort(sortList_Junction_Num);

	for (i = 1; i <= sortList_Junction_Num; i++)
	{
		sortKey_Junction[i] = sortList_Junction[i]->junc->start;
	}
	mergeSort_JunctionSort(sortList_Junction_Num);


	long startPosition, endPosition, juncIndex, tmpCount;
	char chromosome[100];
	string info;
	bool strand;
	startPosition = 0;
	juncIndex = 1;

	chromosome[0] = '\0';
	junctionSupportFile >> chromosome;
	while (chromosome[0] != '\0')
	{
		junctionSupportFile >> startPosition;
		junctionSupportFile >> endPosition;
		junctionSupportFile >> strand;

#ifndef TRANSCRIPTION_DIRECTION
		strand = true;
#endif

		getline(junctionSupportFile, info);

		while (juncIndex <= sortList_Junction_Num && startPosition > sortList_Junction[juncIndex]->junc->start)
		{
			juncIndex++;
		}

		while (juncIndex <= sortList_Junction_Num && startPosition == sortList_Junction[juncIndex]->junc->start && endPosition > sortList_Junction[juncIndex]->junc->end)
		{
			juncIndex++;
		}

		if (juncIndex > sortList_Junction_Num)
		{
			break;
		}

		for (tmpCount = 0; juncIndex + tmpCount <= sortList_Junction_Num && startPosition == sortList_Junction[juncIndex+tmpCount]->junc->start && endPosition == sortList_Junction[juncIndex+tmpCount]->junc->end; ++tmpCount)
		{
			if ((strand == true ? sense : antisense) == sortList_Junction[juncIndex+tmpCount]->junc->transDirection)
			{
				((sortList_Junction[juncIndex+tmpCount]->junc->support)[index]) += 1;
			}
		}

		chromosome[0] = '\0';
		junctionSupportFile >> chromosome;
	}

	junctionSupportFile.close();

	return;
}


void cleanSpliceDirection()
{
	//merge splice junctions with same location but different strands
	if (CLEAN_SPLICE_DIRECTION == false || sortList_Junction_Num < 1)
		return;

	long iLoop;
	int indexLoop;
	double curSupp, nextSupp;

	rangeJunction *junctionList = NULL, *curJunc, *delJunc;
	for (iLoop = sortList_Junction_Num; iLoop >= 1; --iLoop)
	{
		sortList_Junction[iLoop]->next = junctionList;
		junctionList = sortList_Junction[iLoop];
	}

	curJunc = junctionList;
	while (curJunc->next != NULL)
	{
		if (curJunc->junc->start == curJunc->next->junc->start && curJunc->junc->end == curJunc->next->junc->end)
		{
			curSupp = 0.0; nextSupp = 0.0;
			for (indexLoop = 0; indexLoop < SUPPORT_VECTOR_SIZE; ++indexLoop)
			{
				curSupp += curJunc->junc->support[indexLoop];
				nextSupp += curJunc->next->junc->support[indexLoop];
			}

			curJunc->junc->transDirection = curSupp > nextSupp ? curJunc->junc->transDirection : curJunc->next->junc->transDirection;
			for (indexLoop = 0; indexLoop < SUPPORT_VECTOR_SIZE; ++indexLoop)
				curJunc->junc->support[indexLoop] += curJunc->next->junc->support[indexLoop];

			delJunc = curJunc->next;
			curJunc->next = curJunc->next->next;
			delete delJunc;
			--sortList_Junction_Num;
		}
		else
			curJunc = curJunc->next;
	}

	for (curJunc = junctionList, iLoop = 0; curJunc != NULL; curJunc = curJunc->next)
		sortList_Junction[++iLoop] = curJunc;

	if (iLoop != sortList_Junction_Num)
	{
		cout << "Number of splice junction does not match in cleanSpliceDirection" << endl;
		exit(1);
	}

	return;
}


void filterJunction(double thresh_maxValue, int thresh_zeroCnt, double thresh_meanValue)
{
	long iLoop, jLoop;
	int indexLoop, zeroCnt;
	bool keep;
	double total_support;

	for (iLoop = 1; iLoop <= sortList_Junction_Num; iLoop++)
	{
		keep = false;
		zeroCnt = 0;
		total_support = 0.0;

		for (indexLoop = 0; indexLoop < SUPPORT_VECTOR_SIZE; indexLoop++)
		{
			if (sortList_Junction[iLoop]->junc->support[indexLoop] > thresh_maxValue)
			{
				keep = true;
			}
			if (sortList_Junction[iLoop]->junc->support[indexLoop] == 0)
			{
				zeroCnt++;
			}
			total_support += sortList_Junction[iLoop]->junc->support[indexLoop];
		}
		if (keep == false || zeroCnt >= thresh_zeroCnt || total_support/SUPPORT_VECTOR_SIZE < thresh_meanValue)
		{
			//throw away
			delete sortList_Junction[iLoop]->junc;
			delete sortList_Junction[iLoop];

			sortList_Junction_Num--;
			for (jLoop = iLoop; jLoop <= sortList_Junction_Num; jLoop++)
			{
				sortList_Junction[jLoop] = sortList_Junction[jLoop + 1];
			}
			iLoop--;
		}
	}

	return;
}

void filterIntron(long start_sortList_Junction_index)
{
	long iLoop, jLoop;
	bool keep;
	int sampleLoop;
	double totalCoverage;

	for (iLoop = start_sortList_Junction_index; iLoop <= sortList_Junction_Num; ++iLoop)
	{
		if (sortList_Junction[iLoop]->junc->type == frag_retained_intron)
		{
			keep = false;
			totalCoverage = 0.0;
			for (sampleLoop = 0; sampleLoop < SUPPORT_VECTOR_SIZE; ++sampleLoop)
				totalCoverage += sortList_Junction[iLoop]->junc->support[sampleLoop];
			if (totalCoverage/SUPPORT_VECTOR_SIZE >= coverageThreshold_intron)
				keep = true;
		}
		else
			keep = true;

		if (keep == false)
		{
			//throw away
			delete sortList_Junction[iLoop]->junc;
			delete sortList_Junction[iLoop];

			sortList_Junction_Num--;
			for (jLoop = iLoop; jLoop <= sortList_Junction_Num; jLoop++)
			{
				sortList_Junction[jLoop] = sortList_Junction[jLoop + 1];
			}
			iLoop--;
		}
	}

	return;
}

void filterFalseFragments(bool filterStandaloneFragment)
{
	//filter standalone fragments (exons and junctions)
	if (sortList_Junction_Num <= 1)
	{
		//		cout << "Warning: less than 2 fragments." << endl;
		return;
	}

	RangeJunctionList *resultList;
	rangeJunction *curJunc, *prevJunc, *delJunc, *tmpJunc;
	double maxSupport;
	long jLoop = 1;
	bool tobedeleted;

	//first, both ends of a junction must appear. filter false ones
	//filter junctions with no start exon
	for (jLoop = 1; jLoop <= sortList_Junction_Num; jLoop++)
	{
		sortKey_Junction[jLoop] = sortList_Junction[jLoop]->junc->end;
	}
	mergeSort_JunctionSort(sortList_Junction_Num);

	resultList = new RangeJunctionList;	
	for (jLoop = 1; jLoop <= sortList_Junction_Num; jLoop++)
	{
		sortList_Junction[jLoop]->next = resultList->list;
		resultList->list = sortList_Junction[jLoop];
	}

	curJunc = resultList->list;
	prevJunc = NULL;
	while (curJunc != NULL)
	{
		if (curJunc->junc->type == frag_junction)
		{
			tobedeleted = true;
			tmpJunc = curJunc->next;
			while (tmpJunc != NULL && tmpJunc->junc->end >= curJunc->junc->start)
			{
				if ((tmpJunc->junc->type == frag_exon || tmpJunc->junc->type == frag_retained_intron) && tmpJunc->junc->end == curJunc->junc->start)
				{
					tobedeleted = false;
					break;
				}
				tmpJunc = tmpJunc->next;
			}
			if (tobedeleted == true)
			{
				if (prevJunc != NULL)
				{
					prevJunc->next = curJunc->next;
					delete curJunc;
					sortList_Junction_Num--;
					curJunc = prevJunc->next;
				}
				else
				{
					resultList->list = curJunc->next;
					delete curJunc;
					sortList_Junction_Num--;
					curJunc = resultList->list;
				}
			}
			else
			{
				prevJunc = curJunc;
				curJunc = curJunc->next;
			}
		}
		else
		{
			prevJunc = curJunc;
			curJunc = curJunc->next;
		}
	}

	for (jLoop = 1, curJunc = resultList->list; curJunc != NULL; jLoop++, curJunc = curJunc->next)
	{
		sortList_Junction[jLoop] = curJunc;
	}
	if (jLoop != sortList_Junction_Num + 1)
	{
		cout << "Error: filter false junctions" << endl;
		exit(1);
	}
	resultList->list = NULL;

	//filter junctions with no end exon
	for (jLoop = 1; jLoop <= sortList_Junction_Num; jLoop++)
	{
		sortKey_Junction[jLoop] = sortList_Junction[jLoop]->junc->start;
	}
	mergeSort_JunctionSort(sortList_Junction_Num);

	for (jLoop = sortList_Junction_Num; jLoop >= 1; jLoop--)
	{
		sortList_Junction[jLoop]->next = resultList->list;
		resultList->list = sortList_Junction[jLoop];
	}

	curJunc = resultList->list;
	prevJunc = NULL;
	while (curJunc != NULL)
	{
		if (curJunc->junc->type == frag_junction)
		{
			tobedeleted = true;
			tmpJunc = curJunc->next;
			while (tmpJunc != NULL && tmpJunc->junc->start <= curJunc->junc->end)
			{
				if ((tmpJunc->junc->type == frag_exon || tmpJunc->junc->type == frag_retained_intron) && tmpJunc->junc->start == curJunc->junc->end)
				{
					tobedeleted = false;
					break;
				}
				tmpJunc = tmpJunc->next;
			}
			if (tobedeleted == true)
			{
				if (prevJunc != NULL)
				{
					prevJunc->next = curJunc->next;
					delete curJunc;
					sortList_Junction_Num--;
					curJunc = prevJunc->next;
				}
				else
				{
					resultList->list = curJunc->next;
					delete curJunc;
					sortList_Junction_Num--;
					curJunc = resultList->list;
				}
			}
			else
			{
				prevJunc = curJunc;
				curJunc = curJunc->next;
			}
		}
		else
		{
			prevJunc = curJunc;
			curJunc = curJunc->next;
		}
	}

	for (jLoop = 1, curJunc = resultList->list; curJunc != NULL; jLoop++, curJunc = curJunc->next)
	{
		sortList_Junction[jLoop] = curJunc;
	}
	if (jLoop != sortList_Junction_Num + 1)
	{
		cout << "Error: filter false junctions" << endl;
		exit(1);
	}

	resultList->list = NULL;


	//second, filter standalone & low-expression exons
	if (filterStandaloneFragment == true)
	{
		for (jLoop = 1; jLoop <= sortList_Junction_Num; jLoop++)
		{
			sortKey_Junction[jLoop] = sortList_Junction[jLoop]->junc->end;
		}
		mergeSort_JunctionSort(sortList_Junction_Num);

		for (jLoop = 1; jLoop <= sortList_Junction_Num; jLoop++)
		{
			sortKey_Junction[jLoop] = sortList_Junction[jLoop]->junc->start;
		}
		mergeSort_JunctionSort(sortList_Junction_Num);

		for (jLoop = sortList_Junction_Num; jLoop >= 1; jLoop--)
		{
			sortList_Junction[jLoop]->next = resultList->list;
			resultList->list = sortList_Junction[jLoop];
		}	

		jLoop = 1;
		curJunc = resultList->list;
		prevJunc = NULL;
		while (curJunc != NULL)
		{
			tobedeleted = false;
			delJunc = NULL;

			if (prevJunc != NULL && curJunc->next != NULL)
			{
				if (abs(prevJunc->junc->end	- curJunc->junc->start) <= 1 || abs(prevJunc->junc->start - curJunc->junc->start) <= 1 || abs(prevJunc->junc->end - curJunc->junc->end) <= 1 
					|| abs(curJunc->junc->end - curJunc->next->junc->start) <= 1 || abs(curJunc->junc->end - curJunc->next->junc->end) <= 1 || abs(curJunc->junc->start - curJunc->next->junc->start) <= 1)
				{//keep
					prevJunc = curJunc; 
					sortList_Junction[jLoop++] = curJunc;
					curJunc = curJunc->next;
				}
				else
				{
					prevJunc->next = curJunc->next;
					delJunc = curJunc;
					tobedeleted = true;
					sortList_Junction_Num--;
					curJunc = prevJunc->next;
				}
			}
			else if (prevJunc == NULL && curJunc->next != NULL)
			{
				if (abs(curJunc->junc->end - curJunc->next->junc->start) <= 1 || abs(curJunc->junc->end - curJunc->next->junc->end) <= 1 || abs(curJunc->junc->start - curJunc->next->junc->start) <= 1)
				{//keep
					prevJunc = curJunc; 
					sortList_Junction[jLoop++] = curJunc;
					curJunc = curJunc->next;
				}
				else
				{
					resultList->list = curJunc->next;
					delJunc = curJunc;
					tobedeleted = true;
					sortList_Junction_Num--;
					curJunc = resultList->list;
				}
			}
			else if (prevJunc != NULL && curJunc->next == NULL)
			{
				if (abs(prevJunc->junc->end	- curJunc->junc->start) <= 1 || abs(prevJunc->junc->start - curJunc->junc->start) <= 1 || abs(prevJunc->junc->end - curJunc->junc->end) <= 1)
				{//keep
					prevJunc = curJunc; 
					sortList_Junction[jLoop++] = curJunc;
					curJunc = curJunc->next;
				}
				else
				{
					prevJunc->next = curJunc->next;
					delJunc = curJunc;
					tobedeleted = true;
					sortList_Junction_Num--;
					curJunc = prevJunc->next;
				}
			}
			else
			{
				//only one fragment, filter by coverage
				// 				if (maxSupport > MIN_COVERAGE_FILTER_FRAGMENT)
				// 				{//keep
				// 					prevJunc = curJunc; 
				// 					sortList_Junction[jLoop++] = curJunc;
				// 					curJunc = curJunc->next;
				// 				}
				//				else
				{
					resultList->list = curJunc->next;
					delete curJunc;
					sortList_Junction_Num--;
					curJunc = resultList->list;
				}
			}

			if (tobedeleted == true && delJunc != NULL)
			{
				delete delJunc;
			}
		}

		for (jLoop = 1, curJunc = resultList->list; curJunc != NULL; jLoop++, curJunc = curJunc->next)
		{
			sortList_Junction[jLoop] = curJunc;
		}
		if (jLoop != sortList_Junction_Num + 1)
		{
			cout << "Error: filter false exons" << endl;
			exit(1);
		}
	}

	resultList->list = NULL;
	delete resultList;

	return;
}

//sort splice sites
void merge_SpliceSiteSort(long p, long q, long r)
{
	long n1, n2, i, j, k;

	n1 = q - p + 1;
	n2 = r - q;

	for (i = 1; i <= n1; i++)
	{
		mergeSort_Larray[i] = sortKey_spliceSite[p + i - 1];
		mergeSort_LorderedList_spliceSite[i] = sortList_spliceSite[p + i - 1];
	}
	for (j = 1; j <= n2; j++)
	{
		mergeSort_Rarray[j] = sortKey_spliceSite[q + j];
		mergeSort_RorderedList_spliceSite[j] = sortList_spliceSite[q + j];
	}

	mergeSort_Larray[n1 + 1] = MAX_CHR_LENGTH * 2;
	mergeSort_Rarray[n2 + 1] = MAX_CHR_LENGTH * 2;

	i = 1;
	j = 1;

	for (k = p; k <= r; k++)
	{
		if (mergeSort_Larray[i] <= mergeSort_Rarray[j])
		{
			sortKey_spliceSite[k] = mergeSort_Larray[i];
			sortList_spliceSite[k] = mergeSort_LorderedList_spliceSite[i];

			i++;
		} 
		else
		{
			sortKey_spliceSite[k] = mergeSort_Rarray[j];
			sortList_spliceSite[k] = mergeSort_RorderedList_spliceSite[j];

			j++;
		}
	}

	return;
}


void mergeSort_SpliceSiteSort(long sortList_size)
{
	//non-recursive merge sort for sorting junctions
	long m, n, i, r;
	m = 1;
	n = sortList_size;

	while (m <= n)
	{
		i = 1;
		while (i <= n - m)
		{
			r = (i + 2 * m - 1) < n ? (i + 2 * m - 1) : n;
			merge_SpliceSiteSort(i, i + m - 1, r);
			i = i + 2 * m;
		}

		m = m * 2;
	}

	return;
}

void getExons()
{
	spliceSite *curSite;
	fragment *newFragment;
	rangeJunction *newRangeJunc;
	long i;

	if (sortList_Junction_Num == 0)
	{
		newFragment = new fragment;
		newFragment->type = frag_exon;
		newFragment->altersite = 0;
		sprintf(newFragment->name, "Exon_1");
		newFragment->ID = ++fragmentID_Cnt;

		curSite = sortList_spliceSite[1];
		strcpy(newFragment->chromosome_start, resultNamePrefix.c_str());
		strcpy(newFragment->chromosome_end, resultNamePrefix.c_str());
		newFragment->start = CHROMOSOME_START;
		newFragment->end = CHROMOSOME_END;

		//insert exonic fragments to list
		newRangeJunc = new rangeJunction;
		newRangeJunc->junc = newFragment;

		sortList_Junction_Num++;
		sortList_Junction[sortList_Junction_Num] = newRangeJunc;

		return;
	}

	//input splice sites
	sortList_spliceSite_Num = 0;
	for (i = 1; i <= sortList_Junction_Num; i++)
	{
		curSite = new spliceSite;		
		strcpy(curSite->chromosome, sortList_Junction[i]->junc->chromosome_start);
		curSite->position = sortList_Junction[i]->junc->start;
		curSite->directionOut = true;
		sortList_spliceSite_Num++;
		sortList_spliceSite[sortList_spliceSite_Num] = curSite;

		curSite = new spliceSite;		
		strcpy(curSite->chromosome, sortList_Junction[i]->junc->chromosome_end);
		curSite->position = sortList_Junction[i]->junc->end;
		curSite->directionOut = false;
		sortList_spliceSite_Num++;
		sortList_spliceSite[sortList_spliceSite_Num] = curSite;
	}


	//sort splice sites
	for (i = 1; i <= sortList_spliceSite_Num; i++)
	{
		sortKey_spliceSite[i] = sortList_spliceSite[i]->position;
	}
	mergeSort_SpliceSiteSort(sortList_spliceSite_Num);


	//build exonic fragments from adjacent splice sites 
	//exonic if not out-in

	//first exon
	newFragment = new fragment;
	newFragment->type = frag_exon;
	newFragment->altersite = 0;
	sprintf(newFragment->name, "Exon_1");
	newFragment->ID = ++fragmentID_Cnt;

	curSite = sortList_spliceSite[1];
	strcpy(newFragment->chromosome_start, curSite->chromosome);
	strcpy(newFragment->chromosome_end, curSite->chromosome);
	newFragment->start = CHROMOSOME_START;
	newFragment->end = curSite->position;

	//insert exonic fragments to list
	newRangeJunc = new rangeJunction;
	newRangeJunc->junc = newFragment;

	sortList_Junction_Num++;
	sortList_Junction[sortList_Junction_Num] = newRangeJunc;

	for (i = 1; i < sortList_spliceSite_Num; i++)
	{
		if (sortList_spliceSite[i]->position == sortList_spliceSite[i+1]->position)
		{
			curSite = sortList_spliceSite[i];
			delete curSite;
		} 
		else
		{
			if (sortList_spliceSite[i]->directionOut == true && sortList_spliceSite[i+1]->directionOut == false)
			{
				curSite = sortList_spliceSite[i];
				delete curSite;
			} 
			else
			{
				newFragment = new fragment;
				newFragment->type = frag_exon;

				if (sortList_spliceSite[i]->directionOut == false && sortList_spliceSite[i+1]->directionOut == true)
					newFragment->altersite = 0;
				else if (sortList_spliceSite[i]->directionOut == true && sortList_spliceSite[i+1]->directionOut == true)
					newFragment->altersite = 1;
				else if (sortList_spliceSite[i]->directionOut == false && sortList_spliceSite[i+1]->directionOut == false)
					newFragment->altersite = 2;

				sprintf(newFragment->name, "Exon_%ld", i/2 + 2);
				newFragment->ID = ++fragmentID_Cnt;

				curSite = sortList_spliceSite[i];
				strcpy(newFragment->chromosome_start, curSite->chromosome);
				if (curSite->directionOut == true)
				{
					newFragment->start = curSite->position + 1;
				} 
				else
				{
					newFragment->start = curSite->position;
				}
				delete curSite;

				curSite = sortList_spliceSite[i+1];
				strcpy(newFragment->chromosome_end, curSite->chromosome);
				if (curSite->directionOut == true)
				{
					newFragment->end = curSite->position;
				} 
				else
				{
					newFragment->end = curSite->position - 1;
				}

				//insert exonic fragments to list
				newRangeJunc = new rangeJunction;
				newRangeJunc->junc = newFragment;

				sortList_Junction_Num++;
				sortList_Junction[sortList_Junction_Num] = newRangeJunc;
			}
		}
	}

	//last exon
	newFragment = new fragment;
	newFragment->type = frag_exon;
	newFragment->altersite = 0;
	sprintf(newFragment->name, "Exon_%ld", sortList_spliceSite_Num/2 + 2);
	newFragment->ID = ++fragmentID_Cnt;

	curSite = sortList_spliceSite[sortList_spliceSite_Num];
	strcpy(newFragment->chromosome_start, curSite->chromosome);
	strcpy(newFragment->chromosome_end, curSite->chromosome);
	newFragment->start = curSite->position;
	newFragment->end = CHROMOSOME_END;
	delete curSite;

	//insert exonic fragments to list
	newRangeJunc = new rangeJunction;
	newRangeJunc->junc = newFragment;

	sortList_Junction_Num++;
	sortList_Junction[sortList_Junction_Num] = newRangeJunc;

	return;
}


double cutAlterSite(long *coverageVector, long vectorLength, long startIndex, int altersitetype, long &totalSupport, long &cutFragSupport, long &cutFragStart, long &cutFragEnd)
{
	double maxRatio = 0, supp_left, supp_right, curRatio;
	long size_left, size_right, curPos;
	totalSupport = 0; cutFragSupport = 0; cutFragStart = startIndex; cutFragEnd = startIndex + vectorLength - 1;

	if (coverageVector == NULL || altersitetype == 0)
		return 0;
	else if (altersitetype == 1)
	{
		for (size_left = 0, supp_left = 0.0; size_left < vectorLength - 1; ++size_left)
			supp_left += coverageVector[size_left];

		size_right = 1;
		supp_right = coverageVector[vectorLength - 1];

		totalSupport = supp_left + supp_right;

		for (curPos = vectorLength - 1; curPos > 0; --curPos)
		{
			if (supp_left == 0.0)
			{
				if (curPos < MAX_NOCOVERAGE_LENGTH)
				{
					//tolerate this gap
					return maxRatio>1? maxRatio : 0;
				} 
				else
				{
					//return this position
					cutFragStart = startIndex + curPos;
					cutFragSupport = supp_right;
					return COVERAGE_CHANGE_THRESH + 1; //definitely accept
				}				
			}
			else
				curRatio = (supp_right / size_right) / (supp_left / size_left);

			if (curRatio > maxRatio)
			{
				maxRatio = curRatio;
				cutFragStart = startIndex + curPos;				
				cutFragSupport = supp_right;
			}

			--size_left;
			supp_left -= coverageVector[curPos - 1];
			++size_right;
			supp_right += coverageVector[curPos - 1];
		}

		return maxRatio>1? maxRatio : 0;
	}
	else if (altersitetype == 2)
	{
		size_left = 0;
		supp_left = 0.0;

		for (size_right = 0, supp_right = 0.0; size_right < vectorLength; ++size_right)
			supp_right += coverageVector[size_right];

		totalSupport = supp_left + supp_right;

		for (curPos = 0; curPos < vectorLength; ++curPos)
		{
			if (supp_right == 0.0)
			{
				if (vectorLength - curPos <= MAX_NOCOVERAGE_LENGTH)
				{
					//tolerate this gap
					return maxRatio>1? maxRatio : 0;
				} 
				else
				{
					//return this position
					cutFragEnd = startIndex + curPos - 1;
					cutFragSupport = supp_left;
					return COVERAGE_CHANGE_THRESH + 1; //definitely accept
				}				
			}
			else if (size_left == 0)
				curRatio = 0;
			else
				curRatio = (supp_left / size_left) / (supp_right / size_right);

			if (curRatio > maxRatio)
			{
				maxRatio = curRatio;
				cutFragEnd = startIndex + curPos - 1;
				cutFragSupport = supp_left;
			}

			++size_left;
			supp_left += coverageVector[curPos];
			--size_right;
			supp_right -= coverageVector[curPos];
		}

		return maxRatio>1? maxRatio : 0;
	}
	else
	{
		cout << "Abnormal alter site type" << endl;
		exit(1);
	}
}


double windowChangeRatio(long *coverageVector, long vectorLength, long curPos)
{
	double ratio = 1.0;

	if (coverageVector == NULL)
		return 1.0;
	else if (curPos <= 0)
		return 2 * COVERAGE_CHANGE_THRESH;
	else if (curPos >= vectorLength - 1)
		return 0.0;

	double curSupp_left = 0, curSupp_right = 0, maxSupp_left, minSupp_left, maxSupp_right, minSupp_right; //support of the left window and the right window
	int windowSize_left = 0, windowSize_right = 0; //size of the left window and the right window 

	for (maxSupp_left = minSupp_left = coverageVector[curPos - 1]; windowSize_left < COVERAGE_CHANGE_WINDOW && curPos - windowSize_left > 0; ++windowSize_left)
	{
		curSupp_left += coverageVector[curPos - 1 - windowSize_left];
		if (coverageVector[curPos - 1 - windowSize_left] > maxSupp_left)
			maxSupp_left = coverageVector[curPos - 1 - windowSize_left];
		if (coverageVector[curPos - 1 - windowSize_left] < minSupp_left)
			minSupp_left = coverageVector[curPos - 1 - windowSize_left];
	}

	for (maxSupp_right = minSupp_right = coverageVector[curPos]; windowSize_right < COVERAGE_CHANGE_WINDOW && curPos + windowSize_right < vectorLength; ++windowSize_right)
	{
		curSupp_right += coverageVector[curPos + windowSize_right];
		if (coverageVector[curPos + windowSize_right] > maxSupp_right)
			maxSupp_right = coverageVector[curPos + windowSize_right];
		if (coverageVector[curPos + windowSize_right] < minSupp_right)
			minSupp_right = coverageVector[curPos + windowSize_right];
	}

	if (curSupp_left == 0.0)
		ratio = 2 * COVERAGE_CHANGE_THRESH;
	else
		ratio = (curSupp_right / windowSize_right) / (curSupp_left / windowSize_left);

	if (minSupp_left > maxSupp_right && ratio < 1)
		return ratio;
	else if (maxSupp_left < minSupp_right && ratio > 1)
		return ratio;
	else
		return 1.0;
}


void getExonSupport(int index, long junctionListStart)
{
	//get junction support from fragent file
	ifstream exonSupportFile;
	string inputfilename;
	if (index < 0)
	{
		inputfilename = inputPath + "exon_1.txt";
		index = SUPPORT_VECTOR_SIZE;
	} 
	else
	{

#ifdef COVERAGE
		inputfilename = inputPath + "exon_" + itostr(index+1) + ".txt";
#else
		inputfilename = inputPath + "read_" + itostr(index+1) + ".txt";
#endif 

	}
	exonSupportFile.open(inputfilename.c_str());

	long startPosition, endPosition, juncIndex, tmp, tmpVectorLength, iLoop, jLoop, curFragStart, curSupport, read_retainCheck_start, read_retainCheck_end, fragAddCnt, tmpExonAddLength, altersiteFragStart, altersiteFragEnd, altersiteFragSupport;
	double changeRatio;
	fragment *newFragment;
	vector <fragment*> fragBuffer;
	fragBuffer.reserve(1000000);
	rangeJunction *newRangeJunc;
	deletedSites *newDelSite;
	char chromosome[100];
	int spliced;
	bool makingNewFrag, *read_inExon, doCheck, largeGap;
	startPosition = 0;
	juncIndex = junctionListStart;


	chromosome[0] = '\0';
	exonSupportFile >> chromosome;
	while (chromosome[0] != '\0')
	{
		exonSupportFile >> startPosition;
		exonSupportFile >> endPosition;
		endPosition--;
		exonSupportFile >> spliced;

		if (DORETAIN == true)
		{
			if (endPosition - startPosition + 1 + 1 > 10000)
			{
				cout << "error_" << index << "_" << startPosition << "_" << endPosition << " ";
			}
			read_inExon = new bool [endPosition - startPosition + 1 + 1]; //another 1 for sentinel
			for (iLoop = 0; iLoop <= endPosition - startPosition; iLoop++)
			{
				read_inExon[iLoop] = false;
			}
			read_inExon[endPosition - startPosition + 1] = true;
		}

		while (juncIndex <= sortList_Junction_Num && startPosition > sortList_Junction[juncIndex]->junc->end)
		{
			if (DOTRIMMING == true)// && juncIndex != junctionListStart && juncIndex != sortList_Junction_Num)
			{
				//check & separate
				fragAddCnt = 0;
				fragBuffer.clear();

				if (sortList_Junction[juncIndex]->junc->coverage != NULL)
				{
					curSupport = 0;
					makingNewFrag = false;
					curFragStart = sortList_Junction[juncIndex]->junc->start;
					tmpVectorLength = sortList_Junction[juncIndex]->junc->end - sortList_Junction[juncIndex]->junc->start + 1;

					if (sortList_Junction[juncIndex]->junc->type == frag_exon && sortList_Junction[juncIndex]->junc->altersite > 0)
					{
						changeRatio = cutAlterSite(sortList_Junction[juncIndex]->junc->coverage, tmpVectorLength, sortList_Junction[juncIndex]->junc->start, sortList_Junction[juncIndex]->junc->altersite, curSupport, altersiteFragSupport, altersiteFragStart, altersiteFragEnd);

						if (changeRatio < COVERAGE_CHANGE_THRESH)
						{
							//alternative splice site
							tmpExonAddLength = tmpVectorLength;
							if (tmpExonAddLength >= MIN_ALTER_SPLICE_SITE_LENGTH && curSupport >= coverageThreshold_exon * tmpExonAddLength * SUPPORT_VECTOR_SIZE)
							{
								//get an end, make a new fragment
								newFragment = new fragment;
								newFragment->type = sortList_Junction[juncIndex]->junc->type;
								sprintf(newFragment->name, "ExonAdd_%ld", fragmentID_Cnt + 1);
								newFragment->ID = ++fragmentID_Cnt;

								strcpy(newFragment->chromosome_start, sortList_Junction[juncIndex]->junc->chromosome_start);
								strcpy(newFragment->chromosome_end, sortList_Junction[juncIndex]->junc->chromosome_end);
								newFragment->start = sortList_Junction[juncIndex]->junc->start;
								newFragment->end = sortList_Junction[juncIndex]->junc->end;
								newFragment->support[index] += curSupport;

								++fragAddCnt;
								if (fragAddCnt >= fragBuffer.capacity())
									fragBuffer.reserve(fragBuffer.capacity() + 1000000);
								fragBuffer.push_back(newFragment); 
							}
						}
						else
						{
							//alternative start/end
							tmpExonAddLength = altersiteFragEnd - altersiteFragStart + 1;
							if (tmpExonAddLength >= MIN_EXON_LENGTH && altersiteFragSupport >= coverageThreshold_exon * tmpExonAddLength * SUPPORT_VECTOR_SIZE)
							{
								//get an end, make a new fragment
								newFragment = new fragment;
								newFragment->type = sortList_Junction[juncIndex]->junc->type;
								sprintf(newFragment->name, "ExonAdd_%ld", fragmentID_Cnt + 1);
								newFragment->ID = ++fragmentID_Cnt;

								strcpy(newFragment->chromosome_start, sortList_Junction[juncIndex]->junc->chromosome_start);
								strcpy(newFragment->chromosome_end, sortList_Junction[juncIndex]->junc->chromosome_end);
								newFragment->start = altersiteFragStart;
								newFragment->end = altersiteFragEnd;
								newFragment->support[index] += altersiteFragSupport;

								++fragAddCnt;
								if (fragAddCnt >= fragBuffer.capacity())
									fragBuffer.reserve(fragBuffer.capacity() + 1000000);
								fragBuffer.push_back(newFragment); 
							}
						}
					}
					else
					{

						for (iLoop = 0; iLoop <= tmpVectorLength; iLoop++)
						{
							if (iLoop + MAX_NOCOVERAGE_LENGTH - 1 > tmpVectorLength)
							{
								largeGap = true;
							} 
							else
							{
								largeGap = true;
								for (tmp = 0; tmp < MAX_NOCOVERAGE_LENGTH; tmp++)
								{
									if (sortList_Junction[juncIndex]->junc->coverage[iLoop + tmp] > coverageThreshold_exon)
									{
										largeGap = false;
										break;
									}
								}
							}

							changeRatio = windowChangeRatio(sortList_Junction[juncIndex]->junc->coverage, tmpVectorLength, iLoop);

							if (sortList_Junction[juncIndex]->junc->type == frag_exon)
							{
								if (sortList_Junction[juncIndex]->junc->coverage[iLoop] <= coverageThreshold_exon * SUPPORT_VECTOR_SIZE)
								{
									if (makingNewFrag == true && largeGap == true)
									{
										tmpExonAddLength = sortList_Junction[juncIndex]->junc->start + iLoop - curFragStart;
										if (tmpExonAddLength >= MIN_EXON_LENGTH && curSupport >= coverageThreshold_exon * tmpExonAddLength * SUPPORT_VECTOR_SIZE)
										{
											//get an end, make a new fragment
											newFragment = new fragment;
											newFragment->type = sortList_Junction[juncIndex]->junc->type;
											sprintf(newFragment->name, "ExonAdd_%ld", fragmentID_Cnt + 1);
											newFragment->ID = ++fragmentID_Cnt;

											strcpy(newFragment->chromosome_start, sortList_Junction[juncIndex]->junc->chromosome_start);
											strcpy(newFragment->chromosome_end, sortList_Junction[juncIndex]->junc->chromosome_end);
											newFragment->start = curFragStart;
											if (iLoop < tmpVectorLength)
												newFragment->end = sortList_Junction[juncIndex]->junc->start + iLoop;
											else
												newFragment->end = sortList_Junction[juncIndex]->junc->start + iLoop - 1; //sentinel case
											newFragment->support[index] += curSupport;

											++fragAddCnt;
											if (fragAddCnt >= fragBuffer.capacity())
												fragBuffer.reserve(fragBuffer.capacity() + 1000000);
											fragBuffer.push_back(newFragment); 
										}

										curSupport = 0;
										makingNewFrag = false;								
									}
									else
										curSupport += sortList_Junction[juncIndex]->junc->coverage[iLoop];
								}
								else
								{
									if (makingNewFrag == false)// && changeRatio > COVERAGE_CHANGE_THRESH)
										//if (makingNewFrag == false)
									{
										curFragStart = sortList_Junction[juncIndex]->junc->start + iLoop;
										makingNewFrag = true;
									}
									curSupport += sortList_Junction[juncIndex]->junc->coverage[iLoop];
								}
							}
							else if (sortList_Junction[juncIndex]->junc->type == frag_retained_intron)
							{
								if ((changeRatio < 1/COVERAGE_CHANGE_THRESH || sortList_Junction[juncIndex]->junc->coverage[iLoop] <= coverageThreshold_intron * SUPPORT_VECTOR_SIZE) && largeGap == true)
								{
									if (makingNewFrag == true)
									{
										tmpExonAddLength = sortList_Junction[juncIndex]->junc->start + iLoop - curFragStart;
										if (tmpExonAddLength >= MIN_EXON_LENGTH && curSupport >= coverageThreshold_intron * tmpExonAddLength * SUPPORT_VECTOR_SIZE)
											//if (tmpExonAddLength >= MIN_EXON_LENGTH && tmpExonAddLength < tmpVectorLength - 1 && curSupport >= coverageThreshold_intron * tmpExonAddLength * SUPPORT_VECTOR_SIZE)
										{
											//get an end, make a new fragment
											newFragment = new fragment;
											newFragment->type = sortList_Junction[juncIndex]->junc->type;
											sprintf(newFragment->name, "ExonAdd_%ld", fragmentID_Cnt + 1);
											newFragment->ID = ++fragmentID_Cnt;

											strcpy(newFragment->chromosome_start, sortList_Junction[juncIndex]->junc->chromosome_start);
											strcpy(newFragment->chromosome_end, sortList_Junction[juncIndex]->junc->chromosome_end);
											newFragment->start = curFragStart;
											if (iLoop < tmpVectorLength)
												newFragment->end = sortList_Junction[juncIndex]->junc->start + iLoop;
											else
												newFragment->end = sortList_Junction[juncIndex]->junc->start + iLoop - 1; //sentinel case
											newFragment->support[index] += curSupport;

											++fragAddCnt;
											if (fragAddCnt >= fragBuffer.capacity())
												fragBuffer.reserve(fragBuffer.capacity() + 1000000);
											fragBuffer.push_back(newFragment); 
										}

										curSupport = 0;
										makingNewFrag = false;								
									}
								}
								else
								{
									if (makingNewFrag == false && changeRatio > COVERAGE_CHANGE_THRESH) //for real data
										//if (makingNewFrag == false) //for simulated data
									{
										curFragStart = sortList_Junction[juncIndex]->junc->start + iLoop;
										makingNewFrag = true;
									}
									curSupport += sortList_Junction[juncIndex]->junc->coverage[iLoop];
								}
							}
						}
					}
				} 

				if (fragAddCnt > 1)
				{
					//shift
					sortList_Junction_Num += fragAddCnt - 1;
					for (iLoop = sortList_Junction_Num; iLoop > juncIndex + fragAddCnt - 1; iLoop--)
					{
						sortList_Junction[iLoop] = sortList_Junction[iLoop - fragAddCnt + 1];
					}
				}

				delete sortList_Junction[juncIndex]->junc;
				delete sortList_Junction[juncIndex];

				for (iLoop = 0; iLoop < fragAddCnt; iLoop++)
				{
					newRangeJunc = new rangeJunction;
					newRangeJunc->junc = fragBuffer[iLoop];

					sortList_Junction[juncIndex + iLoop] = newRangeJunc;
				}

				if (fragAddCnt == 0)
				{
					//no expression, the exon is deleted
					//shift back
					sortList_Junction_Num--;
					for (iLoop = juncIndex; iLoop <= sortList_Junction_Num; iLoop++)
					{
						sortList_Junction[iLoop] = sortList_Junction[iLoop + 1];
					}
				}

				juncIndex += fragAddCnt;
			}
			else
			{
				juncIndex++;
			}

		}

		if (juncIndex > sortList_Junction_Num)
		{
			break;
		}

		tmp = 0;
		while (juncIndex + tmp <= sortList_Junction_Num && startPosition <= sortList_Junction[juncIndex + tmp]->junc->end && endPosition >= sortList_Junction[juncIndex + tmp]->junc->start)
		{
			if (DOTRIMMING == true)
			{
				if (sortList_Junction[juncIndex + tmp]->junc->coverage == NULL)
				{
					tmpVectorLength = sortList_Junction[juncIndex + tmp]->junc->end - sortList_Junction[juncIndex + tmp]->junc->start + 1;

					sortList_Junction[juncIndex + tmp]->junc->coverage = new long [tmpVectorLength + 1]; //add one sentinel
					for (iLoop = 0; iLoop <= tmpVectorLength; iLoop++)
						sortList_Junction[juncIndex + tmp]->junc->coverage[iLoop] = 0;
				}
			}

			if (startPosition >= sortList_Junction[juncIndex + tmp]->junc->start)
			{
				if (endPosition <= sortList_Junction[juncIndex + tmp]->junc->end)
				{
					((sortList_Junction[juncIndex + tmp]->junc->support)[index]) += 1 * (endPosition - startPosition);
					if (DOTRIMMING == true)
					{
						for (iLoop = startPosition - sortList_Junction[juncIndex + tmp]->junc->start; iLoop <= endPosition - sortList_Junction[juncIndex + tmp]->junc->start; iLoop++)
							sortList_Junction[juncIndex + tmp]->junc->coverage[iLoop] += 1;
					}
					if (DORETAIN == true)
					{
						for (iLoop = 0; iLoop <= endPosition - startPosition; iLoop++)
							read_inExon[iLoop] = true;
					}
				}
				else if (endPosition > sortList_Junction[juncIndex + tmp]->junc->end)
				{
					((sortList_Junction[juncIndex + tmp]->junc->support)[index]) += 1 * (sortList_Junction[juncIndex + tmp]->junc->end - startPosition + 1);
					if (DOTRIMMING == true)
					{
						for (iLoop = startPosition - sortList_Junction[juncIndex + tmp]->junc->start; iLoop <= sortList_Junction[juncIndex + tmp]->junc->end - sortList_Junction[juncIndex + tmp]->junc->start; iLoop++)
							sortList_Junction[juncIndex + tmp]->junc->coverage[iLoop] += 1;	
					}
					if (DORETAIN == true)
					{
						for (iLoop = 0; iLoop <= sortList_Junction[juncIndex + tmp]->junc->end - startPosition; iLoop++)
							read_inExon[iLoop] = true;
					}
				}
			} 
			else
			{
				if (endPosition <= sortList_Junction[juncIndex + tmp]->junc->end)
				{
					((sortList_Junction[juncIndex + tmp]->junc->support)[index]) += 1 * (endPosition - sortList_Junction[juncIndex + tmp]->junc->start);
					if (DOTRIMMING == true)
					{
						for (iLoop = 0; iLoop <= endPosition - sortList_Junction[juncIndex + tmp]->junc->start; iLoop++)
							sortList_Junction[juncIndex + tmp]->junc->coverage[iLoop] += 1;
					}
					if (DORETAIN == true)
					{
						for (iLoop = sortList_Junction[juncIndex + tmp]->junc->start - startPosition; iLoop <= endPosition - startPosition; iLoop++)
							read_inExon[iLoop] = true;
					}
				}
				else if (endPosition > sortList_Junction[juncIndex + tmp]->junc->end)
				{
					((sortList_Junction[juncIndex + tmp]->junc->support)[index]) += 1 * (sortList_Junction[juncIndex + tmp]->junc->end - sortList_Junction[juncIndex + tmp]->junc->start + 1);
					if (DOTRIMMING == true)
					{
						for (iLoop = 0; iLoop <= sortList_Junction[juncIndex + tmp]->junc->end - sortList_Junction[juncIndex + tmp]->junc->start; iLoop++)
							sortList_Junction[juncIndex + tmp]->junc->coverage[iLoop] += 1;	
					}
					if (DORETAIN == true)
					{
						for (iLoop = sortList_Junction[juncIndex + tmp]->junc->start - startPosition; iLoop <= sortList_Junction[juncIndex + tmp]->junc->end - startPosition; iLoop++)
							read_inExon[iLoop] = true;
					}
				}
			}
			tmp++;
		}

		if (DORETAIN == true)
		{
			doCheck = false;
			tmp = 0;
			for (iLoop = 0; iLoop <= endPosition - startPosition + 1; iLoop++)
			{
				if (read_inExon[iLoop] == true && doCheck == true)
				{
					//do check
					read_retainCheck_end = iLoop - 1;

					while (juncIndex + tmp <= sortList_Junction_Num)
					{
						if ((juncIndex + tmp - 1 < junctionListStart || juncIndex + tmp - 1 >= junctionListStart && startPosition + read_retainCheck_start > sortList_Junction[juncIndex + tmp - 1]->junc->end) && startPosition + read_retainCheck_end < sortList_Junction[juncIndex + tmp]->junc->start)
						{
							break;
						}
						tmp++;
					}


					//retained intron

					sortList_Junction_Num++;
					for (jLoop = sortList_Junction_Num; jLoop > juncIndex + tmp; jLoop--)
					{
						sortList_Junction[jLoop] = sortList_Junction[jLoop - 1];
					}

					//make an exon
					newFragment = new fragment;
					newFragment->type = frag_retained_intron;
					sprintf(newFragment->name, "RetainedIntron_&ld", juncIndex);
					newFragment->ID = ++fragmentID_Cnt;

					strcpy(newFragment->chromosome_start, sortList_Junction[juncIndex + tmp -1]->junc->chromosome_start);
					strcpy(newFragment->chromosome_end, sortList_Junction[juncIndex + tmp -1]->junc->chromosome_end);
					if (juncIndex + tmp - 1 >= junctionListStart)
					{
						newFragment->start = sortList_Junction[juncIndex + tmp -1]->junc->end + 1;
						//newFragment->start = sortList_Junction[juncIndex + tmp -1]->junc->end;// + 1;
					} 
					else
					{
						newFragment->start = startPosition;
					}
					if (juncIndex + tmp + 1 <= sortList_Junction_Num)
					{
						newFragment->end = sortList_Junction[juncIndex + tmp +1]->junc->start - 1;
						//newFragment->end = sortList_Junction[juncIndex + tmp +1]->junc->start;// - 1;
					} 
					else
					{
						newFragment->end = endPosition;
					}
					newFragment->support[index] += 1 * (read_retainCheck_end - read_retainCheck_start + 1);

					tmpVectorLength = newFragment->end - newFragment->start + 1;
					newFragment->coverage = new long [tmpVectorLength + 1]; //add one sentinel

					for (jLoop = 0; jLoop <= tmpVectorLength; jLoop++)
						newFragment->coverage[jLoop] = 0;
					for (jLoop = startPosition + read_retainCheck_start - newFragment->start; jLoop <= startPosition + read_retainCheck_end - newFragment->start; jLoop++)
						newFragment->coverage[jLoop] += 1;


					//insert exonic fragments to list
					newRangeJunc = new rangeJunction;
					newRangeJunc->junc = newFragment;

					sortList_Junction[juncIndex + tmp] = newRangeJunc;


					doCheck = false;
				}
				else if (read_inExon[iLoop] == false && doCheck == false)
				{
					read_retainCheck_start = iLoop;
					doCheck = true;
				}
			}
		}


		if (DORETAIN == true)
		{
			delete [] read_inExon;
		}

		chromosome[0] = '\0';
		exonSupportFile >> chromosome;
	}

	exonSupportFile.close();

	while (juncIndex <= sortList_Junction_Num)
	{
		if (DOTRIMMING == true)// && juncIndex != junctionListStart && juncIndex != sortList_Junction_Num)
		{
			//check & separate
			fragAddCnt = 0;
			fragBuffer.clear();

			if (sortList_Junction[juncIndex]->junc->coverage != NULL)
			{
				curSupport = 0;
				makingNewFrag = false;
				curFragStart = sortList_Junction[juncIndex]->junc->start;
				tmpVectorLength = sortList_Junction[juncIndex]->junc->end - sortList_Junction[juncIndex]->junc->start + 1;

				if (sortList_Junction[juncIndex]->junc->type == frag_exon && sortList_Junction[juncIndex]->junc->altersite > 0)
				{
					changeRatio = cutAlterSite(sortList_Junction[juncIndex]->junc->coverage, tmpVectorLength, sortList_Junction[juncIndex]->junc->start, sortList_Junction[juncIndex]->junc->altersite, curSupport, altersiteFragSupport, altersiteFragStart, altersiteFragEnd);

					if (changeRatio < COVERAGE_CHANGE_THRESH)
					{
						//alternative splice site
						tmpExonAddLength = tmpVectorLength;
						if (tmpExonAddLength >= MIN_ALTER_SPLICE_SITE_LENGTH && curSupport >= coverageThreshold_exon * tmpExonAddLength * SUPPORT_VECTOR_SIZE)
						{
							//get an end, make a new fragment
							newFragment = new fragment;
							newFragment->type = sortList_Junction[juncIndex]->junc->type;
							sprintf(newFragment->name, "ExonAdd_%ld", fragmentID_Cnt + 1);
							newFragment->ID = ++fragmentID_Cnt;

							strcpy(newFragment->chromosome_start, sortList_Junction[juncIndex]->junc->chromosome_start);
							strcpy(newFragment->chromosome_end, sortList_Junction[juncIndex]->junc->chromosome_end);
							newFragment->start = sortList_Junction[juncIndex]->junc->start;
							newFragment->end = sortList_Junction[juncIndex]->junc->end;
							newFragment->support[index] += curSupport;

							++fragAddCnt;
							if (fragAddCnt >= fragBuffer.capacity())
								fragBuffer.reserve(fragBuffer.capacity() + 1000000);
							fragBuffer.push_back(newFragment); 
						}
					}
					else
					{
						//alternative start/end
						tmpExonAddLength = altersiteFragEnd - altersiteFragStart + 1;
						if (tmpExonAddLength >= MIN_EXON_LENGTH && altersiteFragSupport >= coverageThreshold_exon * tmpExonAddLength * SUPPORT_VECTOR_SIZE)
						{
							//get an end, make a new fragment
							newFragment = new fragment;
							newFragment->type = sortList_Junction[juncIndex]->junc->type;
							sprintf(newFragment->name, "ExonAdd_%ld", fragmentID_Cnt + 1);
							newFragment->ID = ++fragmentID_Cnt;

							strcpy(newFragment->chromosome_start, sortList_Junction[juncIndex]->junc->chromosome_start);
							strcpy(newFragment->chromosome_end, sortList_Junction[juncIndex]->junc->chromosome_end);
							newFragment->start = altersiteFragStart;
							newFragment->end = altersiteFragEnd;
							newFragment->support[index] += altersiteFragSupport;

							++fragAddCnt;
							if (fragAddCnt >= fragBuffer.capacity())
								fragBuffer.reserve(fragBuffer.capacity() + 1000000);
							fragBuffer.push_back(newFragment); 
						}
					}
				}
				else
				{

					for (iLoop = 0; iLoop <= tmpVectorLength; iLoop++)
					{
						if (iLoop + MAX_NOCOVERAGE_LENGTH - 1 > tmpVectorLength)
						{
							largeGap = true;
						} 
						else
						{
							largeGap = true;
							for (tmp = 0; tmp < MAX_NOCOVERAGE_LENGTH; tmp++)
							{
								if (sortList_Junction[juncIndex]->junc->coverage[iLoop + tmp] > coverageThreshold_exon)
								{
									largeGap = false;
									break;
								}
							}
						}

						changeRatio = windowChangeRatio(sortList_Junction[juncIndex]->junc->coverage, tmpVectorLength, iLoop);

						if (sortList_Junction[juncIndex]->junc->type == frag_exon)
						{
							if (sortList_Junction[juncIndex]->junc->coverage[iLoop] <= coverageThreshold_exon * SUPPORT_VECTOR_SIZE)
							{
								if (makingNewFrag == true && largeGap == true)
								{
									tmpExonAddLength = sortList_Junction[juncIndex]->junc->start + iLoop - curFragStart;
									if (tmpExonAddLength >= MIN_EXON_LENGTH && curSupport >= coverageThreshold_exon * tmpExonAddLength * SUPPORT_VECTOR_SIZE)
									{
										//get an end, make a new fragment
										newFragment = new fragment;
										newFragment->type = sortList_Junction[juncIndex]->junc->type;
										sprintf(newFragment->name, "ExonAdd_%ld", fragmentID_Cnt + 1);
										newFragment->ID = ++fragmentID_Cnt;

										strcpy(newFragment->chromosome_start, sortList_Junction[juncIndex]->junc->chromosome_start);
										strcpy(newFragment->chromosome_end, sortList_Junction[juncIndex]->junc->chromosome_end);
										newFragment->start = curFragStart;
										if (iLoop < tmpVectorLength)
											newFragment->end = sortList_Junction[juncIndex]->junc->start + iLoop;
										else
											newFragment->end = sortList_Junction[juncIndex]->junc->start + iLoop - 1; //sentinel case
										newFragment->support[index] += curSupport;

										++fragAddCnt;
										if (fragAddCnt >= fragBuffer.capacity())
											fragBuffer.reserve(fragBuffer.capacity() + 1000000);
										fragBuffer.push_back(newFragment); 
									}

									curSupport = 0;
									makingNewFrag = false;								
								}
								else
									curSupport += sortList_Junction[juncIndex]->junc->coverage[iLoop];
							}
							else
							{
								if (makingNewFrag == false)// && changeRatio > COVERAGE_CHANGE_THRESH)
									//if (makingNewFrag == false)
								{
									curFragStart = sortList_Junction[juncIndex]->junc->start + iLoop;
									makingNewFrag = true;
								}
								curSupport += sortList_Junction[juncIndex]->junc->coverage[iLoop];
							}
						}
						else if (sortList_Junction[juncIndex]->junc->type == frag_retained_intron)
						{
							if ((changeRatio < 1/COVERAGE_CHANGE_THRESH || sortList_Junction[juncIndex]->junc->coverage[iLoop] <= coverageThreshold_intron * SUPPORT_VECTOR_SIZE) && largeGap == true)
							{
								if (makingNewFrag == true)
								{
									tmpExonAddLength = sortList_Junction[juncIndex]->junc->start + iLoop - curFragStart;
									if (tmpExonAddLength >= MIN_EXON_LENGTH && curSupport >= coverageThreshold_intron * tmpExonAddLength * SUPPORT_VECTOR_SIZE)
										//if (tmpExonAddLength >= MIN_EXON_LENGTH && tmpExonAddLength < tmpVectorLength - 1 && curSupport >= coverageThreshold_intron * tmpExonAddLength * SUPPORT_VECTOR_SIZE)
									{
										//get an end, make a new fragment
										newFragment = new fragment;
										newFragment->type = sortList_Junction[juncIndex]->junc->type;
										sprintf(newFragment->name, "ExonAdd_%ld", fragmentID_Cnt + 1);
										newFragment->ID = ++fragmentID_Cnt;

										strcpy(newFragment->chromosome_start, sortList_Junction[juncIndex]->junc->chromosome_start);
										strcpy(newFragment->chromosome_end, sortList_Junction[juncIndex]->junc->chromosome_end);
										newFragment->start = curFragStart;
										if (iLoop < tmpVectorLength)
											newFragment->end = sortList_Junction[juncIndex]->junc->start + iLoop;
										else
											newFragment->end = sortList_Junction[juncIndex]->junc->start + iLoop - 1; //sentinel case
										newFragment->support[index] += curSupport;

										++fragAddCnt;
										if (fragAddCnt >= fragBuffer.capacity())
											fragBuffer.reserve(fragBuffer.capacity() + 1000000);
										fragBuffer.push_back(newFragment); 
									}

									curSupport = 0;
									makingNewFrag = false;								
								}
							}
							else
							{
								if (makingNewFrag == false && changeRatio > COVERAGE_CHANGE_THRESH)
									//if (makingNewFrag == false)
								{
									curFragStart = sortList_Junction[juncIndex]->junc->start + iLoop;
									makingNewFrag = true;
								}
								curSupport += sortList_Junction[juncIndex]->junc->coverage[iLoop];
							}
						}
					}
				}
			} 


			if (fragAddCnt > 1)
			{
				//shift
				sortList_Junction_Num += fragAddCnt - 1;
				for (iLoop = sortList_Junction_Num; iLoop > juncIndex + fragAddCnt - 1; iLoop--)
				{
					sortList_Junction[iLoop] = sortList_Junction[iLoop - fragAddCnt + 1];
				}
			}

			delete sortList_Junction[juncIndex]->junc;
			delete sortList_Junction[juncIndex];

			for (iLoop = 0; iLoop < fragAddCnt; iLoop++)
			{
				newRangeJunc = new rangeJunction;
				newRangeJunc->junc = fragBuffer[iLoop];

				sortList_Junction[juncIndex + iLoop] = newRangeJunc;
			}

			if (fragAddCnt == 0)
			{
				//no expression, the exon is deleted
				//shift back
				sortList_Junction_Num--;
				for (iLoop = juncIndex; iLoop <= sortList_Junction_Num; iLoop++)
				{
					sortList_Junction[iLoop] = sortList_Junction[iLoop + 1];
				}
			}

			juncIndex += fragAddCnt;
		}
		else
		{
			juncIndex++;
		}

	}

	fragBuffer.clear();
	free_vector(fragBuffer);

	return;
}

void inputData()
{
	string filename;
	int sampleLoopCnt;
	long prevJuncNum;
	fragment *curFrag;

	/************************************************************************/
	/* INPUT ALL JUNCTIONS                                                  */
	/************************************************************************/
	filename = inputPath + "junction_1.txt";
	input_junction(filename);

#ifdef COUTSTEPS
	cout << "input junction done." << endl;
#endif


	/************************************************************************/
	/* INPUT ALL JUNCTIONS                                                  */
	/************************************************************************/

#ifndef JUNCTIONONLY
	for (sampleLoopCnt = 0; sampleLoopCnt < SUPPORT_VECTOR_SIZE; sampleLoopCnt++)
	{
		getJunctionSupport(sampleLoopCnt);
	}
	cleanSpliceDirection();
#endif

#ifdef COUTSTEPS
	cout << "junction support done." << endl;
	cout << "junc " << sortList_Junction_Num << "\t" << flush;
#endif

	/************************************************************************/
	/* INPUT ALL JUNCTIONS                                                  */
	/************************************************************************/

	for (long i = 1; i <= sortList_Junction_Num; i++)
		if (sortList_Junction[i]->junc->type == frag_junction)
		{
			curFrag = sortList_Junction[i]->junc;
			splice_all_file << curFrag->chromosome_start << "\t" << curFrag->start - 19 << "\t" << curFrag->end + 20 << "\talljunction\t1000\t";
			if (curFrag->transDirection == undetermined || curFrag->transDirection == sense)
				splice_all_file << "+\t";
			else
				splice_all_file << "-\t";
			splice_all_file << curFrag->start - 19 << "\t" << curFrag->end + 20 << "\t0\t2\t" << "20,20,\t0," << curFrag->end - curFrag->start + 19 << endl;
		}


#ifdef FILTER_JUNCTION
		filterJunction(MAX_JUNSUPPORT, 16, MEAN_JUNSUPPORT); //throw a junction away if zeroCnt >= 2nd parameter or max support <= 1st parameter
#ifdef COUTSTEPS
		cout << "filtering... " << sortList_Junction_Num << "\t" << flush;
#endif
#endif

		for (long i = 1; i <= sortList_Junction_Num; i++)
			if (sortList_Junction[i]->junc->type == frag_junction)
			{
				curFrag = sortList_Junction[i]->junc;
				splice_filtered_file << curFrag->chromosome_start << "\t" << curFrag->start - 19 << "\t" << curFrag->end + 20 << "\tfilteredjunction\t1000\t";
				if (curFrag->transDirection == undetermined || curFrag->transDirection == sense)
					splice_filtered_file << "+\t";
				else
					splice_filtered_file << "-\t";
				splice_filtered_file << curFrag->start - 19 << "\t" << curFrag->end + 20 << "\t0\t2\t" << "20,20,\t0," << curFrag->end - curFrag->start + 19 << endl;
			}

			prevJuncNum = sortList_Junction_Num;

			getExons();
#ifdef COUTSTEPS
			cout << "get exons done." << endl;
#endif
#ifndef JUNCTIONONLY
			getExonSupport(-1, prevJuncNum + 1);
			DOTRIMMING = false;
			DORETAIN = false;
			for (sampleLoopCnt = 0; sampleLoopCnt < SUPPORT_VECTOR_SIZE; sampleLoopCnt++)
			{
				getExonSupport(sampleLoopCnt, prevJuncNum + 1);
			}
#endif

			/************************************************************************/
			/* CALCULATE COVERAGE                                                   */
			/************************************************************************/
#ifdef COVERAGE
			//calculate read coverage
			for (long i = 1; i <= sortList_Junction_Num; i++)
			{
				if (sortList_Junction[i]->junc->type == frag_exon || sortList_Junction[i]->junc->type == frag_retained_intron)
				{
					for (sampleLoopCnt = 0; sampleLoopCnt < SUPPORT_VECTOR_SIZE; sampleLoopCnt++)
					{
						sortList_Junction[i]->junc->support[sampleLoopCnt] = sortList_Junction[i]->junc->support[sampleLoopCnt] / (sortList_Junction[i]->junc->end - sortList_Junction[i]->junc->start + 1);
					}
				}
			}
#endif

#ifdef NORMALIZE_COVERAGE
			//normalize exon coverage and junction support
			for (sampleLoopCnt = 0; sampleLoopCnt < SUPPORT_VECTOR_SIZE; sampleLoopCnt++)
			{
				if (DatasetReadCount_total[sampleLoopCnt] > 0)
				{
					DatasetNormalizationRatio[sampleLoopCnt] = double(DatasetReadCount_normalization) / DatasetReadCount_total[sampleLoopCnt];
					for (long i = 1; i <= sortList_Junction_Num; i++)
					{
						sortList_Junction[i]->junc->support[sampleLoopCnt] = sortList_Junction[i]->junc->support[sampleLoopCnt] * DatasetNormalizationRatio[sampleLoopCnt];
					}
				}
			}
#endif
#ifdef COUTSTEPS
			cout << "total " << sortList_Junction_Num << "\t" << flush;
#endif
#ifdef FILTER_FRAGMENTS
			filterIntron(prevJuncNum + 1); //filter introns
#ifdef COUTSTEPS
			cout << "filtering intron... " << sortList_Junction_Num << "\t" << flush;
#endif
#endif

#ifdef COUTSTEPS
			cout << "get junctions and exons done." << endl;
#endif

			return;
}

RangeJunctionList* buildOrigList()
{
	RangeJunctionList *resultList;
	alter_junction *curAlterJunc;
	long i;
	int alterJuncCnt;

	//#ifdef FILTER_FRAGMENTS
	filterFalseFragments(false);
	//#endif

#ifdef DO_MERGEALTERSITES
	mergeAlterSpliceJunc();

	//output alternative splice sites to leaf file
	for (i = 1; i <= sortList_Junction_Num; i++)
	{
		if (sortList_Junction[i]->junc->alter != NULL)
		{
			alterJuncCnt = 0;
			sortList_Junction[i]->junc->alterFragCnt = 0;
			curAlterJunc = sortList_Junction[i]->junc->alter;
			while (curAlterJunc != NULL)
			{
				(sortList_Junction[i]->junc->alterFragCnt)++;
				if (curAlterJunc->juncInfo->type == frag_junction)
				{
					alterJuncCnt++;
				}
				curAlterJunc = curAlterJunc->next;
			}
			// if (alterJuncCnt > 1000)
			// 	cout << "weired" << endl;
			leafsizefile << "alter\t" << sortList_Junction[i]->junc->start << '\t' << sortList_Junction[i]->junc->end << '\t' << alterJuncCnt << endl;
		}
	}
#else
	for (i = 1; i <= sortList_Junction_Num; i++)
	{
		sortKey_Junction[i] = sortList_Junction[i]->junc->end;
	}
	mergeSort_JunctionSort(sortList_Junction_Num);

	for (i = 1; i <= sortList_Junction_Num; i++)
	{
		sortKey_Junction[i] = sortList_Junction[i]->junc->start;
	}
	mergeSort_JunctionSort(sortList_Junction_Num);
#endif

	resultList = new RangeJunctionList;
	resultList->rangeLow = CHROMOSOME_START;
	resultList->rangeHigh = CHROMOSOME_END;

	for (i = sortList_Junction_Num; i >= 1; i--)
	{
		sortList_Junction[i]->next = resultList->list;
		resultList->list = sortList_Junction[i];
	}	

	return resultList;
}


RangeJunctionList* separateGene(RangeJunctionList* origList, bool &separable)
{
	//separate genes
	//return a set of lists, each of which corresponds to a gene
	//IMPORTANT: RangeJunctionList is head-inserting, for the convenience of stack

	if (origList == NULL)
	{
		separable = false;
		return NULL;
	}

	RangeJunctionList *resultList, *curList;
	resultList = NULL;
	curList = NULL;

	rangeJunction *curJunc, *curListTail;
	curJunc = NULL;
	curListTail = NULL;

	long endBoard = -2;
	int numIndep = 0; //number of independent regions

	curJunc = origList->list;
	while (curJunc != NULL)
	{
		if (curJunc->junc->start > endBoard+1)
		{
			//get an independent region
			numIndep++;

			//omit the first list, it is empty
			if (curList == NULL)
			{
				//do nothing
			} 
			else
			{
				curList->rangeHigh = endBoard;
			}

			curList = new RangeJunctionList;
			curList->nextList = resultList;
			curList->rangeLow = curJunc->junc->start;

			resultList = curList;

			curListTail = NULL;
		}

		//process current fragment
		if (curListTail == NULL)
		{
			curList->list = curJunc;
			curListTail = curJunc;
		}
		else
		{
			curListTail->next = curJunc;
			curListTail = curJunc;
		}

		if (curJunc->junc->end > endBoard)
		{
			endBoard = curJunc->junc->end;
		}

		curJunc = curJunc->next;
		curListTail->next = NULL;
	}

	if (curList == NULL)
	{
		cout << "curList == NULL" << endl;
		exit(1);
	} 
	else
	{
		//curList->rangeHigh = origList->rangeHigh;
		curList->rangeHigh = endBoard;
	}

	if (numIndep == 1)
		separable = false; //cannot be separated into multiple regions
	else if (numIndep > 1)
		separable = true; //can be separated
	else
	{
		cout << "numGene == 0" << endl;
		exit(1);
	}

	origList->list = NULL;
	delete origList;

	return resultList;
}


// decide whether two adjacent fragments are compatible
bool compatibleJunctions(rangeJunction* junctionA, rangeJunction* junctionB)
{
	bool compatible;
	if (junctionA->junc->end <= junctionB->junc->start && junctionB->junc->start - junctionA->junc->end <= 1)
	{
		compatible = true;
		if (junctionB->junc->start - junctionA->junc->end == 1)
		{
			if (junctionA->junc->type == frag_junction && (junctionB->junc->type == frag_exon || junctionB->junc->type == frag_retained_intron) || junctionB->junc->type == frag_junction && (junctionA->junc->type == frag_exon || junctionA->junc->type == frag_retained_intron))
			{
				compatible = false;
			}
		}
	}
	else if (junctionA->junc->start >= junctionB->junc->end && junctionA->junc->start - junctionB->junc->end <= 1)
	{
		compatible = true;
		if (junctionA->junc->start - junctionB->junc->end == 1)
		{
			if (junctionA->junc->type == frag_junction && (junctionB->junc->type == frag_exon || junctionB->junc->type == frag_retained_intron) || junctionB->junc->type == frag_junction && (junctionA->junc->type == frag_exon || junctionA->junc->type == frag_retained_intron))
			{
				compatible = false;
			}
		}
	}
	else
		compatible = false;

	return compatible;
}


bool compatibleStrand(rangeJunction *junctionA, rangeJunction *junctionB)
{
	bool compatible = true;
	if ((junctionA->junc->transDirection == sense && junctionB->junc->transDirection == antisense) || (junctionA->junc->transDirection == antisense && junctionB->junc->transDirection == sense))
	{
		compatible = false;
	}
	return compatible;
}


void separateJuncGraph_directed(RangeJunctionList* origList)
{
	//construct fragment graph based on given junctions

	JuncGraphVertex *curVertex, *tailVertex, *edgeVertex, *newVertex, *preVertex;
	JuncGraphEdge *newEdge;
	rangeJunction *curJunc;

	//build vertices

	curJunc = origList->list;
	while (curJunc != NULL)
	{
		if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
		{
			curVertex = new JuncGraphVertex;
			curVertex->corresJunc = curJunc;
			curVertex->corresJunc->junc->transDirection = sense;
			newVertex = new JuncGraphVertex;
			newVertex->corresJunc = new rangeJunction;
			newVertex->corresJunc->junc = curJunc->junc->clone();
			newVertex->corresJunc->junc->transDirection = antisense;
			curJunc = curJunc->next;
			curVertex->corresJunc->next = NULL;
			newVertex->corresJunc->next = NULL;
			if (junctionGraph->vertices == NULL)
			{
				junctionGraph->vertices = curVertex;
				curVertex->next = newVertex;
				tailVertex = newVertex;
			} 
			else
			{
				tailVertex->next = curVertex;
				curVertex->next = newVertex;
				tailVertex = newVertex;
			}

		}
		else
		{
			curVertex = new JuncGraphVertex;
			curVertex->corresJunc = curJunc;
			newVertex = new JuncGraphVertex;
			curJunc = curJunc->next;
			curVertex->corresJunc->next = NULL;

			if (junctionGraph->vertices == NULL)
			{
				junctionGraph->vertices = curVertex;
				tailVertex = curVertex;
			} 
			else
			{
				tailVertex->next = curVertex;
				tailVertex = curVertex;
			}
		}
	}

	//build edges

	curVertex = junctionGraph->vertices;
	while (curVertex != NULL)
	{
		edgeVertex = curVertex->next;
		while (edgeVertex != NULL)
		{
			if (compatibleJunctions((curVertex->corresJunc), (edgeVertex->corresJunc)) == true && compatibleStrand((curVertex->corresJunc), (edgeVertex->corresJunc)) == true)
			{
				//create an edge for curVertex
				newEdge = new JuncGraphEdge;
				newEdge->linkedVertex = edgeVertex;

				newEdge->next = curVertex->edges;
				curVertex->edges = newEdge;

				edgeVertex->hasInEdge = true;
				curVertex->hasOutEdge = true;
			} 
			else
			{
				//break; 
			}
			edgeVertex = edgeVertex->next;
		}

		curVertex = curVertex->next;
	}

	preVertex = NULL;
	curVertex = junctionGraph->vertices;
	while (curVertex != NULL)
	{
		if ((curVertex->corresJunc->junc->type == frag_exon || curVertex->corresJunc->junc->type == frag_retained_intron) && curVertex->hasInEdge == false && curVertex->hasOutEdge == false)
		{
			// delete this vertex
			if (preVertex == NULL)
			{
				junctionGraph->vertices = curVertex->next;
				delete curVertex;
				curVertex = junctionGraph->vertices;
			}
			else
			{
				preVertex->next = curVertex->next;
				delete curVertex;
				curVertex = preVertex->next;
			}
		}
		else
		{
			preVertex = curVertex;
			curVertex = curVertex->next;
		}
	}

	return;
}

void DFS_visit(JuncGraphVertex *u)
{
	//visit a vertex u during DFS
	u->traversed = true;

	JuncGraphVertex *v;
	JuncGraphEdge *curEdge;
	curEdge = u->edges;
	while (curEdge != NULL)
	{
		v = curEdge->linkedVertex;
		if (v->traversed == false)
		{
			DFS_visit(v);
		}
		curEdge = curEdge->next;
	}

	//add u into sortList_Junction
	sortList_Junction_Num++;
	sortList_Junction[sortList_Junction_Num] = u->corresJunc;
	u->corresJunc = NULL;

	return;
}

RangeJunctionList* search_conn_comp(int &numCC, long rangeLow, long rangeHigh)
{
	//search all connected component within junctionGraph
	JuncGraphVertex *curVertex;
	RangeJunctionList *resultList, *curList;
	long i;

	numCC = 0; //number of connected component

	resultList = NULL;
	sortList_Junction_Num = 0;
	curVertex = junctionGraph->vertices;
	while (curVertex != NULL)
	{
		if (curVertex->traversed == false)
		{
			DFS_visit(curVertex);

			numCC++;

			curList = new RangeJunctionList;

			//get a connected component
			for (i = 1; i <= sortList_Junction_Num; i++)
			{
				sortKey_Junction[i] = sortList_Junction[i]->junc->end;
			}
			mergeSort_JunctionSort(sortList_Junction_Num);

			curList->rangeHigh = sortList_Junction[sortList_Junction_Num]->junc->end;

			for (i = 1; i <= sortList_Junction_Num; i++)
			{
				sortKey_Junction[i] = sortList_Junction[i]->junc->start;
			}
			mergeSort_JunctionSort(sortList_Junction_Num);

			curList->rangeLow = sortList_Junction[1]->junc->start;

			curList->nextList = resultList;
			resultList = curList;

			curList->listtail = sortList_Junction[sortList_Junction_Num];
			for (i = sortList_Junction_Num; i >= 1; i--)
			{
				sortList_Junction[i]->next = curList->list;
				curList->list = sortList_Junction[i];
			}			

			sortList_Junction_Num = 0;
		}

		curVertex = curVertex->next;
	}

	return resultList;
}

bool ValidFragList(RangeJunctionList *list)
{
	bool valid, noJunc;
	long fragNum;
	rangeJunction *curfrag;

	noJunc = true;
	fragNum = 0;
	curfrag = list->list;
	while(curfrag != NULL)
	{
		if (curfrag->junc->type == frag_junction)
		{
			noJunc = false;
		}
		else
		{
			fragNum++;
		}
		curfrag = curfrag->next;
	}

	valid = true;
	if (noJunc == true && fragNum > 0)
	{
		valid = false;
	}
	return valid;
}


RangeJunctionList* SeparateOverlapGene(RangeJunctionList *geneList, bool &separable_gene)
{
	// pre process the fragment list
	separable_gene = false;
	trans_direction direction = undetermined;
	rangeJunction *curJunc;

	curJunc = geneList->list;
	while (curJunc != NULL)
	{
		if (curJunc->junc->type == frag_junction)
		{
			if (direction == undetermined)
			{
				direction = curJunc->junc->transDirection;
			}
			else
			{
				if (direction != curJunc->junc->transDirection)
				{
					separable_gene = true;
					break;
				}
			}
		}
		curJunc = curJunc->next;
	}

	RangeJunctionList *resultList;
	if (separable_gene == true)
	{
		junctionGraph = new JuncGraph;
		RangeJunctionList *tempList, *prevList, *curList;
		tempList = geneList->clone();
		separateJuncGraph_directed(tempList);

		int numCC = 0;
		resultList = search_conn_comp(numCC, geneList->rangeLow, geneList->rangeHigh);
		bool valid = true;
		prevList = NULL;
		curList = resultList;
		while(curList != NULL)
		{
			valid = ValidFragList(curList);
			if (valid == false)
			{
				if (prevList == NULL)
				{
					resultList = curList->nextList;
					delete curList;
					curList = resultList;
				}
				else
				{
					prevList->nextList = curList->nextList;
					delete curList;
					curList = prevList->nextList;
				}
			}
			else
			{
				prevList = curList;
				curList = curList->nextList;
			}

		}

		delete junctionGraph;
	}
	else
	{
		resultList = geneList->clone();
	}

	return resultList;
}

bool calculateGeneMeanCoverage(RangeJunctionList *fragmentList, double meanCoverageList[SUPPORT_VECTOR_SIZE])
{
	rangeJunction *curJunc;

	//count number of exons in the gene
	long exonCnt = 0, exonCntLoop;
	int sampleLoop;

	curJunc = fragmentList->list;
	while (curJunc != NULL)
	{
		if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
			++exonCnt;
		curJunc = curJunc->next;
	}

	if (exonCnt < 5)
		return false; //no need to do filtering

	//get expression array
	double **coverageArray = new double* [SUPPORT_VECTOR_SIZE];
	for (sampleLoop = 0; sampleLoop < SUPPORT_VECTOR_SIZE; ++sampleLoop)
	{
		coverageArray[sampleLoop] = new double [exonCnt + 1];

		curJunc = fragmentList->list;
		exonCntLoop = 0;
		while (curJunc != NULL)
		{
			if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
				coverageArray[sampleLoop][++exonCntLoop] = curJunc->junc->support[sampleLoop];
			curJunc = curJunc->next;
		}
	}

	//calculate mean coverage in each sample
	long index_q25, index_q75;
	double total_coverage;
	index_q25 = long(floor(0.25 * (exonCnt - 1) + 1));
	index_q75 = long(floor(0.75 * (exonCnt - 1) + 1));

	for (sampleLoop = 0; sampleLoop < SUPPORT_VECTOR_SIZE; ++sampleLoop)
	{
		quicksort(coverageArray[sampleLoop], exonCnt);

		total_coverage = 0.0;
		for (exonCntLoop = index_q25 + 1; exonCntLoop <= index_q75; ++exonCntLoop)
			total_coverage += coverageArray[sampleLoop][exonCntLoop];

		meanCoverageList[sampleLoop] = total_coverage / (index_q75 - index_q25);
	}

	//clean up
	for (sampleLoop = 0; sampleLoop < SUPPORT_VECTOR_SIZE; ++sampleLoop)
	{
		delete [] coverageArray[sampleLoop];
	}
	delete[] coverageArray;

	return true;
}

void filterGeneLowCoverageExon(GTvertex *targetVertex)
{
	double meanCoverageList[SUPPORT_VECTOR_SIZE];
	rangeJunction *curJunc, *prevJunc;
	long iLoop;
	int thresh_lowCovSampleCnt = int(ceil(SUPPORT_VECTOR_SIZE * (1 - coverageThreshold_GeneExon))), lowCovSampleCnt, sampleLoop;

	//find mean coverage of 25% to 75% percentile expression in this gene
	if (calculateGeneMeanCoverage(targetVertex->junctionInRange, meanCoverageList) == false)
		return;

	for (sampleLoop = 0; sampleLoop < SUPPORT_VECTOR_SIZE; ++sampleLoop)
	{
		targetVertex->proportion[sampleLoop] = meanCoverageList[sampleLoop];
	}
	return;
	//filter exons with coverage less than threshold * mean coverage
	curJunc = targetVertex->junctionInRange->list; prevJunc = NULL;
	while (curJunc != NULL)
	{
		if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron || curJunc->junc->type == frag_junction)
		{
			lowCovSampleCnt = 0;
			for (sampleLoop = 0; sampleLoop < SUPPORT_VECTOR_SIZE; ++sampleLoop)
			{
				if (curJunc->junc->support[sampleLoop] < coverageThreshold_GeneExon * meanCoverageList[sampleLoop])
					++lowCovSampleCnt;
			}

			if (lowCovSampleCnt >= thresh_lowCovSampleCnt)
			{
				//filter
				if (prevJunc != NULL)
				{
					prevJunc->next = curJunc->next;
					delete curJunc;
					curJunc = prevJunc->next;
				}
				else
				{
					targetVertex->junctionInRange->list = curJunc->next;
					delete curJunc;
					curJunc = targetVertex->junctionInRange->list;
				}
			}
			else
			{
				prevJunc = curJunc;
				curJunc = curJunc->next;
			}
		}
		else
		{
			prevJunc = curJunc;
			curJunc = curJunc->next;
		}
	}

	//filter junctions connected to nothing
	sortList_Junction_Num = 0;
	curJunc = targetVertex->junctionInRange->list;
	while (curJunc != NULL)
	{
		sortList_Junction[++sortList_Junction_Num] = curJunc;
		curJunc = curJunc->next;
	}
	filterFalseFragments(true);

	//reform the fragment list
	for (iLoop = 1; iLoop <= sortList_Junction_Num; ++iLoop)
	{
		sortKey_Junction[iLoop] = sortList_Junction[iLoop]->junc->end;
	}
	mergeSort_JunctionSort(sortList_Junction_Num);

	for (iLoop = 1; iLoop <= sortList_Junction_Num; ++iLoop)
	{
		sortKey_Junction[iLoop] = sortList_Junction[iLoop]->junc->start;
	}
	mergeSort_JunctionSort(sortList_Junction_Num);

	targetVertex->junctionInRange->list = NULL;
	for (iLoop = sortList_Junction_Num; iLoop >= 1; --iLoop)
	{
		sortList_Junction[iLoop]->next = targetVertex->junctionInRange->list;
		targetVertex->junctionInRange->list = sortList_Junction[iLoop];
	}	

	return;
}

RangeJunctionList* separateDepPath(RangeJunctionList* origList, bool &separable);

void constructGTree(RangeJunctionList* origList)
{
	//construct GTree
	//input: sorted range fragment list
	bool separable_region, separable_path, separable_depPath, separable_gene;
	RangeJunctionList *curList, *curPathList, *gtfPathList, *backupList, *prevPathList, *tempPathList, *curtempList, *curtempList_temp;
	GTvertex *newVertex;
	GTedge *newEdge;

	stack_initial();

	GTvertex *curVertex, *lastVertex;
	curVertex = new GTvertex;
	curVertex->junctionInRange = origList;
	curVertex->level = 0;
	curVertex->rangeLow = CHROMOSOME_START;
	curVertex->rangeHigh = CHROMOSOME_END;
	curVertex->child = NULL;
	curVertex->childType = 1;

	gTree->root = curVertex;

	//separate genes
	separable_region = true;
	curList = NULL;

	curList = separateGene(curVertex->junctionInRange, separable_region);

	if (separable_region == true)
	{
		while (curList != NULL)
		{
			//(1) if two genes are overlapping, e.g., positive strand and negative strand genes on the same locus
			//(2) if there are some standalone stuff, like false region, microRNA, etc.   they are not paths, but may encode some important information
			//therefore, separate them first.
			curPathList = curList;
			curList = curList->nextList;
			curPathList->nextList = NULL;
			separable_path = true;

			///*			
			curPathList = separateIndepPath(curPathList, separable_path);

			if (separable_path == true)
			{
				prevPathList = NULL;
				curtempList_temp = curPathList;
				while (curtempList_temp != NULL)
				{
					tempPathList = SeparateOverlapGene(curtempList_temp, separable_gene);

					curtempList = tempPathList;
					while(curtempList != NULL)
					{
						newVertex = new GTvertex;
						newVertex->level = curVertex->level; //genes, same level as the chromosome, 0
						newVertex->junctionInRange = curtempList;
						newVertex->rangeLow = curtempList->rangeLow;
						newVertex->rangeHigh = curtempList->rangeHigh;
						newVertex->child = NULL;
						newVertex->childType = 1;

						curtempList = curtempList->nextList;
						newVertex->junctionInRange->nextList = NULL;

						newEdge = new GTedge;
						newEdge->linkedVertex = newVertex;
						newEdge->next = curVertex->child;
						curVertex->child = newEdge;

						filterGeneLowCoverageExon(newVertex);
						stack_push(newVertex);						
					}

					curtempList = tempPathList;
					while(curtempList->nextList != NULL)
					{
						curtempList = curtempList->nextList;
					}
					if (prevPathList == NULL)
					{
						curtempList->nextList = curtempList_temp->nextList;
						curPathList = curtempList;
					}
					else
					{
						prevPathList->nextList = tempPathList;
						curtempList->nextList = curtempList_temp->nextList;
					}
					delete curtempList_temp;
					curtempList_temp = curtempList;				
					prevPathList = curtempList_temp;
					curtempList_temp = curtempList_temp->nextList;

				}
			}
			else
				//*/			
			{
				tempPathList = SeparateOverlapGene(curPathList, separable_gene);
				delete curPathList;
				curPathList = tempPathList;
				while(tempPathList != NULL)
				{
					newVertex = new GTvertex;
					newVertex->level = curVertex->level; //genes, same level as the chromosome, 0
					newVertex->junctionInRange = tempPathList;
					newVertex->rangeLow = tempPathList->rangeLow;
					newVertex->rangeHigh = tempPathList->rangeHigh;
					newVertex->child = NULL;
					newVertex->childType = 1;

					tempPathList = tempPathList->nextList;
					newVertex->junctionInRange->nextList = NULL;

					newEdge = new GTedge;
					newEdge->linkedVertex = newVertex;
					newEdge->next = curVertex->child;
					curVertex->child = newEdge;

					filterGeneLowCoverageExon(newVertex);
					stack_push(newVertex);

				}

			}			
		}
	} 
	else
	{
		///*		
		curVertex->junctionInRange = curList;
		curVertex->childType = 1;

		//(1) if two genes are overlapping, e.g., positive strand and negative strand genes on the same locus
		//(2) if there are some standalone stuff, like false region, microRNA, etc.   they are not paths, but may encode some important information
		//therefore, separate them first.
		separable_path = true;
		curList = NULL;

		curList = separateIndepPath(curVertex->junctionInRange, separable_path);

		if (separable_path == true)
		{
			prevPathList = NULL;
			curtempList_temp = curList;
			while (curtempList_temp != NULL)
			{
				tempPathList = SeparateOverlapGene(curtempList_temp, separable_gene);

				curtempList = tempPathList;
				while(curtempList != NULL)
				{
					newVertex = new GTvertex;
					newVertex->level = curVertex->level; //genes, same level as the chromosome, 0
					newVertex->junctionInRange = curtempList;
					newVertex->rangeLow = curtempList->rangeLow;
					newVertex->rangeHigh = curtempList->rangeHigh;
					newVertex->child = NULL;
					newVertex->childType = 1;

					curtempList = curtempList->nextList;
					newVertex->junctionInRange->nextList = NULL;

					newEdge = new GTedge;
					newEdge->linkedVertex = newVertex;
					newEdge->next = curVertex->child;
					curVertex->child = newEdge;

					filterGeneLowCoverageExon(newVertex);
					stack_push(newVertex);					
				}

				curtempList = tempPathList;
				while(curtempList->nextList != NULL)
				{
					curtempList = curtempList->nextList;
				}
				if (prevPathList == NULL)
				{
					curtempList->nextList = curtempList_temp->nextList;
					curList = curtempList;
				}
				else
				{
					prevPathList->nextList = tempPathList;
					curtempList->nextList = curtempList_temp->nextList;
				}
				delete curtempList_temp;
				curtempList_temp = curtempList;				
				prevPathList = curtempList_temp;
				curtempList_temp = curtempList_temp->nextList;

			}

		}
		else
			//*/		
		{
			tempPathList = SeparateOverlapGene(curList, separable_gene);
			delete curList;
			curList = tempPathList;
			if (separable_gene == true)
			{
				while(tempPathList != NULL)
				{
					newVertex = new GTvertex;
					newVertex->level = curVertex->level; //genes, same level as the chromosome, 0
					newVertex->junctionInRange = tempPathList;
					newVertex->rangeLow = tempPathList->rangeLow;
					newVertex->rangeHigh = tempPathList->rangeHigh;
					newVertex->child = NULL;
					newVertex->childType = 1;

					tempPathList = tempPathList->nextList;
					newVertex->junctionInRange->nextList = NULL;

					newEdge = new GTedge;
					newEdge->linkedVertex = newVertex;
					newEdge->next = curVertex->child;
					curVertex->child = newEdge;

					filterGeneLowCoverageExon(newVertex);
					stack_push(newVertex);
				}
			}
			else
			{
				curVertex->junctionInRange = curList;
				filterGeneLowCoverageExon(curVertex);
				stack_push(curVertex);
			}
		}
	}

	//separate ASMs
	curVertex = stack_pop();
	while (curVertex != NULL)
	{
		separable_region = true;
		separable_path = true;
		separable_depPath = true;

		lastVertex = NULL;
		curList = NULL;

		if (curVertex->level < 1)
		{
			curVertex->ID = ++GENEcount;
		}

		if (curVertex->childType == 1)
		{
			//independent regions

			backupList = curVertex->junctionInRange->clone();
			curList = separateIndepRegion(curVertex->junctionInRange, separable_region);

			if (separable_region == true)
			{
				curVertex->junctionInRange = backupList;

				while (curList != NULL)
				{
					newVertex = new GTvertex;
					newVertex->level = curVertex->level + 1;
					newVertex->junctionInRange = curList;
					newVertex->rangeLow = curList->rangeLow;
					newVertex->rangeHigh = curList->rangeHigh;
					newVertex->childType = 2;

					curList = curList->nextList;
					newVertex->junctionInRange->nextList = NULL;

					newEdge = new GTedge;
					newEdge->linkedVertex = newVertex;
					newEdge->next = curVertex->child;
					curVertex->child = newEdge;

					if (lastVertex != NULL)
					{
						lastVertex->prevSibling = newVertex;
						newVertex->nextSibling = lastVertex;
					}
					lastVertex = newVertex;

					stack_push(newVertex);
				}
			} 
			else
			{
				delete backupList;

				curVertex->junctionInRange = curList;
				curVertex->childType = 2;

				stack_push(curVertex);
			}
		} 
		else if (curVertex->childType == 2)
		{
			//independent paths

			/////////////////////////////////////////////////////////////////////////////
			//output gtf tracks for level 1 ASMs 
			//the 2nd condition check whether there are multiple fragments, i.e., an ASM
			if (curVertex->level == 1 && curVertex->junctionInRange->list->next != NULL)
			{		
				curVertex->ID = ++ASMcount;
				gtfPathList = separateDepPath(curVertex->junctionInRange, separable_depPath);
				if (separable_depPath == true)
				{
					//					output_ASMpath_gtf(curVertex, gtfPathList);
					delete gtfPathList;
				}
				else
				{
					//					output_ASMpath_gtf(curVertex, gtfPathList);
					curVertex->junctionInRange = gtfPathList;
				}
			}
			/////////////////////////////////////////////////////////////////////////////


			backupList = curVertex->junctionInRange->clone();
			curList = separateIndepPath(curVertex->junctionInRange, separable_path);

			if (separable_path == true)
			{
				curVertex->junctionInRange = backupList;

				while (curList != NULL)
				{
					newVertex = new GTvertex;
					newVertex->level = curVertex->level + 1;
					newVertex->junctionInRange = curList;
					newVertex->rangeLow = curList->rangeLow;
					newVertex->rangeHigh = curList->rangeHigh;
					newVertex->childType = 1;

					curList = curList->nextList;
					newVertex->junctionInRange->nextList = NULL;

					newEdge = new GTedge;
					newEdge->linkedVertex = newVertex;
					newEdge->next = curVertex->child;
					curVertex->child = newEdge;

					stack_push(newVertex);
				}
			} 
			else
			{
				//reach a leaf
				delete backupList;

				curVertex->junctionInRange = curList;
				curList = separateDepPath(curList, separable_depPath);

				if (separable_depPath == true)
				{
					curVertex->childType = 3;

					while (curList != NULL)
					{
						newVertex = new GTvertex;
						newVertex->level = curVertex->level + 1;
						newVertex->junctionInRange = curList;
						newVertex->rangeLow = curList->rangeLow;
						newVertex->rangeHigh = curList->rangeHigh;
						newVertex->childType = 0;

						curList = curList->nextList;
						newVertex->junctionInRange->nextList = NULL;

						newEdge = new GTedge;
						newEdge->linkedVertex = newVertex;
						newEdge->next = curVertex->child;
						curVertex->child = newEdge;
					}
				} 
				else
				{
					curVertex->junctionInRange = curList;
					curVertex->childType = 0;
				}
			}
		}

		curVertex = stack_pop();
	}

	return;
}


void constructJuncGraph_undirected(RangeJunctionList* origList, bool virtualSE)
{
	//construct fragment graph based on given junctions
	//add virtual start and virtual end if virtualSE is true 

	JuncGraphVertex *curVertex, *tailVertex, *edgeVertex;
	JuncGraphEdge *newEdge;
	rangeJunction *curJunc;

	//build vertices

	curJunc = origList->list;
	while (curJunc != NULL)
	{
		curVertex = new JuncGraphVertex;
		curVertex->corresJunc = curJunc;
		curJunc = curJunc->next;
		curVertex->corresJunc->next = NULL;

		if (junctionGraph->vertices == NULL)
		{
			junctionGraph->vertices = curVertex;
			tailVertex = curVertex;
		} 
		else
		{
			tailVertex->next = curVertex;
			tailVertex = curVertex;
		}
	}

	if (virtualSE == true)
	{
		curVertex = new JuncGraphVertex;
		curVertex->vertexType = virStart;
		curVertex->next = junctionGraph->vertices;
		junctionGraph->vertices = curVertex;

		curVertex = new JuncGraphVertex;
		curVertex->vertexType = virEnd;
		tailVertex->next = curVertex;
		tailVertex = curVertex;
	}

	//build edges

	curVertex = junctionGraph->vertices;
	while (curVertex != NULL)
	{
		if (curVertex->vertexType == normal)
		{
			edgeVertex = curVertex->next;
			while (edgeVertex != NULL)
			{
				if (edgeVertex->vertexType == normal)
				{
					if (compatibleJunctions((curVertex->corresJunc), (edgeVertex->corresJunc)) == true)
					{
						//create an edge for curVertex
						newEdge = new JuncGraphEdge;
						newEdge->linkedVertex = edgeVertex;

						newEdge->next = curVertex->edges;
						curVertex->edges = newEdge;

						//create an edge for edgeVertex
						newEdge = new JuncGraphEdge;
						newEdge->linkedVertex = curVertex;

						newEdge->next = edgeVertex->edges;
						edgeVertex->edges = newEdge;

						edgeVertex->hasInEdge = true;
						curVertex->hasOutEdge = true;
					} 
					else
					{
						//break; 
					}
				}

				edgeVertex = edgeVertex->next;
			}
		}

		curVertex = curVertex->next;
	}

	if (virtualSE == true)
	{
		curVertex = junctionGraph->vertices;
		edgeVertex = curVertex->next;
		while (edgeVertex != NULL)
		{
			if (edgeVertex->vertexType == normal && edgeVertex->hasInEdge == false)
			{
				//create an edge for curVertex
				newEdge = new JuncGraphEdge;
				newEdge->linkedVertex = edgeVertex;

				newEdge->next = curVertex->edges;
				curVertex->edges = newEdge;

				//create an edge for edgeVertex
				newEdge = new JuncGraphEdge;
				newEdge->linkedVertex = curVertex;

				newEdge->next = edgeVertex->edges;
				edgeVertex->edges = newEdge;

			}

			edgeVertex = edgeVertex->next;
		}

		curVertex = tailVertex;
		edgeVertex = junctionGraph->vertices->next;
		while (edgeVertex != NULL)
		{
			if (edgeVertex->vertexType == normal && edgeVertex->hasOutEdge == false)
			{
				//create an edge for curVertex
				newEdge = new JuncGraphEdge;
				newEdge->linkedVertex = edgeVertex;

				newEdge->next = curVertex->edges;
				curVertex->edges = newEdge;

				//create an edge for edgeVertex
				newEdge = new JuncGraphEdge;
				newEdge->linkedVertex = curVertex;

				newEdge->next = edgeVertex->edges;
				edgeVertex->edges = newEdge;

			}

			edgeVertex = edgeVertex->next;
		}
	}


	return;
}


void constructJuncGraph_directed(RangeJunctionList* origList)
{
	//construct fragment graph based on given junctions

	JuncGraphVertex *curVertex, *tailVertex, *edgeVertex;
	JuncGraphEdge *newEdge;
	rangeJunction *curJunc;

	//build vertices

	curJunc = origList->list;
	while (curJunc != NULL)
	{
		curVertex = new JuncGraphVertex;
		curVertex->corresJunc = curJunc;
		curJunc = curJunc->next;
		curVertex->corresJunc->next = NULL;

		if (junctionGraph->vertices == NULL)
		{
			junctionGraph->vertices = curVertex;
			tailVertex = curVertex;
		} 
		else
		{
			tailVertex->next = curVertex;
			tailVertex = curVertex;
		}
	}

	//build edges

	curVertex = junctionGraph->vertices;
	while (curVertex != NULL)
	{
		edgeVertex = curVertex->next;
		while (edgeVertex != NULL)
		{
			if (compatibleJunctions((curVertex->corresJunc), (edgeVertex->corresJunc)) == true)
			{
				//create an edge for curVertex
				newEdge = new JuncGraphEdge;
				newEdge->linkedVertex = edgeVertex;

				newEdge->next = curVertex->edges;
				curVertex->edges = newEdge;

				edgeVertex->hasInEdge = true;
				curVertex->hasOutEdge = true;
			} 
			else
			{
				//break; 
			}
			edgeVertex = edgeVertex->next;
		}

		curVertex = curVertex->next;
	}


	return;
}

//separate independent regions
RangeJunctionList* separateIndepRegion(RangeJunctionList* origList, bool &separable)
{
	//separate independent regions
	//return a set of lists, each of which corresponds to an independent region
	//IMPORTANT: RangeJunctionList is head-inserting, for the convenience of stack

	if (origList == NULL)
	{
		separable = false;
		return NULL;
	}

	junctionGraph = new JuncGraph;
	constructJuncGraph_undirected(origList, true);

	RangeJunctionList *resultList, *curList;
	resultList = NULL;
	curList = NULL;

	rangeJunction *curJunc, *curListTail;
	curJunc = NULL;
	curListTail = NULL;

	JuncGraphVertex *curVertex;
	JuncGraphEdge *curEdge;

	long endBoard = 0, geneEnd = 0;
	int numIndep = 0; //number of independent regions

	//handling alternative start
	curEdge = junctionGraph->vertices->edges;
	while (curEdge != NULL)
	{
		if (curEdge->linkedVertex->corresJunc->junc->start > endBoard)
		{
			endBoard = curEdge->linkedVertex->corresJunc->junc->start;
		}
		curEdge = curEdge->next;
	}

	curVertex = junctionGraph->vertices->next;
	while (curVertex != NULL)
	{
		if (curVertex->vertexType == normal)
		{
			curJunc = curVertex->corresJunc;
			curVertex->corresJunc = NULL;

			if (curJunc->junc->start >= endBoard)
			{
				//omit the first list, it is empty
				if (curList != NULL)
				{
					if (endBoard < MAX_CHR_LENGTH)
						curList->rangeHigh = endBoard;
					else
						curList->rangeHigh = geneEnd;
				}

				curList = NULL;
			}

			if (curList == NULL)
			{
				//get an independent region
				numIndep++;

				curList = new RangeJunctionList;
				curList->nextList = resultList;
				curList->rangeLow = curJunc->junc->start;
				resultList = curList;

				curListTail = NULL;
			}

			//process current fragment
			if (curListTail == NULL)
			{
				curList->list = curJunc;
				curListTail = curJunc;
			}
			else
			{
				curListTail->next = curJunc;
				curListTail = curJunc;
			}
			curListTail->next = NULL;


			if (curVertex->edges->linkedVertex->vertexType == virEnd)
			{
				endBoard = MAX_CHR_LENGTH;
				if (curJunc->junc->end > geneEnd)
				{
					geneEnd = curJunc->junc->end;
				}
			}
			else if (curJunc->junc->end > endBoard)
			{
				endBoard = curJunc->junc->end;
			}
		}

		curVertex = curVertex->next;
	}

	if (curList == NULL)
	{
		cout << "curList == NULL" << endl;
		exit(1);
	} 
	else
	{
		//curList->rangeHigh = origList->rangeHigh;
		if (endBoard < MAX_CHR_LENGTH)
			curList->rangeHigh = endBoard;
		else
			curList->rangeHigh = geneEnd;
	}


	delete junctionGraph;

	if (numIndep == 1)
		separable = false; //cannot be separated into multiple regions
	else if (numIndep > 1)
		separable = true; //can be separated
	else
	{
		cout << "numIndep == 0" << endl;
		exit(1);
	}

	origList->list = NULL;
	delete origList;

	return resultList;
}


RangeJunctionList* separateIndepPath(RangeJunctionList* origList, bool &separable)
{
	//separate independent paths 
	//return a set of lists, each of which corresponds to an independent path
	//return decomposable
	if (origList == NULL)
	{
		separable = false;
		return NULL;
	}

	junctionGraph = new JuncGraph;
	constructJuncGraph_undirected(origList, false);
	origList->list = NULL;

	RangeJunctionList *resultList;
	int numCC = 0;
	resultList = search_conn_comp(numCC, origList->rangeLow, origList->rangeHigh);

	delete junctionGraph;

	if (numCC == 1)
		separable = false; //cannot be separated into multiple paths
	else if (numCC > 1)
		separable = true; //can be separated
	else
	{
		cout << "numCC == 0";
		exit(1);
	}

	origList->list = NULL;
	delete origList;

	return resultList;
}


RangeJunctionList* separateDepPath(RangeJunctionList* origList, bool &separable)
{
	//separate dependent paths 
	//return a set of lists, each of which corresponds to a dependent path
	//return decomposable
	if (origList == NULL)
	{
		separable = false;
		return origList;
	}

	rangeJunction *countJunc;
	long junctionCntinList = 0;
	countJunc = origList->list;
	while (countJunc != NULL)
	{
		if (countJunc->junc->type == frag_junction)
		{
			junctionCntinList++;
		}
		countJunc = countJunc->next;
	}
	if (junctionCntinList > MAXJUNCNUMINDEPPATH)
	{
		unresolvableFile << origList->rangeLow << '\t' << origList->rangeHigh << "\t" << junctionCntinList << endl;
		countJunc = origList->list;
		while (countJunc != NULL)
		{
			if (countJunc->junc->type == frag_exon)
			{
				unresolvableFile << "1\t";
			} 
			else
			{
				unresolvableFile << "0\t";
			}
			unresolvableFile << countJunc->junc->start << "\t" << countJunc->junc->end;
			for (int tmp = 0; tmp < SUPPORT_VECTOR_SIZE; tmp++)
			{
				unresolvableFile << "\t" << countJunc->junc->support[tmp];
			}
			unresolvableFile << endl;

			countJunc = countJunc->next;
		}
		unresolvableFile << endl;

		separable = false;
		return origList;
	}

	RangeJunctionList *backupList;
	backupList = origList->clone();
	junctionGraph = new JuncGraph;
	constructJuncGraph_directed(backupList);
	backupList->list = NULL;

	RangeJunctionList *resultList, *newList;
	rangeJunction *newJunc, *curJunc;
	int numPath = 0;
	long startPosition, endPosition;
	bool pathExtended;

	JuncGraphVertex *curVertex;
	JuncGraphEdge *curEdge;
	JuncGraphPath *curPath, *newPath;
	//find start vertex and end vertex

	startPosition = origList->rangeLow;
	endPosition = origList->rangeHigh;

	resultList = NULL;
	juncPathQueue_initialization();


	//build initial queue
	curVertex = junctionGraph->vertices;
	while (curVertex != NULL)
	{
		//		if (curVertex->corresJunc->junc->start == startPosition)
		if (curVertex->hasInEdge == false)
		{
			newList = new RangeJunctionList;
			newList->rangeLow = curVertex->corresJunc->junc->start;
			newList->rangeHigh = curVertex->corresJunc->junc->end;
			newList->transDirection = curVertex->corresJunc->junc->transDirection;

			newJunc = new rangeJunction;
			newJunc->junc = curVertex->corresJunc->junc;
			newList->list = newJunc;
			newList->listtail = newJunc;

			newPath = new JuncGraphPath;
			newPath->arrivedVertex = curVertex;
			newPath->pathJuncList = newList;

			juncPathQueue_enqueue(newPath);
		}
		curVertex = curVertex->next;
	}

	//main part
	curPath = juncPathQueue_dequeue();
	while (curPath != NULL)
	{
		if (juncPathQueueTail - juncPathQueueHead >= 10 * MAX_NUM_ENUMERATED_PATH || numPath > MAX_NUM_ENUMERATED_PATH)
		{
			separable = false;
			while (curPath = juncPathQueue_dequeue())
			{
				delete curPath;
			}

			break;
		}

		//		if (curPath->arrivedVertex->corresJunc->junc->end == endPosition)
		pathExtended = false;

		if (curPath->arrivedVertex->edges != NULL)
		{
			curEdge = curPath->arrivedVertex->edges;
			while (curEdge != NULL)
			{
				if (curEdge->linkedVertex->corresJunc->junc->end >= curPath->arrivedVertex->corresJunc->junc->end)
				{
					if (curPath->pathJuncList->transDirection == undetermined || curEdge->linkedVertex->corresJunc->junc->transDirection == undetermined
						|| curEdge->linkedVertex->corresJunc->junc->transDirection == curPath->pathJuncList->transDirection)
					{
						pathExtended = true;
						newList = curPath->pathJuncList->clone();

						newJunc = new rangeJunction;
						newJunc->junc = curEdge->linkedVertex->corresJunc->junc;
						//newJunc->next = newList->list;
						//newList->list = newJunc;
						newList->listtail->next = newJunc;
						newList->listtail = newJunc;

						if (newList->transDirection == undetermined && (newJunc->junc->transDirection == sense || newJunc->junc->transDirection == antisense))
						{
							newList->transDirection = newJunc->junc->transDirection;
						}

						if (newJunc->junc->start < newList->rangeLow)
						{
							newList->rangeLow = newJunc->junc->start;
						}
						if (newJunc->junc->end > newList->rangeHigh)
						{
							newList->rangeHigh = newJunc->junc->end;
						}

						newPath = new JuncGraphPath;
						newPath->arrivedVertex = curEdge->linkedVertex;
						newPath->pathJuncList = newList;

						juncPathQueue_enqueue(newPath);
					}
				}

				curEdge = curEdge->next;
			}			
		}

		if (pathExtended == false)
		{
			//current path arrives destination
			numPath++;

			curPath->pathJuncList->nextList = resultList;
			resultList = curPath->pathJuncList;			
		}

		curPath->pathJuncList = NULL;
		delete curPath;

		curPath = juncPathQueue_dequeue();
	}


	delete junctionGraph;

	if (separable == true)
	{
		if (numPath == 1)
			separable = false; //cannot be separated into multiple paths
		else if (numPath > 1)
			separable = true; //can be separated
		else
		{
			cout << origList->rangeLow << " - " << origList->rangeHigh << "  numPath == 0" << endl;
			exit(1);
		}
	}	

	delete backupList;

	if (separable == true)
	{
		return resultList;
	} 
	else
	{
		numPath = 1;
		RangeJunctionList *dellist;
		while(resultList != NULL)
		{
			dellist = resultList;
			resultList = dellist->nextList;
			delete dellist;
		}
		resultList = origList->clone();
		return resultList;
	}
}

bool preCountGTree(GTvertex *rootVertex)
{

	//collect basic counts on GTree for further statistics
	//collect child count and fragment count for every vertex
	//return true for leaves
	GTvertex *curVertex;
	GTedge *curEdge;
	rangeJunction *curJunc;
	double fragSupportSum[SUPPORT_VECTOR_SIZE], exonSupportSum[SUPPORT_VECTOR_SIZE];
	int juncCnt, exonCnt, fragCnt, exonicLength, tmp, iLoop;

	if (rootVertex->child == NULL)
	{
		//leaf

		juncCnt = 0;
		exonCnt = 0;
		fragCnt = 0;
		exonicLength = 0;
		for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
		{
			fragSupportSum[iLoop] = 0;
			exonSupportSum[iLoop] = 0;
		}

		curJunc = rootVertex->junctionInRange->list;
		while (curJunc != NULL)
		{
			fragCnt++;
			if (curJunc->junc->type == frag_junction)
			{
				juncCnt++;
			}
			else if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
			{
				exonCnt++;
				for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
				{
					exonSupportSum[iLoop] += curJunc->junc->support[iLoop];
				}

				exonicLength += curJunc->junc->end - curJunc->junc->start + 1;
			}

			for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
			{
				fragSupportSum[iLoop] += curJunc->junc->support[iLoop];
			}

			curJunc = curJunc->next;
		}
		rootVertex->junctionNum = juncCnt;
		rootVertex->exonNum = exonCnt;
		rootVertex->IndepASMNum = 0;
		rootVertex->exonicLength = exonicLength;

		rootVertex->estimated = false;
		rootVertex->estimate_exonNum = 0; //leaf has nothing to estimate

		for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
		{
			//count the exonic expression
			if (exonCnt > 0)
			{
				rootVertex->support[iLoop] = double(exonSupportSum[iLoop]) / exonCnt;
			} 
			else
			{
				rootVertex->support[iLoop] = 0.0;
			}
		}

		//rootVertex->anovaScore_support = Anova_test(rootVertex->support);

		return true;
	} 
	else
	{
		//extend children
		rootVertex->estimate_exonNum = 0;

		curEdge = rootVertex->child;
		while (curEdge != NULL)
		{
			(rootVertex->childNum)++;
			curVertex = curEdge->linkedVertex;
			if( curVertex->childType == 2 || curVertex->childType == 3)
			{
				// ASM
				rootVertex->IndepASMNum++;
			}
			if (preCountGTree(curVertex) == false)
			{
			}
			if (rootVertex->childType != 3)
			{
				rootVertex->junctionNum += curVertex->junctionNum;
				rootVertex->exonNum += curVertex->exonNum;
				rootVertex->exonicLength += curVertex->exonicLength;

				if (curVertex->estimated == false)
				{
					rootVertex->estimate_exonNum += curVertex->exonNum;
				} 
				else
				{
					rootVertex->estimate_exonNum += 1; //representative exon
				}
			} 

			curEdge = curEdge->next;
		}

		if (rootVertex->childType == 3)
		{
			juncCnt = 0;
			exonCnt = 0;
			fragCnt = 0;
			exonicLength = 0;

			curJunc = rootVertex->junctionInRange->list;
			while (curJunc != NULL)
			{
				fragCnt++;
				if (curJunc->junc->type == frag_junction)
				{
					juncCnt++;
				}
				else if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
				{
					exonCnt++;
					exonicLength += curJunc->junc->end - curJunc->junc->start + 1;
				}
				curJunc = curJunc->next;
			}
			rootVertex->junctionNum = juncCnt;
			rootVertex->exonNum = exonCnt;
			rootVertex->exonicLength = exonicLength;

			rootVertex->estimate_exonNum = exonCnt;

		}

		return false;
	}
}


/************************************************************************/
/* For Assembly */
/************************************************************************/
double WeibullDist(double gaplength, double scale, double genelength)
{
	double prob, shape;
	shape = log(genelength);

	prob = log(shape) - log(scale) + (shape - 1) * log(gaplength/scale) - pow((gaplength/scale), shape);
	if (prob < -10e5)
	{
		prob = -10e5;
	}
	return prob;
}

double prop1, prop2, mean1, mean2, sd1, sd2; // parameters for 2 mixture Gaussian distribution

// likelihood of a 2 mixture normal distribution
double MixNormal(double gaplength)
{
	double prob;
	if (sd1 == 0 && sd2 > 0)
	{
		if (gaplength == mean1)
		{
			prob = prop1 + prop2 * exp((-1) * pow(gaplength - mean2, 2) / (2 * pow(sd2, 2))) / (sd2 * sqrt(2 * pi));
			if (prob > 1)
			{
				prob = 1;
			}
			prob = log(prob);
		}
		else
			prob = log(prop2) + (-1) * pow(gaplength - mean2, 2) / (2 * pow(sd2, 2)) - log(sd2 * sqrt(2 * pi));
	}
	else if (sd2 == 0 && sd1 > 0)
	{
		if (gaplength == mean2)
		{
			prob = prop1 * exp((-1) * pow(gaplength - mean1, 2) / (2 * pow(sd1, 2))) / (sd1 * sqrt(2 * pi)) + prop2;
			if (prob > 1)
			{
				prob = 1;
			}
			prob = log(prob);
		}
		else
			prob = log(prop1) + (-1) * pow(gaplength - mean1, 2) / (2 * pow(sd1, 2)) - log(sd1 * sqrt(2 * pi));
	}
	else if (sd1 == 0 && sd2 == 0)
	{
		if (gaplength == mean1)
		{
			prob = prop1;
		}
		else if (gaplength == mean2)
		{
			prob = prop2;
		}
		else
			prob = epsilon;
		prob = log(prob);
	}
	else
	{
		prob = prop1 * exp((-1) * pow(gaplength - mean1, 2) / (2 * pow(sd1, 2))) / (sd1 * sqrt(2 * pi)) + prop2 * exp((-1) * pow(gaplength - mean2, 2) / (2 * pow(sd2, 2))) / (sd2 * sqrt(2 * pi));
		if (prob > 1)
		{
			prob = 1;
		}
		prob = log(prob);
		if (prob < -MAX/4)
		{
			prob = -MAX/4;
		}
	}

	return prob;
}

// double MixNormal(double gaplength)
// {
// 	double prob, temp1, temp2;
// 	temp1 = (-1) * pow(gaplength - mean1, 2) / (2 * pow(sd1, 2)) - log(sd1 * sqrt(2 * pi));
// 	temp2 = (-1) * pow(gaplength - mean2, 2) / (2 * pow(sd2, 2)) - log(sd2 * sqrt(2 * pi));
// 
// 	prob = temp1 > temp2 ? temp1 : temp2;
// 	return prob;
// }

// likelihood of a normal distribution
double Normal(double gaplength, double mean, double sd)
{
	double prob;
	prob = (-1) * pow(gaplength - mean, 2) / (2 * pow(sd, 2)) - log(sd * sqrt(2 * pi));

	if (prob < -MAX/2)
	{
		prob = -MAX/2;
	}
	return prob;
}

// function used for shrinkage the probability one transcript is true (exp{-cx})
double PenaltyFunc(double prob, long readNum, long transNum)
{
	double value;
	value = prob/exp((double) readNum/transNum);
	value = log(value);
	return value;
}

trans_direction Getdirection(GTvertex *vertex)
{
	trans_direction direction;
	direction = undetermined;
	rangeJunction* frags = vertex->junctionInRange->list;
	while(frags != NULL)
	{
		if (frags->junc->type == frag_junction)
		{
			direction = frags->junc->transDirection;
			break;
		}
		frags = frags->next;
	}
	return direction;
}

// pre-process the read file
void preprocessReadFile(GTvertex *rootVertex, string outputpath, string SAMPath, string SAMFile_prefix)
{
	string filename;
	ifstream inputfile;
	vector<ofstream*> outputfile_SAM;

	GTedge *curEdge, *nextEdge;
	GTvertex *curVertex, *nextVertex;
	long index, OverlapNum;
	trans_direction direction;

	// information in read file
	char cigar[500], field[500], strand[10];
	int flag_read;
	long i, tmp, readlen, startRead, endRead, startPoint, endPoint, readLength;
	string s;

	filename = SAMPath + SAMFile_prefix + ".txt";
	inputfile.open(filename.c_str());
	input_read = true;
	field[0] = '\0';
	inputfile >> field;	

	if ( (rootVertex->child == NULL) || (rootVertex->level == 0 && rootVertex->child->linkedVertex->level > 0) )
	{
		// leaf or root is a gene

		curVertex = rootVertex;
		ofstream *temp = new ofstream;
		if (outputfile_SAM.size() >= outputfile_SAM.capacity())
			outputfile_SAM.reserve(default_dataset_num + outputfile_SAM.capacity());
		outputfile_SAM.push_back(temp);

		direction = Getdirection(curVertex);
		if (direction == antisense)
		{
			filename = outputpath + itostr(curVertex->rangeLow) + "-" + itostr(curVertex->rangeHigh) + "C.txt";
		}
		else
			filename = outputpath + itostr(curVertex->rangeLow) + "-" + itostr(curVertex->rangeHigh) + "W.txt";
		(*outputfile_SAM[outputfile_SAM.size()-1]).open(filename.c_str(), fstream::in | fstream::out | fstream::app);

		while(field[0] != '\0')
		{
			startPoint = atol(field);
			if (input_read == true)
			{
				inputfile >> cigar;
				inputfile >> strand;

				// Get the readlength and the endPoint of this read
				tmp = 0;
				readlen = 0;

				startRead = startPoint;
				endRead = startRead - 1;

				for (i = 0; cigar[i] != '\0'; i++)
				{
					if (cigar[i] == 'M')
					{
						endRead = endRead + tmp;
						readlen += tmp;
						tmp = 0;
					} 
					else if (cigar[i] == 'N')
					{
						endRead = endRead + tmp;
						tmp = 0;
					}
					else if (cigar[i] >= '0' && cigar[i] <= '9')
					{
						tmp = tmp * 10 + cigar[i] - 48;
					}
					else
					{
						tmp = 0;
					}
				}
				endPoint = endRead;
				readLength = readlen;

				getline(inputfile, s);
			}

			if (endPoint == startPoint - 1)
			{
				input_read = true;
			}

			if (endPoint <= curVertex->rangeLow)
			{
				input_read = true;
			}
			else if ((startPoint <= curVertex->rangeLow && endPoint >= curVertex->rangeLow && endPoint <= curVertex->rangeHigh) || (startPoint <= curVertex->rangeLow && endPoint >= curVertex->rangeHigh) || (startPoint >= curVertex->rangeLow && endPoint <= curVertex->rangeHigh) || (startPoint > curVertex->rangeLow && startPoint < curVertex->rangeHigh && endPoint > curVertex->rangeHigh))
			{
				input_read = true;
				index = 0;
				(*outputfile_SAM[index]) << startPoint << "\t" << cigar << "\t" << strand << endl;
			}
			else if (startPoint >= curVertex->rangeHigh)
			{
				input_read = false;
				break;
			}

			if (input_read == true)
			{
				field[0] = '\0';
				inputfile >> field;
			}
		}

		for (index = 0; index < outputfile_SAM.size(); index++)
		{
			(*outputfile_SAM[index]).close();
		}
		outputfile_SAM.clear();
	}
	else
	{
		// extend children
		curEdge = rootVertex->child;
		while(curEdge != NULL)
		{
			curVertex = curEdge->linkedVertex;
			if (curVertex != NULL && curVertex->level == 0)
			{
				// reach a gene
				//				cout << curVertex->rangeLow << "-" << curVertex->rangeHigh << endl;
				ofstream *temp = new ofstream;
				if (outputfile_SAM.size() >= outputfile_SAM.capacity())
					outputfile_SAM.reserve(default_dataset_num + outputfile_SAM.capacity());
				outputfile_SAM.push_back(temp);

				direction = Getdirection(curVertex);
				if (direction == antisense)
				{
					filename = outputpath + itostr(curVertex->rangeLow) + "-" + itostr(curVertex->rangeHigh) + "C.txt";
				}
				else
					filename = outputpath + itostr(curVertex->rangeLow) + "-" + itostr(curVertex->rangeHigh) + "W.txt";
				(*outputfile_SAM[outputfile_SAM.size()-1]).open(filename.c_str(), fstream::in | fstream::out | fstream::app);

				nextEdge = curEdge->next;
				while(nextEdge != NULL)
				{
					nextVertex = nextEdge->linkedVertex;
					if (nextVertex != NULL && nextVertex->level == 0)
					{
						if (nextVertex->rangeLow > curVertex->rangeHigh)
						{
							break;
						}
						else if (nextVertex->rangeLow >= curVertex->rangeLow && nextVertex->rangeLow <= curVertex->rangeHigh)
						{
							ofstream *temp = new ofstream;
							if (outputfile_SAM.size() >= outputfile_SAM.capacity())
								outputfile_SAM.reserve(default_dataset_num + outputfile_SAM.capacity());
							outputfile_SAM.push_back(temp);

							direction = Getdirection(nextVertex);
							if (direction == antisense)
							{
								filename = outputpath + itostr(nextVertex->rangeLow) + "-" + itostr(nextVertex->rangeHigh) + "C.txt";
							}
							else
								filename = outputpath + itostr(nextVertex->rangeLow) + "-" + itostr(nextVertex->rangeHigh) + "W.txt";
							(*outputfile_SAM[outputfile_SAM.size()-1]).open(filename.c_str(), fstream::in | fstream::out | fstream::app);

						}
					}
					nextEdge = nextEdge->next;
				}

				while(field[0] != '\0')
				{
					startPoint = atol(field);
					if (input_read == true)
					{
						inputfile >> cigar;
						inputfile >> strand;

						// Get the readlength and the endPoint of this read
						tmp = 0;
						readlen = 0;

						startRead = startPoint;
						endRead = startRead - 1;

						for (i = 0; cigar[i] != '\0'; i++)
						{
							if (cigar[i] == 'M')
							{
								endRead = endRead + tmp;
								readlen += tmp;
								tmp = 0;
							} 
							else if (cigar[i] == 'N')
							{
								endRead = endRead + tmp;
								tmp = 0;
							}
							else if (cigar[i] >= '0' && cigar[i] <= '9')
							{
								tmp = tmp * 10 + cigar[i] - 48;
							}
							else
							{
								tmp = 0;
							}
						}
						endPoint = endRead;
						readLength = readlen;

						getline(inputfile, s);
					}

					if (endPoint == startPoint - 1)
					{
						input_read = true;
					}

					if (endPoint < curVertex->rangeLow)
					{
						input_read = true;
					}
					else if ((startPoint <= curVertex->rangeLow && endPoint >= curVertex->rangeLow && endPoint <= curVertex->rangeHigh) || (startPoint <= curVertex->rangeLow && endPoint >= curVertex->rangeHigh) || (startPoint >= curVertex->rangeLow && endPoint <= curVertex->rangeHigh) || (startPoint > curVertex->rangeLow && startPoint < curVertex->rangeHigh && endPoint > curVertex->rangeHigh))
					{
						input_read = true;
						index = 0;
						(*outputfile_SAM[index]) << startPoint << "\t" << cigar << "\t" << strand << endl;
						if (outputfile_SAM.size() > 0)
						{
							nextEdge = curEdge->next;
							while(nextEdge != NULL)
							{
								nextVertex = nextEdge->linkedVertex;
								if (nextVertex != NULL && nextVertex->level == 0)
								{
									if (nextVertex->rangeLow > curVertex->rangeHigh)
									{
										break;
									}
									else if (nextVertex->rangeLow >= curVertex->rangeLow && nextVertex->rangeLow <= curVertex->rangeHigh)
									{
										index++;
										if ((startPoint <= nextVertex->rangeLow && endPoint >= nextVertex->rangeLow && endPoint <= nextVertex->rangeHigh) || (startPoint <= nextVertex->rangeLow && endPoint >= nextVertex->rangeHigh) || (startPoint >= nextVertex->rangeLow && endPoint <= nextVertex->rangeHigh) || (startPoint > nextVertex->rangeLow && startPoint < nextVertex->rangeHigh && endPoint > nextVertex->rangeHigh))
										{
											(*outputfile_SAM[index]) << startPoint << "\t" << cigar << "\t" << strand << endl;
										}
									}
								}
								nextEdge = nextEdge->next;
							}
						}

					}
					else if (startPoint >= curVertex->rangeHigh)
					{
						input_read = false;
						break;
					}

					if (input_read == true)
					{
						field[0] = '\0';
						inputfile >> field;
					}
				}

				for (index = 0; index < outputfile_SAM.size(); index++)
				{
					(*outputfile_SAM[index]).close();
				}
				while(outputfile_SAM.empty() != true)
				{
					delete outputfile_SAM[outputfile_SAM.size()-1];
					outputfile_SAM.pop_back();
				}
			}
			curEdge = curEdge->next;
		}
	}

	inputfile.close();

	return;
}


long sourceDegree; // document the maximum of edges should flow out of source (flow in sink)
long MaxReadLength = -1, MeanReadLength;
bool processReadFile(GTvertex *vertex, string filename)
{
	bool exceedMax = false;

	ifstream inputfile;
	char cigar[500], field[500];
	int flag_read, strand;
	long i, tmp, readlen, startRead, endRead, startPoint, endPoint, readLength, junctionNum, junctionIndex;
	bool spliced;
	string s, tmpstr;
	ReadVertex *newVertex;
	rangeJunction *fragList, *curfrag, *newfrag;
	long *readCnt;
	readCnt = new long[GeneEnd - GeneStart + 1];
	for (i = 0; i < GeneEnd - GeneStart + 1; i++)
	{
		readCnt[i] = 0;
	}
	sourceDegree = - MAX;
	MeanReadLength = 0;

	fragList = vertex->junctionInRange->list;

	inputfile.open(filename.c_str());
	field[0] = '\0';
	inputfile >> field;

	while(field[0] != '\0')
	{
		startPoint = atol(field);
		inputfile >> cigar;
		inputfile >> field;
		getline(inputfile, s);

		if (strcmp(field, "+") == 0)
		{
			strand = 1;
		}
		else
		{
			strand = 0;
		}
		spliced = false;

		// Get the readlength and the endPoint of this read
		tmp = 0;
		readlen = 0;

		startRead = startPoint;
		endRead = startRead - 1;

		for (i = 0; cigar[i] != '\0'; i++)
		{
			if (cigar[i] == 'M')
			{
				endRead = endRead + tmp;
				readlen += tmp;
				tmp = 0;
			} 
			else if (cigar[i] == 'N')
			{
				spliced = true;
				endRead = endRead + tmp;
				tmp = 0;
			}
			else if (cigar[i] >= '0' && cigar[i] <= '9')
			{
				tmp = tmp * 10 + cigar[i] - 48;
			}
			else
			{
				tmp = 0;
			}
		}

		endPoint = endRead;
		readLength = readlen;		
		if ((startPoint <= GeneStart && endPoint >= GeneStart && endPoint <= GeneEnd) || (startPoint <= GeneStart && endPoint >= GeneEnd) || (startPoint >= GeneStart && endPoint <= GeneEnd) || (startPoint > GeneStart && startPoint < GeneEnd && endPoint > GeneEnd))
		{
			newVertex = new ReadVertex;
			tmpstr = cigar;
			newVertex->readname = itostr(startPoint) + "_" + tmpstr;
			newVertex->start = startPoint;
			newVertex->end = endPoint;
			newVertex->readLength = readlen;
			newVertex->spliced = spliced;
			//			strcpy(newVertex->cigar, cigar);

			if (strand == 0)
			{
				newVertex->transDirection = antisense;
			}
			else if (strand == 1)
			{
				newVertex->transDirection = sense;
			}

			// find the genomic region where the read aligns		

			bool compatible = false;

			curfrag = fragList;
			junctionIndex = 0;
			if (curfrag != NULL)
			{
				tmp = 0;
				startRead = startPoint;
				endRead = startRead;
				for (i = 0; cigar[i] != '\0'; i++)
				{
					if (cigar[i] == 'M')
					{
						compatible = false;

						endRead = endRead + tmp;
						while(curfrag != NULL)
						{
							if ((startRead >= curfrag->junc->start && (endRead-1) <= curfrag->junc->end) || ((endRead-1) >= curfrag->junc->start && (endRead-1) <= curfrag->junc->end && startRead < curfrag->junc->start))
							{
								if (curfrag->junc->type == frag_exon || curfrag->junc->type == frag_retained_intron)
								{
									newfrag = new rangeJunction;
									newfrag->junc = curfrag->junc->clone();
									newfrag->next = NULL;

									if (newVertex->list == NULL)
									{
										newVertex->list = newfrag;
										newVertex->listTail = newfrag;
									}
									else
									{
										newVertex->listTail->next = newfrag;
										newVertex->listTail = newfrag;
									}
									compatible = true;
									curfrag = curfrag->next;
									if (curfrag != NULL)
									{
										if (curfrag->junc->type == frag_junction)
										{
											junctionIndex++;
										}
									}
									break;
								}
								else
								{
									curfrag = curfrag->next;
									if (curfrag != NULL)
									{
										if (curfrag->junc->type == frag_junction)
										{
											junctionIndex++;
										}
									}
								}
							}
							else if ((startRead >= curfrag->junc->start && startRead <= curfrag->junc->end && (endRead-1) > curfrag->junc->end) || (startRead < curfrag->junc->start && (endRead-1) > curfrag->junc->end))
							{
								if (curfrag->junc->type == frag_exon || curfrag->junc->type == frag_retained_intron)
								{
									newfrag = new rangeJunction;
									newfrag->junc = curfrag->junc->clone();
									newfrag->next = NULL;

									if (newVertex->list == NULL)
									{
										newVertex->list = newfrag;
										newVertex->listTail = newfrag;
									}
									else
									{
										newVertex->listTail->next = newfrag;
										newVertex->listTail = newfrag;
									}
									compatible = false;
									curfrag = curfrag->next;
									if (curfrag != NULL)
									{
										if (curfrag->junc->type == frag_junction)
										{
											junctionIndex++;
										}
									}
								}
								else
								{
									curfrag = curfrag->next;
									if (curfrag != NULL)
									{
										if (curfrag->junc->type == frag_junction)
										{
											junctionIndex++;
										}
									}
								}

							}
							else if (curfrag->junc->end <= startRead)
							{
								curfrag = curfrag->next;
								if (curfrag != NULL)
								{
									if (curfrag->junc->type == frag_junction)
									{
										junctionIndex++;
									}
								}
							}
							else if (curfrag->junc->start >= (endRead-1))
							{
								break;
							}

						}
						if (compatible == false)
						{
							break;
						}
						startRead = endRead;
						tmp = 0;
					} 
					else if (cigar[i] == 'N')
					{
						compatible = false;
						endRead = endRead + tmp;
						while(curfrag != NULL)
						{
							if ((startRead-1) == curfrag->junc->start && endRead == curfrag->junc->end && curfrag->junc->type == frag_junction)
							{
								newfrag = new rangeJunction;
								newfrag->junc = curfrag->junc->clone();
								newfrag->next = NULL;

								newVertex->listTail->next = newfrag;
								newVertex->listTail = newfrag;
								newVertex->junctionsIndex.push_back(junctionIndex - 1);

								curfrag = curfrag->next;
								if (curfrag != NULL)
								{
									if (curfrag->junc->type == frag_junction)
									{
										junctionIndex++;
									}
								}
								compatible = true;
								break;
							}
							else
							{
								curfrag = curfrag->next;
								if (curfrag != NULL)
								{
									if (curfrag->junc->type == frag_junction)
									{
										junctionIndex++;
									}
								}
							}
						}
						if (compatible == false)
						{
							break;
						}
						startRead = endRead;
						tmp = 0;
					}
					else if (cigar[i] >= '0' && cigar[i] <= '9')
					{
						tmp = tmp * 10 + cigar[i] - 48;
					}
					else
					{
						tmp = 0;
					}
				}
			}

			// check one more time for the compatibility
			if (compatible == true)
			{
				if (newVertex->start < newVertex->list->junc->start || newVertex->end > newVertex->listTail->junc->end)
				{
					compatible = false;
				}
				curfrag = newVertex->list;
				while(curfrag->next != NULL)
				{
					if ((curfrag->junc->type == frag_junction && (curfrag->next->junc->type == frag_exon || curfrag->next->junc->type == frag_retained_intron)) || ((curfrag->junc->type == frag_exon || curfrag->junc->type == frag_retained_intron) && curfrag->next->junc->type == frag_junction))
					{
						if (curfrag->junc->end != curfrag->next->junc->start)
						{
							compatible = false;
							break;
						}
					}
					else
					{
						if (abs(curfrag->junc->end - curfrag->next->junc->start) > 1)
						{
							compatible = false;
							break;
						}
					}
					curfrag = curfrag->next;
				}
				if (!( ( abs(newVertex->list->junc->start - newVertex->start) <= 1 || (newVertex->list->junc->start <= newVertex->start) ) && ( abs(newVertex->listTail->junc->end - newVertex->end) <= 1 ) || (newVertex->listTail->junc->end >= newVertex->end) ))
				{
					compatible = false;
				}
			}

			if (compatible == false)
			{
				//				cout << "alert!!!" << endl;

				delete newVertex;
			}
			else
			{
				// add new read to the list
				if (ReadList.size() >= ReadList.capacity())
					ReadList.reserve(default_dataset_num + ReadList.capacity());
				ReadList.push_back(newVertex);
				ReadList[ReadList.size() - 1]->Id = ReadList.size();
				newVertex->start_min = newVertex->start;
				newVertex->end_max = newVertex->end;
				// get the read length distribution on this genomic region
				gapConstraint += readLength;
				MeanReadLength += readLength;
				if (readLength > MaxReadLength)
				{
					MaxReadLength = readLength;
				}

// 				if (ReadList.size() > MAX_DOABLE_READNUM)
// 				{
// 					exceedMax = true;
// 					break;
// 				}
				long index_start = newVertex->start - GeneStart, index_end = newVertex->end - GeneStart;
				if (index_start < 0)
				{
					index_start = 0;
				}
				if (index_end > (GeneEnd - GeneStart))
				{
					index_end = GeneEnd - GeneStart;
				}
				for (i = index_start; i <= index_end; i++)
				{
					readCnt[i]++;
				}	
			}
		}
		else if (startPoint > GeneEnd)
		{
			break;
		}

		field[0] = '\0';
		inputfile >> field;

	}

	inputfile.close();

	sourceDegree = readCnt[0];
	for (i = 1; i < GeneEnd - GeneStart + 1; i++)
	{
		if (readCnt[i] > sourceDegree)
		{
			sourceDegree = readCnt[i];
		}
	}
	delete [] readCnt;

	if (ReadList.size() > 0)
	{
		MeanReadLength /= ReadList.size();
	}
	
	return exceedMax;
}

graphpath* generate_graphpath(rangeJunction *fragstart, rangeJunction *fragend, bool transcriptionEnds);

bool IdenticalJunctions(GTvertex *vertex, ReadVertex *vertex1, ReadVertex *vertex2)
{
	rangeJunction *curfrag = vertex->junctionInRange->list, *startfrag, *endfrag;
	graphpath *paths, *del_path;
	JuncGraphEdge * curedge;
	bool Identical = true, found;
	long start, end;
	if (vertex1->junctionsIndex.size() != vertex2->junctionsIndex.size())
	{
		Identical = false;
	}
	else
	{
		for (int index = 0; index < vertex1->junctionsIndex.size(); index++)
		{
			if (vertex1->junctionsIndex[index] != vertex2->junctionsIndex[index])
			{
				Identical = false;
				break;
			}
		}
	}
	// two reads has to be overlapped
	if (Identical == true)
	{
// 		if (vertex2->end_max < vertex1->start_min || vertex2->start_min > vertex1->end_max)
// 		{
// 			// check whether there is multiple paths between these two reads
// 			Identical = false;
// 			if (vertex2->end_max < vertex1->start_min)
// 			{
// 				start = vertex2->end_max;
// 				end = vertex1->start_min;
// 			}
// 			else if (vertex2->start_min > vertex1->end_max)
// 			{
// 				start = vertex1->end_max;
// 				end = vertex2->start_min;
// 			}
// 			found = false;
// 			while(curfrag != NULL)
// 			{
// 				if (found == true)
// 				{
// 					if ((curfrag->junc->type == frag_exon || curfrag->junc->type == frag_retained_intron) && end >= curfrag->junc->start && end <= curfrag->junc->end)
// 					{
// 						endfrag = curfrag;
// 						break;
// 					}
// 				}
// 				if ((curfrag->junc->type == frag_exon || curfrag->junc->type == frag_retained_intron) && start >= curfrag->junc->start && start <= curfrag->junc->end)
// 				{
// 					startfrag = curfrag;
// 					found = true;
// 				}
// 				curfrag = curfrag->next;
// 			}
// 			if (startfrag != NULL && endfrag != NULL)
// 			{
// 				paths = generate_graphpath(startfrag, endfrag, false);
// 				if (paths != NULL)
// 				{
// 					if (paths->next == NULL)
// 					{
// 						Identical = true;
// 						curedge = paths->edgelist;
// 						while(curedge->next != NULL)
// 						{
// 							if (curedge->linkedVertex->corresJunc->junc->type == frag_junction)
// 							{
// 								Identical = false;
// 								break;
// 							}
// 							curedge = curedge->next;
// 						}
// 					}
// 				}
// 				while(paths != NULL)
// 				{
// 					del_path = paths;
// 					paths = del_path->next;
// 					delete del_path;
// 				}
// 			}
// 
// 		}
		if (vertex2->end_max < vertex1->start_min || vertex2->start_min > vertex1->end_max)
		{
			Identical = false;
		}

	}


	return Identical;
}


// cluster the compatible read fragments within one subgroup
void ClusterReadFragments_subgroup(GTvertex *vertex, long pointer, long pointer_end, int tolerance)
{
	long index, fragmentIndex, iLoop, shift;
	bool found;
	rangeJunction *curfrag, *curfrag_temp, *newfrag, *list, *tail;
	ReadVertex *newVertex;
	shift = tolerance;

// 	ofstream outputfile;
// 	outputfile.open("yan.txt", std::fstream::in | std::fstream::out | std::fstream::app);

	// build the first cluster within this subgroup
	newVertex = ReadList[pointer]->clone();
	newVertex->start_min = newVertex->start;
	newVertex->end_max = newVertex->end;
	if (ReadList_cluster.size() >= ReadList_cluster.capacity())
		ReadList_cluster.reserve(default_dataset_num + ReadList_cluster.capacity());
	ReadList_cluster.push_back(newVertex);
	ReadList_cluster[ReadList_cluster.size() - 1]->includedReadNum = 1;
	index = ReadList_cluster.size() - 1;
	index = (index - shift) > 0? (index - shift) : 0;

//	outputfile << newVertex->readname << "\t" << ReadList_cluster.size() - 1 << "\t" << newVertex->start << "-" << newVertex->end << endl;

	// parse through the fragments within this subgroup
	for (fragmentIndex = pointer + 1; fragmentIndex <= pointer_end; fragmentIndex++)
	{
		found = false;
		for (iLoop = index; iLoop < ReadList_cluster.size(); iLoop++)
		{
			if (IdenticalJunctions(vertex, ReadList[fragmentIndex], ReadList_cluster[iLoop]) == true)
			{
				if (abs(ReadList_cluster[iLoop]->start - ReadList[fragmentIndex]->start) < tolerance && abs(ReadList_cluster[iLoop]->end - ReadList[fragmentIndex]->end) < tolerance)
				{
					found = true;

//					outputfile << ReadList[fragmentIndex]->readname << "\t" << iLoop << "\t" << ReadList[fragmentIndex]->start << "-" << ReadList[fragmentIndex]->end << endl;

					break;
				}
			}
		}
		if(found == true)
		{
			ReadList_cluster[iLoop]->includedReadNum++;
			if (ReadList[fragmentIndex]->end > ReadList_cluster[iLoop]->end_max)
			{
				// check whether need to extend the last rangeJunction
				if (ReadList[fragmentIndex]->listTail->junc->start >= ReadList_cluster[iLoop]->listTail->junc->end)
				{
					curfrag = ReadList_cluster[iLoop]->list;
					while(curfrag->next != NULL)
					{
						curfrag = curfrag->next;
					}

					if (ReadList[fragmentIndex]->start_min > ReadList_cluster[iLoop]->end_max)
					{
						rangeJunction *genefrag = vertex->junctionInRange->list;
						while(genefrag != NULL)
						{
							if (genefrag->junc->start == curfrag->junc->start && genefrag->junc->end == curfrag->junc->end)
							{
								break;
							}
							genefrag = genefrag->next;
						}
						genefrag = genefrag->next;
						curfrag_temp = ReadList[fragmentIndex]->list;
						while(genefrag != NULL)
						{
							if (genefrag->junc->start == curfrag_temp->junc->start && genefrag->junc->end == curfrag_temp->junc->end)
							{
								break;
							}
							newfrag = new rangeJunction;
							newfrag->junc = genefrag->junc->clone();
							curfrag->next = newfrag;
							curfrag = newfrag;
							ReadList_cluster[iLoop]->listTail = newfrag;
							genefrag = genefrag->next;
						}
					}
					else if (abs(ReadList[fragmentIndex]->start_min - ReadList_cluster[iLoop]->end_max) <= 1)
					{
						curfrag_temp = ReadList[fragmentIndex]->list;
					}
					else
					{
						curfrag_temp = ReadList[fragmentIndex]->list;
						while(curfrag_temp != NULL)
						{
							if (curfrag_temp->junc->start == curfrag->junc->start && curfrag_temp->junc->end == curfrag->junc->end)
							{
								break;
							}
							curfrag_temp = curfrag_temp->next;
						}
						curfrag_temp = curfrag_temp->next;
					}
					while(curfrag_temp != NULL)
					{
						newfrag = new rangeJunction;
						newfrag->junc = curfrag_temp->junc->clone();
						curfrag->next = newfrag;
						curfrag = newfrag;
						ReadList_cluster[iLoop]->listTail = newfrag;
						curfrag_temp = curfrag_temp->next;
					}
				}
				ReadList_cluster[iLoop]->end_max = ReadList[fragmentIndex]->end;
			}
			if (ReadList[fragmentIndex]->start < ReadList_cluster[iLoop]->start_min)
			{
				// check whether need to extend the first rangeJunction
				if (ReadList[fragmentIndex]->list->junc->end <= ReadList_cluster[iLoop]->list->junc->start)
				{
					list = tail = NULL;
					curfrag_temp = ReadList[fragmentIndex]->list;
					while(curfrag_temp->next != NULL)
					{
						newfrag = new rangeJunction;
						newfrag->junc = curfrag_temp->junc->clone();
						if (list == NULL)
						{
							list = newfrag;
							tail = newfrag;
						}
						else
						{
							tail->next = newfrag;
							tail = newfrag;
						}
						curfrag_temp = curfrag_temp->next;
					}
					curfrag = ReadList_cluster[iLoop]->list;
					if (ReadList[fragmentIndex]->end_max < ReadList_cluster[iLoop]->start_min)
					{
						rangeJunction *genefrag = vertex->junctionInRange->list;
						bool found = false;
						while(genefrag != NULL)
						{
							if (genefrag->junc->start == curfrag_temp->junc->start && genefrag->junc->end == curfrag_temp->junc->end)
							{
								found = true;
							}
							if (found == true)
							{
								newfrag = new rangeJunction;
								newfrag->junc = genefrag->junc->clone();
								if (list == NULL)
								{
									list = newfrag;
									tail = newfrag;
								}
								else
								{
									tail->next = newfrag;
									tail = newfrag;
								}
							}
							if (genefrag->junc->start == curfrag->junc->start && genefrag->junc->end == curfrag->junc->end)
							{
								break;
							}
							genefrag = genefrag->next;
						}
						curfrag = curfrag->next;
					}
					else if (abs(ReadList[fragmentIndex]->end_max - ReadList_cluster[iLoop]->start_min) <= 1)
					{
						newfrag = new rangeJunction;
						newfrag->junc = curfrag_temp->junc->clone();
						if (list == NULL)
						{
							list = newfrag;
							tail = newfrag;
						}
						else
						{
							tail->next = newfrag;
							tail = newfrag;
						}
					}
					else
					{
						while(curfrag != NULL)
						{
							if (curfrag->junc->start == curfrag_temp->junc->start && curfrag->junc->end == curfrag_temp->junc->end)
							{
								break;
							}
							curfrag = curfrag->next;
						}
					}
					while(curfrag != NULL)
					{
						newfrag = new rangeJunction;
						newfrag->junc = curfrag->junc->clone();
						if (list == NULL)
						{
							list = newfrag;
							tail = newfrag;
						}
						else
						{
							tail->next = newfrag;
							tail = newfrag;
						}
						curfrag = curfrag->next;
					}
					while(ReadList_cluster[iLoop]->list != NULL)
					{
						curfrag = ReadList_cluster[iLoop]->list;
						ReadList_cluster[iLoop]->list = curfrag->next;
						delete curfrag;
					}
					ReadList_cluster[iLoop]->list = list;
					ReadList_cluster[iLoop]->listTail = tail;
				}
				ReadList_cluster[iLoop]->start_min = ReadList[fragmentIndex]->start;
			}
			if (ReadList[fragmentIndex]->end < ReadList_cluster[iLoop]->end)
			{
				ReadList_cluster[iLoop]->end = ReadList[fragmentIndex]->end;
			}
			if (ReadList[fragmentIndex]->start > ReadList_cluster[iLoop]->start)
			{
				ReadList_cluster[iLoop]->start = ReadList[fragmentIndex]->start;
			}
		}
		else
		{
			newVertex = ReadList[fragmentIndex]->clone();
			newVertex->start_min = newVertex->start;
			newVertex->end_max = newVertex->end;
			if (ReadList_cluster.size() >= ReadList_cluster.capacity())
				ReadList_cluster.reserve(default_dataset_num + ReadList_cluster.capacity());
			ReadList_cluster.push_back(newVertex);
			ReadList_cluster[ReadList_cluster.size() - 1]->includedReadNum = 1;

//			outputfile << ReadList[fragmentIndex]->readname << "\t" << ReadList_cluster.size() - 1 << "\t" << ReadList[fragmentIndex]->start << "-" << ReadList[fragmentIndex]->end << endl;
		}
	}

//	outputfile.close();

	return;
}


// Get the indices of the fragment clusters that need special process
ReadIndex* GetSpecialEndInidces(rangeJunction *transEndfraglist)
{
	rangeJunction *curfrag;
	long index, iLoop, specialEndNum, vertex1_Id, vertex2_Id;
	ReadIndex *transEndIndices, *newindex, *tailindex;
	vector<long> *indices;
	bool notspecial;
	transEndIndices = NULL;

	// get the number of alternative ends need special process
	specialEndNum = 0;
	curfrag = transEndfraglist->next;
	while(curfrag != NULL)
	{
		specialEndNum++;
		curfrag = curfrag->next;
	}
	indices = new vector<long>[specialEndNum];

	for (index = 0; index < ReadList_cluster.size(); index++)
	{
		iLoop = 0;
		curfrag = transEndfraglist->next;
		while(curfrag != NULL)
		{
			if (ReadList_cluster[index]->end_max >= curfrag->junc->start && ReadList_cluster[index]->end_max <= curfrag->junc->end)
			{
				if (indices[iLoop].size() >= indices[iLoop].capacity())
					indices[iLoop].reserve(default_dataset_num + indices[iLoop].capacity());
				indices[iLoop].push_back(index);
				break;
			}
			iLoop++;
			curfrag = curfrag->next;
		}
	}
	for (index = 0; index < specialEndNum; index++)
	{
		for (iLoop = 0; iLoop < indices[index].size(); iLoop++)
		{
			vertex1_Id = indices[index][iLoop];
			notspecial = false;
			for (long jLoop = indices[index].size() - 1; jLoop > iLoop; jLoop--)
			{
				vertex2_Id = indices[index][jLoop];
				if (ReadList_cluster[vertex2_Id]->start > ReadList_cluster[vertex1_Id]->end)
				{
					notspecial = true;
					break;
				}
			}
			if (notspecial == false)
			{
				newindex = new ReadIndex;
				newindex->readIndex = indices[index][iLoop];
				if (transEndIndices == NULL)
				{
					transEndIndices = newindex;
					tailindex = newindex;
				}
				else
				{
					tailindex->next = newindex;
					tailindex = newindex;
				}
			}
		}
	}


	for (index = 0; index < specialEndNum; index++)
	{
		indices[index].clear();
	}
	delete [] indices;

	return transEndIndices;
}


// cluster the similar fragments into smaller groups and build representative for each group
void ClusterReadFragments(GTvertex *vertex, int tolerance)
{
	long pointer, pointer_end, index;

	long TotalReadNum = ReadList.size();
	void** sortlist_read = new void* [TotalReadNum + 2];
	double* sortkey_read = new double [TotalReadNum + 2];

	// 1. sort the read fragments by the end points
	for (index = 1; index <= TotalReadNum; ++index)
	{
		sortlist_read[index] = (void*)ReadList[index - 1];
		sortkey_read[index] = ReadList[index - 1]->end;
	}
	mergeSort_general(sortlist_read, sortkey_read, TotalReadNum);
	for (index = 1; index <= TotalReadNum; ++index)
	{
		ReadList[index - 1] = (ReadVertex*) sortlist_read[index];
	}

	// 2. sort the read fragments by the start points
	for (index = 1; index <= TotalReadNum; ++index)
	{
		sortlist_read[index] = (void*)ReadList[index - 1];
		sortkey_read[index] = ReadList[index - 1]->start;
	}
	mergeSort_general(sortlist_read, sortkey_read, TotalReadNum);

	for (index = 1; index <= TotalReadNum; ++index)
	{
		ReadList[index - 1] = (ReadVertex*) sortlist_read[index];
	}

	delete [] sortlist_read;
	delete [] sortkey_read;

	// 3. linearly cluster the read fragments
	pointer = 0;
	while(pointer <= TotalReadNum - 1)
	{
		for (index = pointer + 1; index < TotalReadNum; index++)
		{
			if (abs(ReadList[index]->end - ReadList[pointer]->end) > tolerance || abs(ReadList[index]->start - ReadList[pointer]->start) > tolerance)
			{
				break;
			}
		}
		pointer_end = index - 1;

		// cluster read fragment from pointer to pointer_end
		ClusterReadFragments_subgroup(vertex, pointer, pointer_end, tolerance);

		pointer = pointer_end + 1;
	}

// 	ofstream outputfile;
// 	outputfile.open("yan_clusters.txt");
// 	for (index = 0; index < ReadList_cluster.size(); index++)
// 	{
// 		outputfile << index+1 << "\t" << ReadList_cluster[index]->readname << "\t" << ReadList_cluster[index]->includedReadNum << "\t" << ReadList_cluster[index]->start_min << "-" << ReadList_cluster[index]->end_max << endl;
// 	}
// 	outputfile.close();

//	transEndIndices = GetSpecialEndInidces(transEndfraglist);

	return;
}


long GetMinGraphpathLen(graphpath *paths, ReadVertex *vertex1, ReadVertex *vertex2, long &pathIndex)
{
	long length = MAX, temp, index = 0;
	graphpath *curpath;

	if (vertex1->end_max > vertex2->start_min)
	{
		length = 0;
	}
	else
	{
		curpath = paths;
		while(curpath != NULL)
		{
			if (curpath->edgelist == NULL)
			{
				temp = vertex2->start_min - vertex1->end_max;
			}
			else
			{
				temp = curpath->pathlength + vertex1->listTail->junc->end - vertex1->end_max + 1 + vertex2->start_min - vertex2->list->junc->start;
			}
			if (temp < length)
			{
				length = temp;
				pathIndex = index;
			}
			index++;
			curpath = curpath->next;
		}
	}

	if (length < 0)
	{
		length = 0;
	}
	return length;
}

// Construct graph on each gene 
void construct_splice_graph(RangeJunctionList *origList)
{
	//construct fragment graph based on given junctions

	JuncGraphVertex *curVertex, *tailVertex, *edgeVertex;
	JuncGraphEdge *newEdge;
	rangeJunction *curJunc;

	//build vertices

	curJunc = origList->list;
	while (curJunc != NULL)
	{
		curVertex = new JuncGraphVertex;
		curVertex->corresJunc = curJunc;
		curJunc = curJunc->next;
		curVertex->corresJunc->next = NULL;

		if (spliceGraph->vertices == NULL)
		{
			spliceGraph->vertices = curVertex;
			tailVertex = curVertex;
		} 
		else
		{
			tailVertex->next = curVertex;
			tailVertex = curVertex;
		}
	}

	//build edges

	curVertex = spliceGraph->vertices;
	while (curVertex != NULL)
	{
		edgeVertex = curVertex->next;
		while (edgeVertex != NULL)
		{
			if (compatibleJunctions((curVertex->corresJunc), (edgeVertex->corresJunc)) == true)
			{
				//create an edge for curVertex
				newEdge = new JuncGraphEdge;
				newEdge->linkedVertex = edgeVertex;

				newEdge->next = curVertex->edges;
				curVertex->edges = newEdge;

				edgeVertex->hasInEdge = true;
				curVertex->hasOutEdge = true;
			} 
			else
			{
				//break; 
			}
			edgeVertex = edgeVertex->next;
		}

		curVertex = curVertex->next;
	}

	return;
}


graphpath* generate_graphpath(rangeJunction *fragstart, rangeJunction *fragend, bool transcriptionEnds)
{
	//compute possible paths from fstart to fend
	//return the path list

	// 	RangeJunctionList *backupList;
	// 	backupList = origList->clone();
	// 	spliceGraph = new JuncGraph;
	// 	construct_splice_graph(backupList);
	// 	backupList->list = NULL;

	// 	if (fragstart->junc->start == 19952881 && fragstart->junc->end == 19956316 && fragend->junc->start == 19983359 && fragend->junc->end == 19984950)
	// 	{
	// 		cout << endl;
	// 	}
	graphpath *cur_path = NULL, *result_path = NULL, *new_path = NULL;
	JuncGraphVertex *fstart, *fend, *curVertex;
	JuncGraphEdge *cur_edge, *new_edge;
	long cnt_all_path_in_queue = 0, cnt_result_path = 0;
	int newChangeOfDir = 0;
	pathCluster *cluster;
	rangeJunction *newfrag;

	cluster = new pathCluster;
	newfrag = new rangeJunction;
	newfrag->junc = fragstart->junc->clone();
	cluster->fragstart = newfrag;
	newfrag = new rangeJunction;
	newfrag->junc = fragend->junc->clone();
	cluster->fragend = newfrag;

	result_path = NULL;
	queue<graphpath*> pathqueue;

	// find the corresponding vertices of the start/end fragments
	curVertex = spliceGraph->vertices;
	while(curVertex != NULL)
	{
		if (curVertex->corresJunc->junc->start == fragstart->junc->start && curVertex->corresJunc->junc->end == fragstart->junc->end)
		{
			fstart = curVertex;
		}
		if(curVertex->corresJunc->junc->start == fragend->junc->start && curVertex->corresJunc->junc->end == fragend->junc->end)
		{
			fend = curVertex;
		}

		curVertex = curVertex->next;
	}

	//build initial queue
	new_path = new graphpath;
	new_path->vertex_start = fstart; 
	new_path->vertex_end = fstart;
	new_path->direction = fstart->corresJunc->junc->transDirection;

	pathqueue.push(new_path);

	//main part
	if (transcriptionEnds == false)
	{
		while(pathqueue.empty() != true)
		{
			cur_path = pathqueue.front();
			pathqueue.pop();

			if (cur_path->pathlength <= gapConstraint)
			{
				if (cur_path->vertex_end == fend)
				{
					//currentPath achieves the destination
					++cnt_result_path;

					cur_path->next = result_path;
					result_path = cur_path;
					result_path->pathindex = cnt_result_path;
				} 
				else
				{
					cur_edge = cur_path->vertex_end->edges;
					while (cur_edge != NULL)
					{
						if (cur_path->edgelist_tail == NULL || cur_path->edgelist_tail->linkedVertex->corresJunc->junc->type != frag_junction || cur_edge->linkedVertex->corresJunc->junc->type != frag_junction)
						{
							//first, do not allow two junctions in a row

							if (cur_path->direction == undetermined || cur_edge->linkedVertex->corresJunc->junc->transDirection == undetermined || cur_path->direction == cur_edge->linkedVertex->corresJunc->junc->transDirection)
							{
								//second, must have consistent direction
								if (cur_edge->linkedVertex->corresJunc->junc->type != frag_junction && cur_edge->linkedVertex->corresJunc->junc->start > fragend->junc->start)
								{
									break;
								}
								new_path = cur_path->clone();

								if (new_path->edgelist_tail != NULL && (new_path->edgelist_tail->linkedVertex->corresJunc->junc->type == frag_exon || new_path->edgelist_tail->linkedVertex->corresJunc->junc->type == frag_retained_intron))
									new_path->pathlength += new_path->edgelist_tail->linkedVertex->corresJunc->junc->end - new_path->edgelist_tail->linkedVertex->corresJunc->junc->start + 1;

								new_edge = new JuncGraphEdge;
								new_edge->linkedVertex = cur_edge->linkedVertex;
								new_path->add_edge(new_edge);
								if (new_path->direction == undetermined && (new_edge->linkedVertex->corresJunc->junc->transDirection == sense || new_edge->linkedVertex->corresJunc->junc->transDirection == antisense))
									new_path->direction = new_edge->linkedVertex->corresJunc->junc->transDirection;

								new_path->vertex_end = cur_edge->linkedVertex;

								pathqueue.push(new_path);
							}
						}

						cur_edge = cur_edge->next;
					}

					delete cur_path;
				}
			}
			else
				delete cur_path;
		}
	}
	else
	{
		while(pathqueue.empty() != true)
		{
			cur_path = pathqueue.front();
			pathqueue.pop();

			if (cur_path->vertex_end == fend)
			{
				//currentPath achieves the destination
				++cnt_result_path;

				cur_path->next = result_path;
				result_path = cur_path;
				result_path->pathindex = cnt_result_path;
			} 
			else
			{
				cur_edge = cur_path->vertex_end->edges;
				while (cur_edge != NULL)
				{
					if (cur_path->edgelist_tail == NULL || cur_path->edgelist_tail->linkedVertex->corresJunc->junc->type != frag_junction || cur_edge->linkedVertex->corresJunc->junc->type != frag_junction)
					{
						//first, do not allow two junctions in a row

						if (cur_path->direction == undetermined || cur_edge->linkedVertex->corresJunc->junc->transDirection == undetermined || cur_path->direction == cur_edge->linkedVertex->corresJunc->junc->transDirection)
						{
							//second, must have consistent direction
							if (cur_edge->linkedVertex->corresJunc->junc->type != frag_junction && cur_edge->linkedVertex->corresJunc->junc->start > fragend->junc->start)
							{
								break;
							}
							new_path = cur_path->clone();

							if (new_path->edgelist_tail != NULL && (new_path->edgelist_tail->linkedVertex->corresJunc->junc->type == frag_exon || new_path->edgelist_tail->linkedVertex->corresJunc->junc->type == frag_retained_intron))
								new_path->pathlength += new_path->edgelist_tail->linkedVertex->corresJunc->junc->end - new_path->edgelist_tail->linkedVertex->corresJunc->junc->start + 1;

							new_edge = new JuncGraphEdge;
							new_edge->linkedVertex = cur_edge->linkedVertex;
							new_path->add_edge(new_edge);
							if (new_path->direction == undetermined && (new_edge->linkedVertex->corresJunc->junc->transDirection == sense || new_edge->linkedVertex->corresJunc->junc->transDirection == antisense))
								new_path->direction = new_edge->linkedVertex->corresJunc->junc->transDirection;

							new_path->vertex_end = cur_edge->linkedVertex;

							pathqueue.push(new_path);
						}
					}

					cur_edge = cur_edge->next;
				}

				delete cur_path;
			}
		}
	}

	//	delete spliceGraph;
	//	delete backupList;

	cluster->paths = result_path;

	if (FragpathCluster.size() >= FragpathCluster.capacity())
		FragpathCluster.reserve(default_dataset_num + FragpathCluster.capacity());
	FragpathCluster.push_back(cluster);

	return result_path;
}

// return valid paths between two vertices
graphpath* ValidatePaths(ReadVertex *vertex1, ReadVertex *vertex2, graphpath *paths)
{
	graphpath *resultpaths, *newpath, *curpath;
	JuncGraphEdge *curedge;
	RangeJunctionList *list, *reflist;
	rangeJunction *newfrag, *curfrag, *fragtail;
	int relation;
	bool found;
	resultpaths = NULL;

	curpath = paths;
	while(curpath != NULL)
	{
		list = new RangeJunctionList;
		list->list = NULL;
		curfrag = vertex1->list;
		while(curfrag != NULL)
		{
			newfrag = new rangeJunction;
			newfrag->junc = curfrag->junc->clone();
			if (list->list == NULL)
			{
				list->list = newfrag;
				fragtail = newfrag;
			}
			else
			{
				fragtail->next = newfrag;
				fragtail = newfrag;
			}
			curfrag = curfrag->next;
		}
		curedge = curpath->edgelist;
		if (curedge != NULL)
		{
			while(curedge->next != NULL)
			{
				newfrag = new rangeJunction;
				newfrag->junc = curedge->linkedVertex->corresJunc->junc->clone();
				if (list->list == NULL)
				{
					list->list = newfrag;
					fragtail = newfrag;
				}
				else
				{
					fragtail->next = newfrag;
					fragtail = newfrag;
				}
				curedge = curedge->next;
			}
			curfrag = vertex2->list;
			while(curfrag != NULL)
			{
				newfrag = new rangeJunction;
				newfrag->junc = curfrag->junc->clone();
				if (list->list == NULL)
				{
					list->list = newfrag;
					fragtail = newfrag;
				}
				else
				{
					fragtail->next = newfrag;
					fragtail = newfrag;
				}
				curfrag = curfrag->next;
			}
		}
		else
		{
			curfrag = vertex2->list->next;
			while(curfrag != NULL)
			{
				newfrag = new rangeJunction;
				newfrag->junc = curfrag->junc->clone();
				if (list->list == NULL)
				{
					list->list = newfrag;
					fragtail = newfrag;
				}
				else
				{
					fragtail->next = newfrag;
					fragtail = newfrag;
				}
				curfrag = curfrag->next;
			}
		}

		found = false;
		MergeAdjExons(list);
		for (long transIndex = 0; transIndex < TranscriptNum; transIndex++)
		{
			reflist = TransVec[transIndex]->transcript->clone();
			MergeAdjExons(reflist);
			relation = IdenticalTrans(reflist, list);
			if (relation == 1 || relation == 5)
			{
				found = true;
				break;
			}
			delete reflist;
		}
		delete list;
		if (found == true)
		{
			newpath = curpath->clone();
			newpath->next = resultpaths;
			resultpaths = newpath;
		}
		curpath = curpath->next;
	}

	return resultpaths;
}


void fillIngaps(connectedRead *cluster, graphpath *paths, ReadVertex *vertex1, ReadVertex *vertex2, bool transcriptionEnds, long *lengthlist)
{
	graphpath *curpath;
	pathlength_record *newRecord, *tailRecord;

	long length_min, length_max, length;
	double temp, gap = MAX/4;

	int i;
	if (vertex1->end_max > vertex2->start_min)
	{
		// could only occur for internal nodes
		length_min = 0;
		curpath = paths;
		while(curpath != NULL)
		{
			newRecord = new pathlength_record;
			if (curpath->edgelist == NULL)
			{
				length_max = vertex2->start - vertex1->end;
			}
			else
			{
				length_max = curpath->pathlength + vertex1->listTail->junc->end - vertex1->end + 1 + vertex2->start - vertex2->list->junc->start;
			}

			length_max = length_max > gapConstraint ? gapConstraint : length_max;
			for (length = length_min; length <= length_max; length++)
			{
				temp = -MixNormal(length);
				if (temp < gap)
				{
					lengthlist[length]++;
				}
			}
			newRecord->lowerbound = length_min;
			newRecord->higherbound = length_max;
			if (cluster->lengthlist == NULL)
			{
				cluster->lengthlist = newRecord;
				tailRecord = newRecord;
			}
			else
			{
				tailRecord->next = newRecord;
				tailRecord = newRecord;
			}

			curpath = curpath->next;
		}
	}
	else
	{
		curpath = paths;
		while(curpath != NULL)
		{
			newRecord = new pathlength_record;
			if (curpath->edgelist == NULL)
			{
				length_min = vertex2->start_min - vertex1->end_max;
				length_max = vertex2->start - vertex1->end;
				if (transcriptionEnds == true)
				{
					if (abs(length_min - 1 ) <= 1)
					{
						newRecord->lowerbound = 0;
						newRecord->higherbound = 0;
						newRecord->next = cluster->lengthlist;
						cluster->lengthlist = newRecord;
						break;
					}
				}
			}
			else
			{
				length_min = curpath->pathlength + vertex1->listTail->junc->end - vertex1->end_max + 1 + vertex2->start_min - vertex2->list->junc->start;
				length_max = curpath->pathlength + vertex1->listTail->junc->end - vertex1->end + 1 + vertex2->start - vertex2->list->junc->start;
			}
			if (length_min <= gapConstraint)
			{
				length_max = length_max > gapConstraint ? gapConstraint : length_max;
				for (length = length_min; length <= min(length_max, gapConstraint); length++)
				{
					temp = -MixNormal(length);
					if (temp < gap)
					{
						lengthlist[length]++;
					}
				}
			}
			newRecord->lowerbound = length_min;
			newRecord->higherbound = length_max;
			if (cluster->lengthlist == NULL)
			{
				cluster->lengthlist = newRecord;
				tailRecord = newRecord;
			}
			else
			{
				tailRecord->next = newRecord;
				tailRecord = newRecord;
			}
			curpath = curpath->next;
		}
	}

	return;
}

double bandWidth(long* lengthlist, long vectorsize)
{
	// Silverman Rule-of-Thumb
	// Epanechnikov kernel smoothing is adopted
	// It is a second order kernel

	double sd, mean, constant, h;
	long index, totalCnt;

	mean = 0;
	totalCnt = 0;
	if (vectorsize > 0)
	{
		for (index = 0; index < vectorsize; index++)
		{
			mean += lengthlist[index] * index;
			totalCnt += lengthlist[index];
		}
	}

	if (totalCnt > 1)
	{
		mean = mean/totalCnt;
		sd = 0;
		for (index = 0; index < vectorsize; index++)
		{
			sd += pow(((double)index - mean),2) * (double) lengthlist[index];
		}
		sd = sd/(totalCnt - 1);
		sd = sqrt(sd);

		constant = 2.34;
		h = sd * constant * (1/pow(totalCnt, 1.0/3.0));
	}
	else
	{
		sd = 0;
		h = 0;
	}

	return h;
}


double EmpiricalDist(long x, double h, long* lengthlist, long vectorsize)
{
	// Get the empirical distribution of the gaps
	// Epanechnikov kernel smoothing is adopted

	/*	sort(gapList.begin(), gapList.end());*/
	double pdf;
	long index, totalCnt;

	pdf = 0;
	totalCnt = 0;
	for (index = 0; index < vectorsize; index++)
	{
		if (abs((double)(index - x)/h) <= 1)
		{
			//			cout << (double)(gapList[index] - x)/h << endl;
			pdf += 3.0/4.0 * (1 - pow(((double)(index - x)/h), 2)) * (double) lengthlist[index];
		}
		totalCnt += lengthlist[index];
	}
	pdf = pdf/(h*totalCnt);

	pdf = log(pdf);
	if (pdf < -MAX/4)
	{
		pdf = -MAX/4;
	}
	return pdf;
}

/************************************************************************/
/* Date: 09/21/2013 */
// function control false positives
/************************************************************************/

// get the covered frags from two reads and the connecting sequences
ReadInfo* GetCoveredFrags_extension(GTvertex *genevertex, ReadVertex *vertex1, ReadVertex *vertex2, graphpath *connectedpaths)
{
	rangeJunction *curfrag = genevertex->junctionInRange->list, *readfrag1 = vertex1->list, *readfrag2 = vertex2->list->next, *tempfrag;
	ReadInfo *readinfoList = NULL, *newreadinfo, *curreadinfo, *tailinfo;
	graphpath *curpath;
	JuncGraphEdge *curedge;
	vector<int> prevreadinfo, postreadinfo;

	int fragindex = 1, tempindex;
	while(curfrag != NULL && readfrag1 != NULL)
	{
		if (curfrag->junc->start == readfrag1->junc->start && curfrag->junc->end == readfrag1->junc->end)
		{
			prevreadinfo.push_back(fragindex);
			readfrag1 = readfrag1->next;
		}
		fragindex++;
		curfrag = curfrag->next;
	}
	tempfrag = curfrag;
	tempindex = fragindex;

	if (connectedpaths != NULL)
	{
		curpath = connectedpaths;
		while(curpath != NULL)
		{
			newreadinfo = new ReadInfo;
			for (int index = 0; index < prevreadinfo.size(); index++)
			{
				newreadinfo->FragIndexCovered.push_back(prevreadinfo[index]);
			}
			if (curpath->edgelist == NULL)
			{
				// do nothing
			}
			else
			{
				curedge = curpath->edgelist;
				curfrag = tempfrag;
				fragindex = tempindex;
				while(curfrag != NULL && curedge != NULL)
				{
					if (curfrag->junc->start == curedge->linkedVertex->corresJunc->junc->start && curfrag->junc->end == curedge->linkedVertex->corresJunc->junc->end)
					{
						newreadinfo->FragIndexCovered.push_back(fragindex);
						curedge = curedge->next;
					}
					curfrag = curfrag->next;
					fragindex++;
				}
			}
			if (readinfoList == NULL)
			{
				readinfoList = newreadinfo;
				tailinfo = newreadinfo;
			}
			else
			{
				tailinfo->nextreadInfo = newreadinfo;
				tailinfo = newreadinfo;
			}

			curpath = curpath->next;
		}
	}

	curfrag = tempfrag;
	fragindex = tempindex;
	while(curfrag != NULL && readfrag2 != NULL)
	{
		if (curfrag->junc->start == readfrag2->junc->start && curfrag->junc->end == readfrag2->junc->end)
		{
			postreadinfo.push_back(fragindex);
			readfrag2 = readfrag2->next;
		}
		fragindex++;
		curfrag = curfrag->next;
	}

	curreadinfo = readinfoList;
	while(curreadinfo != NULL)
	{
		for (int index = 0; index < postreadinfo.size(); index++)
		{
			curreadinfo->FragIndexCovered.push_back(postreadinfo[index]);
		}
		curreadinfo = curreadinfo->nextreadInfo;
	}

	prevreadinfo.clear();
	postreadinfo.clear();

	return readinfoList;
}

// Calculate the likelihood that the transcript piece constituted by two vertices is true
bool isDeterministicPath(ReadInfo *readinfo_ref, ReadInfo *readinfo)
{
	bool isDeterministic = true, overlap = false;
	int index_ref, index;
	if (readinfo_ref->FragIndexCovered.size() > readinfo->FragIndexCovered.size())
	{
		// fragments in readinfo should be totally contained in readinfo_ref
		index_ref = 0;
		while(index_ref < readinfo_ref->FragIndexCovered.size())
		{
			if (readinfo_ref->FragIndexCovered[index_ref] == readinfo->FragIndexCovered[0])
			{
				overlap = true;
				break;
			}
			index_ref++;
		}
		if (overlap == true)
		{
			index = 0;
			while(index_ref < readinfo_ref->FragIndexCovered.size() && index < readinfo->FragIndexCovered.size())
			{
				if (readinfo_ref->FragIndexCovered[index_ref] != readinfo->FragIndexCovered[index])
				{
					isDeterministic = false;
					break;
				}
				index++;
				index_ref++;
			}
			if (index_ref == readinfo_ref->FragIndexCovered.size() && index < readinfo->FragIndexCovered.size())
			{
				isDeterministic == false;
			}
		}
		else
			isDeterministic = false;
	}
	else if (readinfo_ref->FragIndexCovered.size() == readinfo->FragIndexCovered.size())
	{
		// fragments should be exactly the same
		for (index = 0; index < readinfo->FragIndexCovered.size(); index++)
		{
			if (readinfo->FragIndexCovered[index] != readinfo_ref->FragIndexCovered[index])
			{
				isDeterministic = false;
				break;
			}
		}
	}
	else
	{
		// fragments in readinfo_ref should be totally contained in readinfo
		index = 0;
		while(index < readinfo->FragIndexCovered.size())
		{
			if (readinfo->FragIndexCovered[index] == readinfo_ref->FragIndexCovered[0])
			{
				overlap = true;
				break;
			}
			index++;
		}
		if (overlap == true)
		{
			index_ref = 0;
			while(index_ref < readinfo_ref->FragIndexCovered.size() && index < readinfo->FragIndexCovered.size())
			{
				if (readinfo_ref->FragIndexCovered[index_ref] != readinfo->FragIndexCovered[index])
				{
					isDeterministic = false;
					break;
				}
				index++;
				index_ref++;
			}
			if (index == readinfo->FragIndexCovered.size() && index_ref < readinfo_ref->FragIndexCovered.size())
			{
				isDeterministic = false;
			}
		}
		else
			isDeterministic = false;
	}
	
	return isDeterministic;
}

int IdenticalTrans_ASM(vector <int> ref, vector <int> readbuilt)
{
	int relation = -1; 
	// 0 indicates trans 2 is different from trans 1; 1 indicates trans1 contains trans2; 2 indicates trans1 is contained within trans1
	// 3 indicates overlapping (trans1 is downstream); 4 indicates overlapping (trans2 is downstream); 5 indicates identical.

	bool overlap = false;
	int index_ref, index;
	index_ref = index = 0;
	if (ref[0] < readbuilt[0])
	{
		while(ref[index_ref] < readbuilt[0])
		{
			if (ref[index_ref] > readbuilt[0])
			{
				break;
			}
			else if (ref[index_ref] == readbuilt[0])
			{
				overlap = true;
				break;
			}
			index_ref++;
		}
	}
	else if (ref[0] > readbuilt[0])
	{
		while(readbuilt[index] < ref[0])
		{
			if (readbuilt[index] > ref[0])
			{
				break;
			}
			else if (readbuilt[index] == ref[0])
			{
				overlap = true;
				break;
			}
			index++;
		}
	}
	else
	{
		overlap = true;
	}

	if (overlap == true)
	{
		while(index_ref < ref.size() && index < readbuilt.size())
		{
			if (ref[index_ref] != readbuilt[index])
			{
				relation = 0;
				break;
			}
			index++;
			index_ref++;
		}
		if (relation != 0)
		{
			if (index_ref < ref.size())
			{
				if (ref[0] < readbuilt[0])
				{
					relation = 3;
				}
				else
					relation = 2;
			}
			else if (index < readbuilt.size())
			{
				if (readbuilt[0] < ref[0])
				{
					relation = 4;
				}
				else
					relation = 1;
			}
			else
			{
				if (readbuilt[0] < ref[0])
				{
					relation = 2;
				}
				else if (readbuilt[0] > ref[0])
				{
					relation = 1;
				}
				else
				{
					relation = 5;
				}
			}
		}
	}
	else
	{
		relation = 0;
	}

	return relation;
}


double CalculateTransProb_ASMlevel(ReadInfo *readinfo, TransInfo *transInfo_ASM, long meanReadLength, long HubLength)
{
	int relation;
	double prob, Optimalprob = 0;
	TransInfo *curTransinfo = transInfo_ASM;
	while(curTransinfo != NULL)
	{
		relation = IdenticalTrans_ASM(curTransinfo->FragIndexCovered, readinfo->FragIndexCovered);
		if (relation != 0 && relation != -1)
		{
			if (relation == 0 || relation == 1 || relation == 2)
			{
				Optimalprob = 1-double(meanReadLength - HubLength - 1)/(curTransinfo->transLength - meanReadLength + 1);
				break;
			}
			else
			{
				prob = 1-double(meanReadLength - HubLength - 1)/(curTransinfo->transLength - meanReadLength + 1);
				if (prob > Optimalprob)
				{
					Optimalprob = prob;
				}
			}
			
		}
		curTransinfo = curTransinfo->nextTransinfo;
	}

	return Optimalprob;
}


void TransExistProb(GTvertex *genevertex, AdjASMsInfo *adjASMinfoList, ReadVertex *vertex1, ReadVertex *vertex2, graphpath *connectedpaths)
{
	double prob;
	AdjASMsInfo *curASMsinfo = adjASMinfoList;
	ReadInfo* readinfoList, *delinfo, *curinfo, *curinfo_ref, *transinfo_read;
	bool builtList = false, isDeterministic;
	int transindex;

	TransProb_ASM.clear();
	readinfoList = NULL;
	while(curASMsinfo != NULL)
	{
		if (vertex1->start_min < curASMsinfo->Hubstart && vertex2->end_max > curASMsinfo->Hubend)
		{
			// the genomic region spanned by two vertices cover this hub
			// 1. construct the fragment lists covered by the two vertices
			if (builtList == false)
			{
				readinfoList = GetCoveredFrags_extension(genevertex, vertex1, vertex2, connectedpaths);
				builtList = true;
				curinfo = readinfoList;
				while(curinfo != NULL)
				{
					TransProb_ASM.push_back(1);
					curinfo = curinfo->nextreadInfo;
				}
			}
			
			// 2. calculate the probability on these two adjacent ASMs
			transindex = 0;
			curinfo = readinfoList;
			while(curinfo != NULL)
			{
				isDeterministic = false;
				curinfo_ref = curASMsinfo->deterministicReadinfo;
				while(curinfo_ref != NULL)
				{
					if (isDeterministicPath(curinfo_ref, curinfo) == true)
					{
						// do nothing
						isDeterministic = true;
						break;
					}
					curinfo_ref = curinfo_ref->nextreadInfo;
				}
				if (isDeterministic == false)
				{
					// calculate the probability
					prob = CalculateTransProb_ASMlevel(curinfo, curASMsinfo->ASM_Transinfo, MeanReadLength, curASMsinfo->Hublength);
					TransProb_ASM[transindex] *= prob;

				}
				transindex++;
				curinfo = curinfo->nextreadInfo;
			}
			
		}
		curASMsinfo = curASMsinfo->nextinfo;
	}

	while(readinfoList != NULL)
	{
		delinfo = readinfoList;
		readinfoList = delinfo->nextreadInfo;
		delete delinfo;
	}
	return;
}

double getOptimalProb(GTvertex *genevertex, AdjASMsInfo *adjASMinfoList, connectedRead *cluster, void** sortlist, double* sortkey, long dimension, graphpath *connectedpaths)
{
	long index, transindex;
	double optimalProb = -MAX, prob, temp;
	ReadVertex *vertex1, *vertex2;
	pathlength_record *record;
	graphpath *curpath;
	vertex1 = ReadList_cluster[cluster->readlowIndex];
	vertex2 = ReadList_cluster[cluster->readhighIndex];

	// 1. check whether the edge weight between these two vertices need to be adjusted
	TransExistProb(genevertex, adjASMinfoList, vertex1, vertex2, connectedpaths);
	if (TransProb_ASM.size() > 0)
	{
		// need adjust (span two ASMs)
		optimalProb = -MAX;
		record = cluster->lengthlist;
		curpath = connectedpaths;
		transindex = 0;
		while(record != NULL)
		{
			if (TransProb_ASM[transindex] == 1)
			{
				optimalProb = sortkey[dimension - 1] + constant_Penalty;
				cluster->pathselectedIndex = curpath->pathindex;
				break;
			}
			else
			{
				if (record->lowerbound <= gapConstraint)
				{
					prob = PenaltyFunc(TransProb_ASM[transindex], ReadNum, TransNum);
					if (prob < -MAX/4)
					{
						prob = -MAX/4;
					}
					for (index = dimension; index >= 1; index--)
					{
						if ((long)sortlist[index] >= record->lowerbound && (long)sortlist[index] <= record->higherbound)
						{
							temp = sortkey[index];
							break;
						}
					}
					temp += prob;
					if (temp > optimalProb)
					{
						optimalProb = temp;
						cluster->pathselectedIndex = curpath->pathindex;
					}
				}
			}
			transindex++;
			record = record->next;
			curpath = curpath->next;
		}
	}
	else
	{
		// do not need adjust
		for (index = dimension - 1; index >= 1; index--)
		{
			record = cluster->lengthlist;
			curpath = connectedpaths;
			while(record != NULL)
			{
				if ((long)sortlist[index] >= record->lowerbound && (long)sortlist[index] <= record->higherbound)
				{
					optimalProb = sortkey[index];
					cluster->pathselectedIndex = curpath->pathindex;
					break;
				}
				curpath = connectedpaths->next;
				record = record->next;
			}
			
		}
	}
	TransProb_ASM.clear();

	return optimalProb;
}

long getOptimalPathIndex(graphpath *paths, ReadVertex *vertex1, ReadVertex *vertex2, long length, bool transcriptionEnds);

// Initialize flow graph capability and cost 
void InitializeFlowGraph(GTvertex *genevertex, long **cap, double **cost, vector< vector<int> > &connectedIndex, long source, long sink, long transStartNum, long transEndNum, long &maxSup, long totalReadNum, AdjASMsInfo *adjASMinfoList, long **IndexLocator_fragcluster, long **IndexLocator_connectedReadscluter)
{
	bool compatible, connected, pathsConstructed, found, transcriptionEnds;
	long rowIndex, colIndex, index, fragpathIndex, pathIndex, length, maxSup1, maxSup2, *lengthlist, dimension;
	graphpath *subpaths, *curpath, *paths;
	graphpathCluster *cluster;
	ReadVertex *vertex1, *vertex2;
	connectedRead *connectedreadcluster;
	pathlength_record *newRecord;

	dimension = gapConstraint + 1;
	lengthlist = new long[dimension];
	for (rowIndex = 0; rowIndex < dimension; rowIndex++)
	{
		lengthlist[rowIndex] = 0;
	}
	subpaths = NULL;

	cap[source][sink] = MAX;
	cost[source][sink] = MAX/4;

	// Build edges from source to transcription starts and transcription ends to sink
	for (rowIndex = 1; rowIndex <= transStartNum; rowIndex++)
	{
		cap[source][rowIndex] = totalReadNum; // the capacity equals the total number of reads
		cost[source][rowIndex] = 0;
	}
	rowIndex = sink - 1; 
	index = 0;
	while(index < transEndNum)
	{
		cap[rowIndex][sink] = totalReadNum;
		cost[rowIndex][sink] = 0;
		rowIndex--;
		index++;
	}

	// Build edges from virtual transcription start/end to other internal vertices
	transcriptionEnds = true;
	maxSup1 = 0;
	for (rowIndex = source + 1; rowIndex <= transStartNum; rowIndex++)
	{
		vertex1 = ReadList_cluster[ReadVertexIndex[rowIndex - 1]];
		for (colIndex = transStartNum + 1; colIndex < sink - transEndNum; colIndex = colIndex + 2)
		{
			vertex2 = ReadList_cluster[ReadVertexIndex[(colIndex - 1 - transStartNum)/2 + transStartNum]];
			pathsConstructed = false;
			for (fragpathIndex = 0; fragpathIndex < FragpathCluster.size(); fragpathIndex++)
			{
				if (FragpathCluster[fragpathIndex]->fragstart->junc->start == vertex1->list->junc->start && FragpathCluster[fragpathIndex]->fragstart->junc->end == vertex1->list->junc->end && FragpathCluster[fragpathIndex]->fragend->junc->start == vertex2->list->junc->start && FragpathCluster[fragpathIndex]->fragend->junc->end == vertex2->list->junc->end)
				{
					subpaths = FragpathCluster[fragpathIndex]->paths;
					pathsConstructed = true;
					IndexLocator_fragcluster[rowIndex - 1][(colIndex - 1 - transStartNum)/2 + transStartNum] = fragpathIndex;
					break;
				}
			}
			if (pathsConstructed == false)
			{
				subpaths = generate_graphpath(vertex1->list, vertex2->list, transcriptionEnds);
				IndexLocator_fragcluster[rowIndex - 1][(colIndex - 1 - transStartNum)/2 + transStartNum] = FragpathCluster.size() - 1;
			}

			if (subpaths != NULL)
			{
				length = GetMinGraphpathLen(subpaths, vertex1, vertex2, pathIndex);
				connectedreadcluster = new connectedRead;
				connectedreadcluster->readlowIndex = rowIndex - 1;
				connectedreadcluster->readhighIndex = (colIndex - 1 - transStartNum)/2 + transStartNum;
				connectedreadcluster->minLength_inbetween = length;

				if (length <= gapConstraint)
				{
					maxSup1 += vertex2->includedReadNum;
					fillIngaps(connectedreadcluster, subpaths, vertex1, vertex2, true, lengthlist);
				}
				else
				{					
					connectedreadcluster->pathselectedIndex = getOptimalPathIndex(subpaths, vertex1, vertex2, length, true);
					newRecord = new pathlength_record;
					newRecord->lowerbound = length;
					newRecord->higherbound = length;
					newRecord->next = connectedreadcluster->lengthlist;
					connectedreadcluster->lengthlist = newRecord;
				}
				if (ConnectedReadCluster.size() >= ConnectedReadCluster.capacity())
					ConnectedReadCluster.reserve(default_dataset_num + ConnectedReadCluster.capacity());
				ConnectedReadCluster.push_back(connectedreadcluster);
				IndexLocator_connectedReadscluter[rowIndex - 1][(colIndex - 1 - transStartNum)/2 + transStartNum] = ConnectedReadCluster.size() - 1;

				cap[rowIndex][colIndex] = vertex2->includedReadNum;				
			}	
		}
	}

	maxSup2 = 0;
	for (rowIndex = sink - 1; rowIndex >= sink - transEndNum; rowIndex--)
	{
		vertex1 = ReadList_cluster[ReadVertexIndex[transStartNum - 1 + (vertexNum - 2 - transStartNum - transEndNum)/2 + rowIndex - (sink - transEndNum) + 1]];
		for (colIndex = sink - transEndNum - 1; colIndex > transStartNum; colIndex = colIndex - 2)
		{
			vertex2 = ReadList_cluster[ReadVertexIndex[(colIndex - 1 - transStartNum)/2 + transStartNum]];
			pathsConstructed = false;
			for (fragpathIndex = 0; fragpathIndex < FragpathCluster.size(); fragpathIndex++)
			{
				if (FragpathCluster[fragpathIndex]->fragstart->junc->start == vertex2->listTail->junc->start && FragpathCluster[fragpathIndex]->fragstart->junc->end == vertex2->listTail->junc->end && FragpathCluster[fragpathIndex]->fragend->junc->start == vertex1->list->junc->start && FragpathCluster[fragpathIndex]->fragend->junc->end == vertex1->list->junc->end)
				{
					subpaths = FragpathCluster[fragpathIndex]->paths;
					pathsConstructed = true;
					IndexLocator_fragcluster[(colIndex - 1 - transStartNum)/2 + transStartNum][transStartNum - 1 + (vertexNum - 2 - transStartNum - transEndNum)/2 + rowIndex - (sink - transEndNum) + 1] = fragpathIndex;
					break;
				}
			}
			if (pathsConstructed == false)
			{
				subpaths = generate_graphpath(vertex2->listTail, vertex1->list, transcriptionEnds);
				IndexLocator_fragcluster[(colIndex - 1 - transStartNum)/2 + transStartNum][transStartNum - 1 + (vertexNum - 2 - transStartNum - transEndNum)/2 + rowIndex - (sink - transEndNum) + 1] = FragpathCluster.size() - 1;
			}
			if (subpaths != NULL)
			{
				connectedreadcluster = new connectedRead;				
				connectedreadcluster->readlowIndex = (colIndex - 1 - transStartNum)/2 + transStartNum;
				connectedreadcluster->readhighIndex = transStartNum - 1 + (vertexNum - 2 - transStartNum - transEndNum)/2 + rowIndex - (sink - transEndNum) + 1;

				length = GetMinGraphpathLen(subpaths, vertex2, vertex1, pathIndex);
				connectedreadcluster->minLength_inbetween = length;
				if (length <= gapConstraint)
				{
					maxSup2 += vertex2->includedReadNum;
					fillIngaps(connectedreadcluster,  subpaths, vertex2, vertex1, true, lengthlist);
				}
				else
				{				
					newRecord = new pathlength_record;
					newRecord->lowerbound = length;
					newRecord->higherbound = length;
					newRecord->next = connectedreadcluster->lengthlist;
					connectedreadcluster->lengthlist = newRecord;
					connectedreadcluster->pathselectedIndex = getOptimalPathIndex(subpaths, vertex2, vertex1, length, true);
				}
				if (ConnectedReadCluster.size() >= ConnectedReadCluster.capacity())
					ConnectedReadCluster.reserve(default_dataset_num + ConnectedReadCluster.capacity());
				ConnectedReadCluster.push_back(connectedreadcluster);
				IndexLocator_connectedReadscluter[(colIndex - 1 - transStartNum)/2 + transStartNum][transStartNum - 1 + (vertexNum - 2 - transStartNum - transEndNum)/2 + rowIndex - (sink - transEndNum) + 1] = ConnectedReadCluster.size() - 1;

				cap[colIndex][rowIndex] = vertex2->includedReadNum;				
			}	
		}
	}

	// Build edges among the internal vertices
	rangeJunction *frag1, *frag2;
	transcriptionEnds = false;
	for (rowIndex = transStartNum + 1; rowIndex < sink - transEndNum; rowIndex = rowIndex + 2)
	{
		vertex1 = ReadList_cluster[ReadVertexIndex[(rowIndex - 1 - transStartNum)/2 + transStartNum]];
		cost[rowIndex][rowIndex + 1] = epsilon;
		frag1 = vertex1->list;
		while(frag1 != NULL)
		{
			if (vertex1->end >= frag1->junc->start && vertex1->end <= frag1->junc->end && frag1->junc->type != frag_junction)
			{
				break;
			}
			frag1 = frag1->next;
		}
		for (index = rowIndex + 2; index < sink - transEndNum; index = index + 2)
		{
			vertex2 = ReadList_cluster[ReadVertexIndex[(index - 1 - transStartNum)/2 + transStartNum]];	
			if (vertex2->start - vertex1->end >= 0)
			{
				// find the two rangeJunctions we want to build paths in between
				frag2 = vertex2->list;
				while(frag2 != NULL)
				{
					if (vertex2->start >= frag2->junc->start && vertex2->start <= frag2->junc->end && frag2->junc->type != frag_junction)
					{
						break;
					}
					frag2 = frag2->next;
				}
				pathsConstructed = false;
				for (fragpathIndex = 0; fragpathIndex < FragpathCluster.size(); fragpathIndex++)
				{
					if (FragpathCluster[fragpathIndex]->fragstart->junc->start == frag1->junc->start && FragpathCluster[fragpathIndex]->fragstart->junc->end == frag1->junc->end && FragpathCluster[fragpathIndex]->fragend->junc->start == frag2->junc->start && FragpathCluster[fragpathIndex]->fragend->junc->end == frag2->junc->end)
					{
						subpaths = FragpathCluster[fragpathIndex]->paths;
						pathsConstructed = true;
						IndexLocator_fragcluster[(rowIndex - 1 - transStartNum)/2 + transStartNum][(index - 1 - transStartNum)/2 + transStartNum] = fragpathIndex;
						break;
					}
				}
				if (pathsConstructed == false)
				{
					subpaths = generate_graphpath(frag1, frag2, transcriptionEnds);
					IndexLocator_fragcluster[(rowIndex - 1 - transStartNum)/2 + transStartNum][(index - 1 - transStartNum)/2 + transStartNum] = FragpathCluster.size() - 1;
				}
				if (subpaths != NULL)
				{
					connectedreadcluster = new connectedRead;				
					connectedreadcluster->readlowIndex = (rowIndex - 1 - transStartNum)/2 + transStartNum;
					connectedreadcluster->readhighIndex = (index - 1 - transStartNum)/2 + transStartNum;

					length = GetMinGraphpathLen(subpaths, vertex1, vertex2, pathIndex);
					connectedreadcluster->minLength_inbetween = length;
					if (length < gapConstraint)
					{
						cap[rowIndex + 1][index] = min(vertex1->includedReadNum, vertex2->includedReadNum);
						fillIngaps(connectedreadcluster, subpaths, vertex1, vertex2, false, lengthlist);
					}
					else
					{
						newRecord = new pathlength_record;
						newRecord->lowerbound = length;
						newRecord->higherbound = length;
						newRecord->next = connectedreadcluster->lengthlist;
						connectedreadcluster->lengthlist = newRecord;
					}
					if (ConnectedReadCluster.size() >= ConnectedReadCluster.capacity())
						ConnectedReadCluster.reserve(default_dataset_num + ConnectedReadCluster.capacity());
					ConnectedReadCluster.push_back(connectedreadcluster);
					IndexLocator_connectedReadscluter[(rowIndex - 1 - transStartNum)/2 + transStartNum][(index - 1 - transStartNum)/2 + transStartNum] = ConnectedReadCluster.size() - 1;

				}
			}
		}

	}
	maxSup = max(maxSup1, maxSup2);

	double width = bandWidth(lengthlist, dimension);

	// speed up process
	double *LikelihoodDiff;
	LikelihoodDiff = new double[dimension];
	for (index = 0; index < dimension; index++)
	{
		if (width > 0)
		{
			LikelihoodDiff[index] = MixNormal(index) - EmpiricalDist(index, width, lengthlist, dimension);
		}
		else
		{
			LikelihoodDiff[index] = MixNormal(index);
		}
	}
	delete [] lengthlist;
	// sort LikelihoodDiff from largest to smallest
	void** sortlist;
	double* sortkey;
	sortlist = new void* [dimension + 2];
	sortkey = new double [dimension + 2];
	for (long tmpIndex = 0; tmpIndex < dimension; tmpIndex++)
	{
		sortlist[tmpIndex + 1] = (void*)tmpIndex;
		sortkey[tmpIndex + 1] = LikelihoodDiff[tmpIndex];
	}	
	mergeSort_general(sortlist, sortkey, dimension);

	double temp, gap = -MAX, maxgap = -MAX, **optimalweight;
	optimalweight = new double*[vertexNum];
	for (rowIndex = 0; rowIndex < vertexNum; rowIndex++)
	{
		optimalweight[rowIndex] = new double[vertexNum];
		for (colIndex = 0; colIndex < vertexNum; colIndex++)
		{
			optimalweight[rowIndex][colIndex] = -MAX;
		}
	}
	for (rowIndex = source + 1; rowIndex <= transStartNum; rowIndex++)
	{
		vertex1 = ReadList_cluster[ReadVertexIndex[rowIndex - 1]];
		for (colIndex = transStartNum + 1; colIndex < sink - transEndNum; colIndex = colIndex + 2)
		{
			vertex2 = ReadList_cluster[ReadVertexIndex[(colIndex - 1 - transStartNum)/2 + transStartNum]];
			gap = -MAX, length = MAX;
			// choose the one length with maximal weight
			subpaths = FragpathCluster[IndexLocator_fragcluster[rowIndex - 1][(colIndex - 1 - transStartNum)/2 + transStartNum]]->paths;
			if (subpaths != NULL)
			{
				connectedreadcluster = ConnectedReadCluster[IndexLocator_connectedReadscluter[rowIndex - 1][(colIndex - 1 - transStartNum)/2 + transStartNum]];
				if (connectedreadcluster->lengthlist->lowerbound == 0 && connectedreadcluster->lengthlist->higherbound == 0)
				{
					cost[rowIndex][colIndex] = epsilon;
					optimalweight[rowIndex][colIndex] = epsilon;
					length = 0;
					connectedreadcluster->pathselectedIndex = getOptimalPathIndex(subpaths, vertex1, vertex2, length, true);
				}
				else
				{
					if (connectedreadcluster->minLength_inbetween >= gapConstraint)
					{
						cost[rowIndex][colIndex] = connectedreadcluster->lengthlist->lowerbound;
					}
					else
					{
						gap = getOptimalProb(genevertex, adjASMinfoList, connectedreadcluster, sortlist, sortkey, dimension, subpaths);
						optimalweight[rowIndex][colIndex] = gap;
						if (gap > maxgap)
						{
							maxgap = gap;
						}
					}

				}
			}
		}
	}
	for (rowIndex = sink - 1; rowIndex >= sink - transEndNum; rowIndex--)
	{
		vertex1 = ReadList_cluster[ReadVertexIndex[transStartNum - 1 + (vertexNum - 2 - transStartNum - transEndNum)/2 + rowIndex - (sink - transEndNum) + 1]];
		for (colIndex = sink - transEndNum - 1; colIndex > transStartNum; colIndex = colIndex - 2)
		{
			vertex2 = ReadList_cluster[ReadVertexIndex[(colIndex - 1 - transStartNum)/2 + transStartNum]];
			gap = -MAX, length = MAX;
			subpaths = FragpathCluster[IndexLocator_fragcluster[(colIndex - 1 - transStartNum)/2 + transStartNum][transStartNum - 1 + (vertexNum - 2 - transStartNum - transEndNum)/2 + rowIndex - (sink - transEndNum) + 1]]->paths;
			if (subpaths != NULL)
			{
				connectedreadcluster = ConnectedReadCluster[IndexLocator_connectedReadscluter[(colIndex - 1 - transStartNum)/2 + transStartNum][transStartNum - 1 + (vertexNum - 2 - transStartNum - transEndNum)/2 + rowIndex - (sink - transEndNum) + 1]];
				if (connectedreadcluster->lengthlist->lowerbound == 0 && connectedreadcluster->lengthlist->higherbound == 0)
				{
					cost[colIndex][rowIndex] = epsilon;
					optimalweight[colIndex][rowIndex] = epsilon;
					length = 0;
					connectedreadcluster->pathselectedIndex = getOptimalPathIndex(subpaths, vertex2, vertex1, length, true);
				}
				else
				{
					if (connectedreadcluster->minLength_inbetween >= gapConstraint)
					{
						cost[colIndex][rowIndex] = connectedreadcluster->lengthlist->lowerbound;
					}
					else
					{
						gap = getOptimalProb(genevertex, adjASMinfoList, connectedreadcluster, sortlist, sortkey, dimension, subpaths);
						optimalweight[colIndex][rowIndex] = gap;
						if (gap > maxgap)
						{
							maxgap = gap;
						}
					}

				}
			}

		}
	}

	for (rowIndex = transStartNum + 1; rowIndex < sink - transEndNum; rowIndex = rowIndex + 2)
	{
		vertex1 = ReadList_cluster[ReadVertexIndex[(rowIndex - 1 - transStartNum)/2 + transStartNum]];
		frag1 = vertex1->list;
		while(frag1 != NULL)
		{
			if (vertex1->end >= frag1->junc->start && vertex1->end <= frag1->junc->end && frag1->junc->type != frag_junction)
			{
				break;
			}
			frag1 = frag1->next;
		}
		for (colIndex = rowIndex + 2; colIndex < sink - transEndNum; colIndex = colIndex + 2)
		{
			vertex2 = ReadList_cluster[ReadVertexIndex[(colIndex - 1 - transStartNum)/2 + transStartNum]];
			if (vertex2->start - vertex1->end >= 0)
			{
				gap = -MAX, length = MAX;

				subpaths = FragpathCluster[IndexLocator_fragcluster[(rowIndex - 1 - transStartNum)/2 + transStartNum][(colIndex - 1 - transStartNum)/2 + transStartNum]]->paths;
				if (subpaths != NULL)
				{
					connectedreadcluster = ConnectedReadCluster[IndexLocator_connectedReadscluter[(rowIndex - 1 - transStartNum)/2 + transStartNum][(colIndex - 1 - transStartNum)/2 + transStartNum]];
					if (connectedreadcluster->minLength_inbetween >= gapConstraint)
					{
						optimalweight[rowIndex + 1][colIndex] = connectedreadcluster->lengthlist->lowerbound;
					}
					else
					{
						gap = getOptimalProb(genevertex, adjASMinfoList, connectedreadcluster, sortlist, sortkey, dimension, subpaths);
						optimalweight[rowIndex + 1][colIndex] = gap;
						if (gap > maxgap)
						{
							maxgap = gap;
						}
					}
				}
			}
		}
	}

	delete [] LikelihoodDiff;
	delete [] sortlist;
	delete [] sortkey;

	// set the weight cost on edges
	sortlist = new void* [vertexNum + 2];
	sortkey = new double [vertexNum + 2];
	double* likelihood = new double[vertexNum];
	for (index = 0; index < vertexNum; index++)
	{
		likelihood[index] = -1;
	}

	index = 1;
	for (rowIndex = source + 1; rowIndex <= transStartNum; rowIndex++)
	{
		for (colIndex = transStartNum + 1; colIndex < sink - transEndNum; colIndex = colIndex + 2)
		{
			if (optimalweight[rowIndex][colIndex] > -MAX)
			{
				if (cost[rowIndex][colIndex] < MAX)
				{
					// do nothing
				}
				else
				{
					if (maxgap >= 0)
					{
						cost[rowIndex][colIndex] = -(optimalweight[rowIndex][colIndex]) + maxgap + 1; 
					}
					else
						cost[rowIndex][colIndex] = -optimalweight[rowIndex][colIndex];
				}
			}
			if (likelihood[colIndex] < 0)
			{
				likelihood[colIndex] = cost[rowIndex][colIndex];
			}
			else
			{
				if (likelihood[colIndex] > cost[rowIndex][colIndex])
				{
					likelihood[colIndex] = cost[rowIndex][colIndex];
				}
			}
		}
	}
	likelihood[sink] = cost[source][sink];
	index = 1;
	for (long tmpIndex = 0; tmpIndex < vertexNum; tmpIndex++)
	{
		if (likelihood[tmpIndex] > 0)
		{
			sortlist[index] = (void*)tmpIndex;
			sortkey[index] = likelihood[tmpIndex];
			index++;
		}
	}	
	mergeSort_general(sortlist, sortkey, index - 1);
	for (long tmpIndex = 1; tmpIndex <= index - 1; tmpIndex++)
	{
		connectedIndex[source].push_back((long)sortlist[tmpIndex]);
	}

	for (index = 0; index < vertexNum; index++)
	{
		likelihood[index] = -1;
	}
	for (rowIndex = sink - 1; rowIndex >= sink - transEndNum; rowIndex--)
	{
		for (colIndex = sink - transEndNum - 1; colIndex > transStartNum; colIndex = colIndex - 2)
		{
			if (optimalweight[colIndex][rowIndex] > -MAX)
			{
				if (cost[colIndex][rowIndex] < MAX)
				{
					// do nothing
				}
				else
				{
					if (maxgap >= 0)
					{
						cost[colIndex][rowIndex] = -(optimalweight[colIndex][rowIndex]) + maxgap + 1; 
					}
					else
						cost[colIndex][rowIndex] = -optimalweight[colIndex][rowIndex];
				}
			}
			if (likelihood[colIndex] < 0)
			{
				likelihood[colIndex] = cost[colIndex][rowIndex];
			}
			else
			{
				if (likelihood[colIndex] > cost[colIndex][rowIndex])
				{
					likelihood[colIndex] = cost[colIndex][rowIndex];
				}
			}
		}
	}

	for (rowIndex = transStartNum + 1; rowIndex < sink - transEndNum; rowIndex = rowIndex + 2)
	{
		colIndex = 1;
		if (likelihood[rowIndex + 1] > 0)
		{
			sortlist[colIndex] = (void*)sink;
			sortkey[colIndex] = likelihood[rowIndex + 1];
			colIndex++;
		}
		for (index = rowIndex + 2; index < sink - transEndNum; index = index + 2)
		{
			if (optimalweight[rowIndex + 1][index] > -MAX)
			{
				if (maxgap >= 0)
				{
					cost[rowIndex + 1][index] = -(optimalweight[rowIndex + 1][index]) + maxgap + 1; 
				}
				else
					cost[rowIndex + 1][index] = -optimalweight[rowIndex + 1][index];
				sortlist[colIndex] = (void*)index;
				sortkey[colIndex] = cost[rowIndex + 1][index];
				colIndex++;
			}
		}
		mergeSort_general(sortlist, sortkey, colIndex - 1);
		for (long tmpIndex = 1; tmpIndex <= colIndex - 1; tmpIndex++)
		{
			connectedIndex[rowIndex + 1].push_back((long)sortlist[tmpIndex]);
		}
	}

	for (rowIndex = 0; rowIndex < vertexNum; rowIndex++)
	{
		delete [] optimalweight[rowIndex];
	}
	delete [] optimalweight;
	delete [] likelihood;
	delete [] sortlist;
	delete [] sortkey;


	return;
}


// Initialize flow graph capability and cost 
void InitializeFlowGraph_Weibull(GTvertex *genevertex, long **cap, double **cost, vector< vector<int> > &connectedIndex, long source, long sink, long transStartNum, long transEndNum, long &maxSup, long totalReadNum, AdjASMsInfo *adjASMinfoList, long **IndexLocator_fragcluster, long **IndexLocator_connectedReadscluter, double scale)
{
	bool compatible, connected, pathsConstructed, found, transcriptionEnds;
	long rowIndex, colIndex, index, fragpathIndex, pathIndex, length, maxSup1, maxSup2, *lengthlist, dimension, genelength;
	graphpath *subpaths, *curpath, *paths;
	graphpathCluster *cluster;
	ReadVertex *vertex1, *vertex2;
	connectedRead *connectedreadcluster;
	pathlength_record *newRecord;

	genelength = 0;
	rangeJunction *curfrag = genevertex->junctionInRange->list;
	while(curfrag != NULL)
	{
		if (curfrag->junc->type != frag_junction)
		{
			genelength += curfrag->junc->end - curfrag->junc->start + 1;
		}
		curfrag = curfrag->next;
	}

	dimension = gapConstraint + 1;
	lengthlist = new long[dimension];
	for (rowIndex = 0; rowIndex < dimension; rowIndex++)
	{
		lengthlist[rowIndex] = 0;
	}
	subpaths = NULL;

	cap[source][sink] = MAX;
	//	cost[source][sink] = MAX/4;
	cost[source][sink] = 0;

	// Build edges from source to transcription starts and transcription ends to sink
	for (rowIndex = 1; rowIndex <= transStartNum; rowIndex++)
	{
		cap[source][rowIndex] = totalReadNum; // the capacity equals the total number of reads
		cost[source][rowIndex] = 0;
	}
	rowIndex = sink - 1; 
	index = 0;
	while(index < transEndNum)
	{
		cap[rowIndex][sink] = totalReadNum;
		cost[rowIndex][sink] = 0;
		rowIndex--;
		index++;
	}

	// Build edges from virtual transcription start/end to other internal vertices
	transcriptionEnds = true;
	maxSup1 = 0;
	for (rowIndex = source + 1; rowIndex <= transStartNum; rowIndex++)
	{
		vertex1 = ReadList_cluster[ReadVertexIndex[rowIndex - 1]];
		for (colIndex = transStartNum + 1; colIndex < sink - transEndNum; colIndex = colIndex + 2)
		{
			vertex2 = ReadList_cluster[ReadVertexIndex[(colIndex - 1 - transStartNum)/2 + transStartNum]];
			pathsConstructed = false;
			for (fragpathIndex = 0; fragpathIndex < FragpathCluster.size(); fragpathIndex++)
			{
				if (FragpathCluster[fragpathIndex]->fragstart->junc->start == vertex1->list->junc->start && FragpathCluster[fragpathIndex]->fragstart->junc->end == vertex1->list->junc->end && FragpathCluster[fragpathIndex]->fragend->junc->start == vertex2->list->junc->start && FragpathCluster[fragpathIndex]->fragend->junc->end == vertex2->list->junc->end)
				{
					subpaths = FragpathCluster[fragpathIndex]->paths;
					pathsConstructed = true;
					IndexLocator_fragcluster[rowIndex - 1][(colIndex - 1 - transStartNum)/2 + transStartNum] = fragpathIndex;
					break;
				}
			}
			if (pathsConstructed == false)
			{
				subpaths = generate_graphpath(vertex1->list, vertex2->list, transcriptionEnds);
				IndexLocator_fragcluster[rowIndex - 1][(colIndex - 1 - transStartNum)/2 + transStartNum] = FragpathCluster.size() - 1;
			}

			if (subpaths != NULL)
			{
				length = GetMinGraphpathLen(subpaths, vertex1, vertex2, pathIndex);
				connectedreadcluster = new connectedRead;
				connectedreadcluster->readlowIndex = rowIndex - 1;
				connectedreadcluster->readhighIndex = (colIndex - 1 - transStartNum)/2 + transStartNum;
				connectedreadcluster->minLength_inbetween = length;

				if (length <= gapConstraint)
				{
					maxSup1 += vertex2->includedReadNum;
					fillIngaps(connectedreadcluster, subpaths, vertex1, vertex2, true, lengthlist);
				}
				else
				{					
					connectedreadcluster->pathselectedIndex = getOptimalPathIndex(subpaths, vertex1, vertex2, length, true);
					newRecord = new pathlength_record;
					newRecord->lowerbound = length;
					newRecord->higherbound = length;
					newRecord->next = connectedreadcluster->lengthlist;
					connectedreadcluster->lengthlist = newRecord;
				}
				if (ConnectedReadCluster.size() >= ConnectedReadCluster.capacity())
					ConnectedReadCluster.reserve(default_dataset_num + ConnectedReadCluster.capacity());
				ConnectedReadCluster.push_back(connectedreadcluster);
				IndexLocator_connectedReadscluter[rowIndex - 1][(colIndex - 1 - transStartNum)/2 + transStartNum] = ConnectedReadCluster.size() - 1;

				cap[rowIndex][colIndex] = vertex2->includedReadNum;				
			}	
		}
	}

	maxSup2 = 0;
	for (rowIndex = sink - 1; rowIndex >= sink - transEndNum; rowIndex--)
	{
		vertex1 = ReadList_cluster[ReadVertexIndex[transStartNum - 1 + (vertexNum - 2 - transStartNum - transEndNum)/2 + rowIndex - (sink - transEndNum) + 1]];
		for (colIndex = sink - transEndNum - 1; colIndex > transStartNum; colIndex = colIndex - 2)
		{
			vertex2 = ReadList_cluster[ReadVertexIndex[(colIndex - 1 - transStartNum)/2 + transStartNum]];
			pathsConstructed = false;
			for (fragpathIndex = 0; fragpathIndex < FragpathCluster.size(); fragpathIndex++)
			{
				if (FragpathCluster[fragpathIndex]->fragstart->junc->start == vertex2->listTail->junc->start && FragpathCluster[fragpathIndex]->fragstart->junc->end == vertex2->listTail->junc->end && FragpathCluster[fragpathIndex]->fragend->junc->start == vertex1->list->junc->start && FragpathCluster[fragpathIndex]->fragend->junc->end == vertex1->list->junc->end)
				{
					subpaths = FragpathCluster[fragpathIndex]->paths;
					pathsConstructed = true;
					IndexLocator_fragcluster[(colIndex - 1 - transStartNum)/2 + transStartNum][transStartNum - 1 + (vertexNum - 2 - transStartNum - transEndNum)/2 + rowIndex - (sink - transEndNum) + 1] = fragpathIndex;
					break;
				}
			}
			if (pathsConstructed == false)
			{
				subpaths = generate_graphpath(vertex2->listTail, vertex1->list, transcriptionEnds);
				IndexLocator_fragcluster[(colIndex - 1 - transStartNum)/2 + transStartNum][transStartNum - 1 + (vertexNum - 2 - transStartNum - transEndNum)/2 + rowIndex - (sink - transEndNum) + 1] = FragpathCluster.size() - 1;
			}
			if (subpaths != NULL)
			{
				connectedreadcluster = new connectedRead;				
				connectedreadcluster->readlowIndex = (colIndex - 1 - transStartNum)/2 + transStartNum;
				connectedreadcluster->readhighIndex = transStartNum - 1 + (vertexNum - 2 - transStartNum - transEndNum)/2 + rowIndex - (sink - transEndNum) + 1;

				length = GetMinGraphpathLen(subpaths, vertex2, vertex1, pathIndex);
				connectedreadcluster->minLength_inbetween = length;
				if (length <= gapConstraint)
				{
					maxSup2 += vertex2->includedReadNum;
					fillIngaps(connectedreadcluster,  subpaths, vertex2, vertex1, true, lengthlist);
				}
				else
				{				
					newRecord = new pathlength_record;
					newRecord->lowerbound = length;
					newRecord->higherbound = length;
					newRecord->next = connectedreadcluster->lengthlist;
					connectedreadcluster->lengthlist = newRecord;
					connectedreadcluster->pathselectedIndex = getOptimalPathIndex(subpaths, vertex2, vertex1, length, true);
				}
				if (ConnectedReadCluster.size() >= ConnectedReadCluster.capacity())
					ConnectedReadCluster.reserve(default_dataset_num + ConnectedReadCluster.capacity());
				ConnectedReadCluster.push_back(connectedreadcluster);
				IndexLocator_connectedReadscluter[(colIndex - 1 - transStartNum)/2 + transStartNum][transStartNum - 1 + (vertexNum - 2 - transStartNum - transEndNum)/2 + rowIndex - (sink - transEndNum) + 1] = ConnectedReadCluster.size() - 1;

				cap[colIndex][rowIndex] = vertex2->includedReadNum;				
			}	
		}
	}

	// Build edges among the internal vertices
	rangeJunction *frag1, *frag2;
	transcriptionEnds = false;
	for (rowIndex = transStartNum + 1; rowIndex < sink - transEndNum; rowIndex = rowIndex + 2)
	{
		vertex1 = ReadList_cluster[ReadVertexIndex[(rowIndex - 1 - transStartNum)/2 + transStartNum]];
		cost[rowIndex][rowIndex + 1] = epsilon;
		frag1 = vertex1->list;
		while(frag1 != NULL)
		{
			if (vertex1->end >= frag1->junc->start && vertex1->end <= frag1->junc->end && frag1->junc->type != frag_junction)
			{
				break;
			}
			frag1 = frag1->next;
		}
		for (index = rowIndex + 2; index < sink - transEndNum; index = index + 2)
		{
			vertex2 = ReadList_cluster[ReadVertexIndex[(index - 1 - transStartNum)/2 + transStartNum]];	
			if (vertex2->start - vertex1->end >= 0)
			{
				// find the two rangeJunctions we want to build paths in between
				frag2 = vertex2->list;
				while(frag2 != NULL)
				{
					if (vertex2->start >= frag2->junc->start && vertex2->start <= frag2->junc->end && frag2->junc->type != frag_junction)
					{
						break;
					}
					frag2 = frag2->next;
				}
				pathsConstructed = false;
				for (fragpathIndex = 0; fragpathIndex < FragpathCluster.size(); fragpathIndex++)
				{
					if (FragpathCluster[fragpathIndex]->fragstart->junc->start == frag1->junc->start && FragpathCluster[fragpathIndex]->fragstart->junc->end == frag1->junc->end && FragpathCluster[fragpathIndex]->fragend->junc->start == frag2->junc->start && FragpathCluster[fragpathIndex]->fragend->junc->end == frag2->junc->end)
					{
						subpaths = FragpathCluster[fragpathIndex]->paths;
						pathsConstructed = true;
						IndexLocator_fragcluster[(rowIndex - 1 - transStartNum)/2 + transStartNum][(index - 1 - transStartNum)/2 + transStartNum] = fragpathIndex;
						break;
					}
				}
				if (pathsConstructed == false)
				{
					subpaths = generate_graphpath(frag1, frag2, transcriptionEnds);
					IndexLocator_fragcluster[(rowIndex - 1 - transStartNum)/2 + transStartNum][(index - 1 - transStartNum)/2 + transStartNum] = FragpathCluster.size() - 1;
				}
				if (subpaths != NULL)
				{
					connectedreadcluster = new connectedRead;				
					connectedreadcluster->readlowIndex = (rowIndex - 1 - transStartNum)/2 + transStartNum;
					connectedreadcluster->readhighIndex = (index - 1 - transStartNum)/2 + transStartNum;

					length = GetMinGraphpathLen(subpaths, vertex1, vertex2, pathIndex);
					connectedreadcluster->minLength_inbetween = length;
					if (length < gapConstraint)
					{
						cap[rowIndex + 1][index] = min(vertex1->includedReadNum, vertex2->includedReadNum);
						fillIngaps(connectedreadcluster, subpaths, vertex1, vertex2, false, lengthlist);
					}
					else
					{
						newRecord = new pathlength_record;
						newRecord->lowerbound = length;
						newRecord->higherbound = length;
						newRecord->next = connectedreadcluster->lengthlist;
						connectedreadcluster->lengthlist = newRecord;
					}
					if (ConnectedReadCluster.size() >= ConnectedReadCluster.capacity())
						ConnectedReadCluster.reserve(default_dataset_num + ConnectedReadCluster.capacity());
					ConnectedReadCluster.push_back(connectedreadcluster);
					IndexLocator_connectedReadscluter[(rowIndex - 1 - transStartNum)/2 + transStartNum][(index - 1 - transStartNum)/2 + transStartNum] = ConnectedReadCluster.size() - 1;

				}
			}
		}

	}
	maxSup = max(maxSup1, maxSup2);

	// speed up process
	double *LikelihoodDiff;
	LikelihoodDiff = new double[dimension];
	for (index = 0; index < dimension; index++)
	{
		LikelihoodDiff[index] = WeibullDist(index, scale, genelength);
	}
	delete [] lengthlist;
	// sort LikelihoodDiff from largest to smallest
	void** sortlist;
	double* sortkey;
	sortlist = new void* [dimension + 2];
	sortkey = new double [dimension + 2];
	for (long tmpIndex = 0; tmpIndex < dimension; tmpIndex++)
	{
		sortlist[tmpIndex + 1] = (void*)tmpIndex;
		sortkey[tmpIndex + 1] = LikelihoodDiff[tmpIndex];
	}	
	mergeSort_general(sortlist, sortkey, dimension);

	double temp, gap = -MAX, maxgap = -MAX, **optimalweight;
	optimalweight = new double*[vertexNum];
	for (rowIndex = 0; rowIndex < vertexNum; rowIndex++)
	{
		optimalweight[rowIndex] = new double[vertexNum];
		for (colIndex = 0; colIndex < vertexNum; colIndex++)
		{
			optimalweight[rowIndex][colIndex] = -MAX;
		}
	}
	for (rowIndex = source + 1; rowIndex <= transStartNum; rowIndex++)
	{
		vertex1 = ReadList_cluster[ReadVertexIndex[rowIndex - 1]];
		for (colIndex = transStartNum + 1; colIndex < sink - transEndNum; colIndex = colIndex + 2)
		{
			vertex2 = ReadList_cluster[ReadVertexIndex[(colIndex - 1 - transStartNum)/2 + transStartNum]];
			gap = -MAX, length = MAX;
			// choose the one length with maximal weight
			subpaths = FragpathCluster[IndexLocator_fragcluster[rowIndex - 1][(colIndex - 1 - transStartNum)/2 + transStartNum]]->paths;
			if (subpaths != NULL)
			{
				connectedreadcluster = ConnectedReadCluster[IndexLocator_connectedReadscluter[rowIndex - 1][(colIndex - 1 - transStartNum)/2 + transStartNum]];
				if (connectedreadcluster->lengthlist->lowerbound == 0 && connectedreadcluster->lengthlist->higherbound == 0)
				{
					cost[rowIndex][colIndex] = epsilon;
					optimalweight[rowIndex][colIndex] = epsilon;
					length = 0;
					connectedreadcluster->pathselectedIndex = getOptimalPathIndex(subpaths, vertex1, vertex2, length, true);
				}
				else
				{
					if (connectedreadcluster->minLength_inbetween >= gapConstraint)
					{
						cost[rowIndex][colIndex] = connectedreadcluster->lengthlist->lowerbound;
					}
					else
					{
						gap = getOptimalProb(genevertex, adjASMinfoList, connectedreadcluster, sortlist, sortkey, dimension, subpaths);
						optimalweight[rowIndex][colIndex] = gap;
						if (gap > maxgap)
						{
							maxgap = gap;
						}
					}

				}
			}
		}
	}
	for (rowIndex = sink - 1; rowIndex >= sink - transEndNum; rowIndex--)
	{
		vertex1 = ReadList_cluster[ReadVertexIndex[transStartNum - 1 + (vertexNum - 2 - transStartNum - transEndNum)/2 + rowIndex - (sink - transEndNum) + 1]];
		for (colIndex = sink - transEndNum - 1; colIndex > transStartNum; colIndex = colIndex - 2)
		{
			vertex2 = ReadList_cluster[ReadVertexIndex[(colIndex - 1 - transStartNum)/2 + transStartNum]];
			gap = -MAX, length = MAX;
			subpaths = FragpathCluster[IndexLocator_fragcluster[(colIndex - 1 - transStartNum)/2 + transStartNum][transStartNum - 1 + (vertexNum - 2 - transStartNum - transEndNum)/2 + rowIndex - (sink - transEndNum) + 1]]->paths;
			if (subpaths != NULL)
			{
				connectedreadcluster = ConnectedReadCluster[IndexLocator_connectedReadscluter[(colIndex - 1 - transStartNum)/2 + transStartNum][transStartNum - 1 + (vertexNum - 2 - transStartNum - transEndNum)/2 + rowIndex - (sink - transEndNum) + 1]];
				if (connectedreadcluster->lengthlist->lowerbound == 0 && connectedreadcluster->lengthlist->higherbound == 0)
				{
					cost[colIndex][rowIndex] = epsilon;
					optimalweight[colIndex][rowIndex] = epsilon;
					length = 0;
					connectedreadcluster->pathselectedIndex = getOptimalPathIndex(subpaths, vertex2, vertex1, length, true);
				}
				else
				{
					if (connectedreadcluster->minLength_inbetween >= gapConstraint)
					{
						cost[colIndex][rowIndex] = connectedreadcluster->lengthlist->lowerbound;
					}
					else
					{
						gap = getOptimalProb(genevertex, adjASMinfoList, connectedreadcluster, sortlist, sortkey, dimension, subpaths);
						optimalweight[colIndex][rowIndex] = gap;
						if (gap > maxgap)
						{
							maxgap = gap;
						}
					}

				}
			}

		}
	}

	for (rowIndex = transStartNum + 1; rowIndex < sink - transEndNum; rowIndex = rowIndex + 2)
	{
		vertex1 = ReadList_cluster[ReadVertexIndex[(rowIndex - 1 - transStartNum)/2 + transStartNum]];
		frag1 = vertex1->list;
		while(frag1 != NULL)
		{
			if (vertex1->end >= frag1->junc->start && vertex1->end <= frag1->junc->end && frag1->junc->type != frag_junction)
			{
				break;
			}
			frag1 = frag1->next;
		}
		for (colIndex = rowIndex + 2; colIndex < sink - transEndNum; colIndex = colIndex + 2)
		{
			vertex2 = ReadList_cluster[ReadVertexIndex[(colIndex - 1 - transStartNum)/2 + transStartNum]];
			if (vertex2->start - vertex1->end >= 0)
			{
				gap = -MAX, length = MAX;

				subpaths = FragpathCluster[IndexLocator_fragcluster[(rowIndex - 1 - transStartNum)/2 + transStartNum][(colIndex - 1 - transStartNum)/2 + transStartNum]]->paths;
				if (subpaths != NULL)
				{
					connectedreadcluster = ConnectedReadCluster[IndexLocator_connectedReadscluter[(rowIndex - 1 - transStartNum)/2 + transStartNum][(colIndex - 1 - transStartNum)/2 + transStartNum]];
					if (connectedreadcluster->minLength_inbetween >= gapConstraint)
					{
						optimalweight[rowIndex + 1][colIndex] = connectedreadcluster->lengthlist->lowerbound;
					}
					else
					{
						gap = getOptimalProb(genevertex, adjASMinfoList, connectedreadcluster, sortlist, sortkey, dimension, subpaths);
						optimalweight[rowIndex + 1][colIndex] = gap;
						if (gap > maxgap)
						{
							maxgap = gap;
						}
					}
				}
			}
		}
	}

	delete [] LikelihoodDiff;
	delete [] sortlist;
	delete [] sortkey;

	// set the weight cost on edges
	sortlist = new void* [vertexNum + 2];
	sortkey = new double [vertexNum + 2];
	double* likelihood = new double[vertexNum];
	for (index = 0; index < vertexNum; index++)
	{
		likelihood[index] = -1;
	}

	index = 1;
	for (rowIndex = source + 1; rowIndex <= transStartNum; rowIndex++)
	{
		for (colIndex = transStartNum + 1; colIndex < sink - transEndNum; colIndex = colIndex + 2)
		{
			if (optimalweight[rowIndex][colIndex] > -MAX)
			{
				if (cost[rowIndex][colIndex] < MAX)
				{
					// do nothing
				}
				else
				{
					if (maxgap >= 0)
					{
						cost[rowIndex][colIndex] = -(optimalweight[rowIndex][colIndex]) + maxgap + 1; 
					}
					else
						cost[rowIndex][colIndex] = -optimalweight[rowIndex][colIndex];
				}
			}
			if (likelihood[colIndex] < 0)
			{
				likelihood[colIndex] = cost[rowIndex][colIndex];
			}
			else
			{
				if (likelihood[colIndex] > cost[rowIndex][colIndex])
				{
					likelihood[colIndex] = cost[rowIndex][colIndex];
				}
			}
		}
	}
	likelihood[sink] = cost[source][sink];
	index = 1;
	for (long tmpIndex = 0; tmpIndex < vertexNum; tmpIndex++)
	{
		if (likelihood[tmpIndex] > 0)
		{
			sortlist[index] = (void*)tmpIndex;
			sortkey[index] = likelihood[tmpIndex];
			index++;
		}
	}	
	mergeSort_general(sortlist, sortkey, index - 1);
	for (long tmpIndex = 1; tmpIndex <= index - 1; tmpIndex++)
	{
		connectedIndex[source].push_back((long)sortlist[tmpIndex]);
	}

	for (index = 0; index < vertexNum; index++)
	{
		likelihood[index] = -1;
	}
	for (rowIndex = sink - 1; rowIndex >= sink - transEndNum; rowIndex--)
	{
		for (colIndex = sink - transEndNum - 1; colIndex > transStartNum; colIndex = colIndex - 2)
		{
			if (optimalweight[colIndex][rowIndex] > -MAX)
			{
				if (cost[colIndex][rowIndex] < MAX)
				{
					// do nothing
				}
				else
				{
					if (maxgap >= 0)
					{
						cost[colIndex][rowIndex] = -(optimalweight[colIndex][rowIndex]) + maxgap + 1; 
					}
					else
						cost[colIndex][rowIndex] = -optimalweight[colIndex][rowIndex];
				}
			}
			if (likelihood[colIndex] < 0)
			{
				likelihood[colIndex] = cost[colIndex][rowIndex];
			}
			else
			{
				if (likelihood[colIndex] > cost[colIndex][rowIndex])
				{
					likelihood[colIndex] = cost[colIndex][rowIndex];
				}
			}
		}
	}

	for (rowIndex = transStartNum + 1; rowIndex < sink - transEndNum; rowIndex = rowIndex + 2)
	{
		colIndex = 1;
		if (likelihood[rowIndex + 1] > 0)
		{
			sortlist[colIndex] = (void*)sink;
			sortkey[colIndex] = likelihood[rowIndex + 1];
			colIndex++;
		}
		for (index = rowIndex + 2; index < sink - transEndNum; index = index + 2)
		{
			if (optimalweight[rowIndex + 1][index] > -MAX)
			{
				if (maxgap >= 0)
				{
					cost[rowIndex + 1][index] = -(optimalweight[rowIndex + 1][index]) + maxgap + 1; 
				}
				else
					cost[rowIndex + 1][index] = -optimalweight[rowIndex + 1][index];
				sortlist[colIndex] = (void*)index;
				sortkey[colIndex] = cost[rowIndex + 1][index];
				colIndex++;
			}
		}
		mergeSort_general(sortlist, sortkey, colIndex - 1);
		for (long tmpIndex = 1; tmpIndex <= colIndex - 1; tmpIndex++)
		{
			connectedIndex[rowIndex + 1].push_back((long)sortlist[tmpIndex]);
		}
	}

	for (rowIndex = 0; rowIndex < vertexNum; rowIndex++)
	{
		delete [] optimalweight[rowIndex];
	}
	delete [] optimalweight;
	delete [] likelihood;
	delete [] sortlist;
	delete [] sortkey;


	return;
}


long getOptimalPathIndex(graphpath *paths, ReadVertex *vertex1, ReadVertex *vertex2, long length, bool transcriptionEnds)
{
	long pathIndex, length_min, length_max;
	graphpath *curpath;

	if (vertex1->end_max > vertex2->start_min)
	{
		// could only occur for internal nodes
		length_min = 0;
		curpath = paths;
		while(curpath != NULL)
		{
			if (curpath->edgelist == NULL)
			{
				length_max = vertex2->start - vertex1->end;
			}
			else
			{
				length_max = curpath->pathlength + vertex1->listTail->junc->end - vertex1->end + 1 + vertex2->start - vertex2->list->junc->start;
			}
			if (length >= length_min && length <= length_max)
			{
				pathIndex = curpath->pathindex;
				break;
			}

			curpath = curpath->next;
		}
	}
	else
	{
		curpath = paths;
		while(curpath != NULL)
		{
			if (curpath->edgelist == NULL)
			{
				length_min = vertex2->start_min - vertex1->end_max;
				length_max = vertex2->start - vertex1->end;
				if (transcriptionEnds == true)
				{
					if (abs(length_min - 1 ) <= 1)
					{
						pathIndex = curpath->pathindex;
						break;
					}
				}

			}
			else
			{
				length_min = curpath->pathlength + vertex1->listTail->junc->end - vertex1->end_max + 1 + vertex2->start_min - vertex2->list->junc->start;
				length_max = curpath->pathlength + vertex1->listTail->junc->end - vertex1->end + 1 + vertex2->start - vertex2->list->junc->start;
			}

			if (length >= length_min && length <= length_max)
			{
				pathIndex = curpath->pathindex;
				break;
			}
			curpath = curpath->next;
		}
	}
	return pathIndex;
}


// Update the gaps collected from one assembly result (debug use)
void updateGapList(long **flownet, long source, long sink, long transStartNum, long transEndNum)
{
	long rowIndex, index, iLoop, fragpathIndex, length, vertex1_Id, vertex2_Id;
	graphpath *graphpathread, *connectedPath;
	rangeJunction *fragstart, *fragend;

	// deal with source and sink first
	for (rowIndex = 1; rowIndex <= transStartNum; rowIndex++)
	{
		vertex1_Id = ReadVertexIndex[rowIndex - 1];
		fragstart = ReadList[vertex1_Id]->list;
		for (index = transStartNum + 1; index < sink - transEndNum; index = index + 2)
		{
			vertex2_Id = ReadVertexIndex[(index - 1 - transStartNum)/2 + transStartNum];
			fragend = ReadList[vertex2_Id]->list;
			if (flownet[rowIndex][index] == 1)
			{
				for (fragpathIndex = 0; fragpathIndex < FragpathCluster.size(); fragpathIndex++)
				{
					if (FragpathCluster[fragpathIndex]->fragstart->junc->start == fragstart->junc->start && FragpathCluster[fragpathIndex]->fragstart->junc->end == fragstart->junc->end && FragpathCluster[fragpathIndex]->fragend->junc->start == fragend->junc->start && FragpathCluster[fragpathIndex]->fragend->junc->end == fragend->junc->end)
					{
						graphpathread = FragpathCluster[fragpathIndex]->paths;
						if (graphpathread != NULL)
						{
							// add the gap length to gapList
							for (iLoop = 0; iLoop < ConnectedReadCluster.size(); iLoop++)
							{
								if (vertex1_Id == ConnectedReadCluster[iLoop]->readlowIndex && vertex2_Id == ConnectedReadCluster[iLoop]->readhighIndex)
								{
									break;
								}
							}
							connectedPath = graphpathread;
							while(connectedPath != NULL)
							{
								if (connectedPath->pathindex == ConnectedReadCluster[iLoop]->pathselectedIndex)
								{
									if (connectedPath->edgelist == NULL)
									{
										length = ReadList[vertex2_Id]->start - fragstart->junc->start;
										if (abs(length - 1) > 1)
										{
											if (gapList_partial.size() >= gapList_partial.capacity())
												gapList_partial.reserve(default_dataset_num + gapList_partial.capacity());
											gapList_partial.push_back(length);
										}	
									}
									else
									{
										if (gapList_partial.size() >= gapList_partial.capacity())
											gapList_partial.reserve(default_dataset_num + gapList_partial.capacity());
										length = connectedPath->pathlength + fragstart->junc->end - fragstart->junc->start + 1 + ReadList[vertex2_Id]->start - fragend->junc->start;
										gapList_partial.push_back(length);
									}
									break;
								}
								connectedPath = connectedPath->next;
							}

						}

					}
				}
			}
		}
	}

	for (rowIndex = sink - 1; rowIndex >= sink - transEndNum; rowIndex--)
	{
		vertex1_Id = ReadVertexIndex[transStartNum - 1 + (vertexNum - 2 - transStartNum - transEndNum)/2 + rowIndex - (sink - transEndNum) + 1];
		fragend = ReadList[vertex1_Id]->list;
		for (index = sink - transEndNum - 1; index > transStartNum; index = index - 2)
		{
			vertex2_Id = ReadVertexIndex[(index - 1 - transStartNum)/2 + transStartNum];
			fragstart = ReadList[vertex2_Id]->listTail;
			if (flownet[index][rowIndex] == 1)
			{
				for (fragpathIndex = 0; fragpathIndex < FragpathCluster.size(); fragpathIndex++)
				{
					if (FragpathCluster[fragpathIndex]->fragstart->junc->start == fragstart->junc->start && FragpathCluster[fragpathIndex]->fragstart->junc->end == fragstart->junc->end && FragpathCluster[fragpathIndex]->fragend->junc->start == fragend->junc->start && FragpathCluster[fragpathIndex]->fragend->junc->end == fragend->junc->end)
					{
						graphpathread = FragpathCluster[fragpathIndex]->paths;
						if (graphpathread != NULL)
						{
							// add the gap length to gapList
							for (iLoop = 0; iLoop < ConnectedReadCluster.size(); iLoop++)
							{
								if (vertex2_Id == ConnectedReadCluster[iLoop]->readlowIndex && vertex1_Id == ConnectedReadCluster[iLoop]->readhighIndex)
								{
									break;
								}
							}
							connectedPath = graphpathread;
							while(connectedPath != NULL)
							{
								if (connectedPath->pathindex == ConnectedReadCluster[iLoop]->pathselectedIndex)
								{
									if (connectedPath->edgelist == NULL)
									{
										length = fragend->junc->end - ReadList[vertex2_Id]->end;
										if (abs(length - 1) > 1)
										{
											if (gapList_partial.size() >= gapList_partial.capacity())
												gapList_partial.reserve(default_dataset_num + gapList_partial.capacity());
											gapList_partial.push_back(length);
										}	
									}
									else
									{
										if (gapList_partial.size() >= gapList_partial.capacity())
											gapList_partial.reserve(default_dataset_num + gapList_partial.capacity());
										length = connectedPath->pathlength + fragend->junc->end - fragend->junc->start + 1 + fragstart->junc->end - ReadList[vertex2_Id]->end;
										gapList_partial.push_back(length);
									}
									break;
								}
								connectedPath = connectedPath->next;
							}
						}

						break;
					}
				}
			}
		}
	}

	// deal with other edges on the flow
	for (rowIndex = transStartNum + 1; rowIndex < sink - transEndNum; rowIndex = rowIndex + 2)
	{
		vertex1_Id = ReadVertexIndex[(rowIndex - 1 - transStartNum)/2 + transStartNum];
		fragstart = ReadList[vertex1_Id]->listTail;
		for (index = rowIndex + 2; index < sink - transEndNum; index = index + 2)
		{
			vertex2_Id = ReadVertexIndex[(index - 1 - transStartNum)/2 + transStartNum];
			fragend = ReadList[vertex2_Id]->list;
			if (flownet[rowIndex + 1][index] == 1)
			{
				for (fragpathIndex = 0; fragpathIndex < FragpathCluster.size(); fragpathIndex++)
				{
					if (FragpathCluster[fragpathIndex]->fragstart->junc->start == fragstart->junc->start && FragpathCluster[fragpathIndex]->fragstart->junc->end == fragstart->junc->end && FragpathCluster[fragpathIndex]->fragend->junc->start == fragend->junc->start && FragpathCluster[fragpathIndex]->fragend->junc->end == fragend->junc->end)
					{
						graphpathread = FragpathCluster[fragpathIndex]->paths;
						if (graphpathread != NULL)
						{
							// add the gap length to gapList
							for (iLoop = 0; iLoop < ConnectedReadCluster.size(); iLoop++)
							{
								if (vertex1_Id == ConnectedReadCluster[iLoop]->readlowIndex && vertex2_Id == ConnectedReadCluster[iLoop]->readhighIndex)
								{
									break;
								}
							}
							connectedPath = graphpathread;
							while(connectedPath != NULL)
							{
								if (connectedPath->pathindex == ConnectedReadCluster[iLoop]->pathselectedIndex)
								{
									if (connectedPath->edgelist == NULL)
									{
										length = ReadList[vertex2_Id]->start - ReadList[vertex1_Id]->end - 1;	
									}
									else
									{
										length = connectedPath->pathlength + fragstart->junc->end - ReadList[vertex1_Id]->end + ReadList[vertex2_Id]->start - fragend->junc->start;
									}
									if (gapList_partial.size() >= gapList_partial.capacity())
										gapList_partial.reserve(default_dataset_num + gapList_partial.capacity());
									gapList_partial.push_back(length);

									break;
								}
								connectedPath = connectedPath->next;
							}
						}

						break;
					}
				}
			}

		}
	}

	return;
}

// Compare two read paths
bool CompareReadPaths(ReadIndicesPath *path1, ReadIndicesPath *path2)
{
	bool thesame = true;
	ReadIndex *curpath1, *curpath2;

	curpath1 = path1->readIndices;
	curpath2 = path2->readIndices;
	while(curpath1 != NULL && curpath2 != NULL)
	{
		if (curpath1->readIndex != curpath2->readIndex)
		{
			thesame = false;
			break;
		}
		curpath1 = curpath1->next;
		curpath2 = curpath2->next;
	}
	if (curpath1 != NULL || curpath2 != NULL)
	{
		thesame = false;
	}

	return thesame;
}

// From minimum cost maximum flow to path cover
ReadIndicesPath* FindMinPathCover(long **flownet, long source, long sink, long transStartNum, long transEndNum)
{
	bool compatible;
	long rowIndex, index, readpathIndex, flowIndex, vertex1_Id, vertex2_Id;
	ReadIndicesPath *curPath, *newPath, *curPathread, *resultreadpaths;
	ReadIndex *newindex;
	graphpath *graphpathread, *curgraphpath;
	JuncGraphEdge *graphedge;
	curPath = NULL;

	resultreadpaths = NULL;
	vector<ReadIndicesPath *> readpathQueue;

	// Transcription start
	for (rowIndex = 1; rowIndex <= transStartNum; rowIndex++)
	{
		for (index = transStartNum + 1; index < sink - transEndNum; index = index + 2)
		{
			vertex1_Id = ReadVertexIndex[(index - 1 - transStartNum)/2 + transStartNum];
			if (flownet[rowIndex][index] > 0)
			{
				flowIndex = 0;
				while(flowIndex < flownet[rowIndex][index])
				{
					newPath = new ReadIndicesPath;
					newindex = new ReadIndex;
					newindex->readIndex = rowIndex - 1;
					newPath->readIndices = newindex;
					newindex = new ReadIndex;
					newindex->readIndex = vertex1_Id;
					newPath->readIndices->next = newindex;
					newPath->lastreadIndex = newindex; 
					newPath->copyNum = 1;
					if (ReadList_cluster[vertex1_Id]->spliced == true)
					{
						newPath->transDirection = ReadList_cluster[vertex1_Id]->transDirection;
					}
					if (readpathQueue.size() >= readpathQueue.capacity())
						readpathQueue.reserve(default_dataset_num + readpathQueue.capacity());
					readpathQueue.push_back(newPath);

					flowIndex++;
				}
			}
		}
	}

	// Internal edges
	for (rowIndex = transStartNum + 1; rowIndex < sink - transEndNum; rowIndex = rowIndex + 2)
	{
		vertex1_Id = ReadVertexIndex[(rowIndex - 1 - transStartNum)/2 + transStartNum];
		for (index = rowIndex + 2; index < sink - transEndNum; index = index + 2)
		{
			if (flownet[rowIndex + 1][index] > 0)
			{
				flowIndex = 0;
				vertex2_Id = ReadVertexIndex[(index - 1 - transStartNum)/2 + transStartNum];

				readpathIndex = 0;
				while(readpathIndex < readpathQueue.size())
				{
					curPath = readpathQueue[readpathIndex];

					if (curPath->lastreadIndex->readIndex == vertex1_Id)
					{		
						newindex = new ReadIndex;
						newindex->readIndex = vertex2_Id;
						curPath->lastreadIndex->next = newindex;
						curPath->lastreadIndex = newindex;
						curPath->lastreadIndex->next = NULL;
						if (curPath->transDirection == undetermined)
						{
							if (ReadList_cluster[vertex2_Id]->spliced == true)
							{
								curPath->transDirection = ReadList_cluster[vertex2_Id]->transDirection;
							}
						}
						flowIndex++;

						if (flowIndex >= flownet[rowIndex + 1][index])
						{
							break;
						}
					}

					readpathIndex++;
				}
				if(flowIndex < flownet[rowIndex + 1][index])
				{
					newPath = new ReadIndicesPath;
					newindex = new ReadIndex;
					newindex->readIndex = vertex1_Id;
					newPath->readIndices = newindex;
					newindex = new ReadIndex;
					newindex->readIndex = vertex2_Id;
					newPath->readIndices->next = newindex;
					newPath->lastreadIndex = newindex; 
					newPath->copyNum = 1;
					if (ReadList_cluster[vertex1_Id]->spliced == true)
					{
						newPath->transDirection = ReadList_cluster[vertex1_Id]->transDirection;
					}
					else if (ReadList_cluster[vertex2_Id]->spliced == true)
					{
						newPath->transDirection = ReadList_cluster[vertex2_Id]->transDirection;
					}
					if (readpathQueue.size() >= readpathQueue.capacity())
						readpathQueue.reserve(default_dataset_num + readpathQueue.capacity());
					readpathQueue.push_back(newPath);
				}

			}
		}
	}

	// transcription end
	bool found;
	int i;
	for (readpathIndex = 0; readpathIndex < readpathQueue.size(); readpathIndex++)
	{
		found = false;
		for (rowIndex = sink - 1; rowIndex >= sink - transEndNum; rowIndex--)
		{
			vertex2_Id = ReadVertexIndex[transStartNum - 1 + (vertexNum - 2 - transStartNum - transEndNum)/2 + rowIndex - (sink - transEndNum) + 1];
			if (flownet[(readpathQueue[readpathIndex]->lastreadIndex->readIndex - transStartNum + 1) * 2 + transStartNum][rowIndex] > 0)
			{
				found = true;
				newindex = new ReadIndex;
				newindex->readIndex = vertex2_Id;
				readpathQueue[readpathIndex]->lastreadIndex->next = newindex;
				readpathQueue[readpathIndex]->lastreadIndex = newindex; 
				readpathQueue[readpathIndex]->lastreadIndex->next = NULL;
				break;
			}
		}
		if (found == false)
		{
			cout << "This path can not be constructed! Please check...\t" << readpathIndex << endl;
			cin >> i;
		}
	}

	// return the set of assembled fragments
	readpathIndex = 0;
	while(readpathIndex < readpathQueue.size())
	{
		curPathread = readpathQueue[readpathIndex];

		if (resultreadpaths == NULL)
		{
			curPathread->next = resultreadpaths;
			resultreadpaths = curPathread;
		}
		else
		{
			found = false;
			curPath = resultreadpaths;
			while(curPath != NULL)
			{
				found = CompareReadPaths(curPath, curPathread);
				if (found == true)
				{
					break;
				}
				curPath = curPath->next;
			}
			if (found == false)
			{
				curPathread->next = resultreadpaths;
				resultreadpaths = curPathread;
			}
			else
			{
				curPath->copyNum++;
			}
		}
		readpathIndex++;
	}

	readpathQueue.clear();

	return resultreadpaths;
}

long NumFragsInPaths(long **flownet, long transStartNum, long transEndNum, long source, long sink)
{
	long fragNum = 0, rowIndex, index;

	// Internal edges
	for (rowIndex = transStartNum + 1; rowIndex < sink - transEndNum; rowIndex = rowIndex + 2)
	{
		if (flownet[rowIndex][rowIndex + 1] == 1)
		{
			fragNum++;
		}
	}

	return fragNum;
}


// For debug use only
// output the fragment paths
void OutputFragPaths(ReadPath *fragpaths, string Destination)
{
	ofstream outputfile;
	ReadPath *curpath;
	ReadVertex *curvertex;
	graphpath *graphpathread, *connectedPath;

	string filename;
	filename = Destination + "_readpaths.txt";
	long totalReadNum = 0, iLoop;

	outputfile.open(filename.c_str());
	curpath = fragpaths;
	while(curpath != NULL)
	{
		curvertex = curpath->pathVertices;
		while(curvertex->next != NULL)
		{
			totalReadNum++;
			// output curvertex and gap
			outputfile << curvertex->readname << "(" << curvertex->start << "-" << curvertex->end << ")" << "\t";
			for (long index = 0; index < FragpathCluster.size(); index++)
			{
				if (FragpathCluster[index]->fragstart->junc->start == curvertex->listTail->junc->start && FragpathCluster[index]->fragstart->junc->end == curvertex->listTail->junc->end && FragpathCluster[index]->fragend->junc->start == curvertex->next->list->junc->start && FragpathCluster[index]->fragend->junc->end == curvertex->next->list->junc->end)
				{
					graphpathread = FragpathCluster[index]->paths;
					if (graphpathread == NULL)
					{
						// do nothing
					}
					else
					{
						for (iLoop = 0; iLoop < ConnectedReadCluster.size(); iLoop++)
						{
							if ((curvertex->Id-1) == ConnectedReadCluster[iLoop]->readlowIndex && (curvertex->next->Id-1) == ConnectedReadCluster[iLoop]->readhighIndex)
							{
								break;
							}
						}
						int pathIndex = 0;
						connectedPath = graphpathread;
						while(connectedPath != NULL)
						{
							if (pathIndex == ConnectedReadCluster[iLoop]->pathselectedIndex)
							{
								if (connectedPath->edgelist == NULL)
								{
									outputfile << curvertex->next->start - curvertex->end - 1 << "\t";
								}
								else
								{
									outputfile << connectedPath->pathlength + curvertex->listTail->junc->end - curvertex->end + curvertex->next->start - curvertex->next->list->junc->start << "\t";
								}
								break;
							}
							pathIndex++;
							connectedPath = connectedPath->next;
						}

					}

					break;
				}
			}
			curvertex = curvertex->next;
		}
		totalReadNum++;
		outputfile << curvertex->readname << "(" << curvertex->start << "-" << curvertex->end << ")" << endl;
		curpath = curpath->next;
	}
	outputfile.close();
	//	cout << "Total number of reads:" << totalReadNum << endl;

	filename = Destination + "_gaps.txt";
	outputfile.open(filename.c_str());
	for (iLoop = 0; iLoop < gapList_partial.size(); iLoop++)
	{
		outputfile << gapList_partial[iLoop] << endl;
	}
	outputfile.close();
	gapList_partial.clear();

	return;
}


// get the supply or demand for each vertex
void GetSupply(long *supply, long source, long sink, long transStartNum, long transEndNum, long degree)
{
	long index;
	ReadVertex *vertex;
	for (index = transStartNum + 1; index < sink - transEndNum; index = index + 2)
	{
		vertex = ReadList_cluster[ReadVertexIndex[(index - 1 - transStartNum)/2 + transStartNum]];
		supply[index] = -vertex->includedReadNum;
		supply[index + 1] = vertex->includedReadNum;
	}
	supply[source] = degree;
	supply[sink] = -degree;

	return;
}

// get the total flow and the cost from the flow network
long ReverseTransform(long *supply, long **flownet, double **cost, long source, long sink, long transStartNum, long transEndNum, double &fcost)
{
	long flow = 0, vertex_Id;
	for (long rowIndex = source + 1; rowIndex <= source + transStartNum; rowIndex++)
	{
		if (flownet[source][rowIndex] > 0 )
		{
			flow += flownet[source][rowIndex];
		}
	}
	for (long rowIndex = transStartNum + 1; rowIndex < sink - transEndNum; rowIndex = rowIndex + 2)
	{
		vertex_Id = ReadVertexIndex[(rowIndex - 1 - transStartNum)/2 + transStartNum];
		flownet[rowIndex][rowIndex + 1] = ReadList_cluster[vertex_Id]->includedReadNum;
	}
	fcost = 0;
	for (long rowIndex = 0; rowIndex < vertexNum; rowIndex++)
	{
		for (long colIndex = 0; colIndex < vertexNum; colIndex++)
		{
			if (flownet[rowIndex][colIndex] > 0)
			{
				fcost += flownet[rowIndex][colIndex] * cost[rowIndex][colIndex];
			}
		}
	}

	return flow;
}

Transcript* CountCopyNum(ReadIndicesPath *pathList, long transStartNum, long transEndNum, long source, long sink, long **flownet, long **IndexLocator_fragcluster, long **IndexLocator_connectedReadscluter);


// Calculate the first degree derivative 
double firstDerivative(long maxSup, long minSup, long point, double *fcost_vec, bool *changed_vec, long **cap, double **cost, long *supply,  long **flownet, vector< vector<int> > connectedIndex, long source, long sink, long transStartNum, long transEndNum, vector< vector<int> > Indexassist)
{
	double derivative, fcost;
	long neighbor_point, flow;
	bool flowexist;
	neighbor_point = point + 1;
	if (neighbor_point > maxSup)
	{
		neighbor_point = point - 1;
	}
	if (fcost_vec[neighbor_point - minSup] > 0)
	{
		// do nothing
	}
	else
	{
		GetSupply(supply, source, sink, transStartNum, transEndNum, neighbor_point);
		flowexist = ShortestPath_MinCostFlow(vertexNum, source, sink, cost, cap, supply, flownet, connectedIndex, Indexassist);
		if (flowexist == true)
		{
			flow = ReverseTransform(supply, flownet, cost, source, sink, transStartNum, transEndNum, fcost);
			fcost_vec[neighbor_point - minSup] = fcost;
			changed_vec[neighbor_point - minSup] = true;
		}
	}

	derivative = (fcost_vec[point - minSup] - fcost_vec[neighbor_point - minSup]) / (point - neighbor_point);
	return derivative;
}

double round(double d)
{
	return floor(d + 0.5);
}

// update stepsize using Backtracking line search
int CaculateStepsize(int stepsize_prev, bool &changable, long supplyNum, double beta, double *fcost_vec, bool *changed_vec, long maxSup, long minSup, long **cap, double **cost, long *supply,  long **flownet, vector< vector<int> > connectedIndex, long source, long sink, long transStartNum, long transEndNum, vector< vector<int> > Indexassist)
{
	bool flowexist;
	double stepsize;
	long supplyNum_change, flow;
	double fcost;
	double derivative = firstDerivative(maxSup, minSup, supplyNum, fcost_vec, changed_vec, cap, cost, supply,  flownet, connectedIndex, source, sink, transStartNum, transEndNum, Indexassist);
	supplyNum_change = supplyNum - round(derivative);

	if (supplyNum_change >= minSup && supplyNum_change <= maxSup)
	{
		if (fcost_vec[supplyNum_change - minSup] > 0)
		{
			// do nothing
		}
		else
		{
			GetSupply(supply, source, sink, transStartNum, transEndNum, supplyNum_change);
			flowexist = ShortestPath_MinCostFlow(vertexNum, source, sink, cost, cap, supply, flownet, connectedIndex, Indexassist);
			if (flowexist == true)
			{
				flow = ReverseTransform(supply, flownet, cost, source, sink, transStartNum, transEndNum, fcost);
				fcost_vec[supplyNum_change - minSup] = fcost;
				changed_vec[supplyNum_change - minSup] = true;
			}
		}
		stepsize = stepsize_prev;
		if (fcost_vec[supplyNum_change - minSup] > fcost_vec[supplyNum - minSup] - stepsize/2*pow(derivative,2) && stepsize_prev == 1)
		{
			changable = false;
			return 1;
		}
		else
		{
			changable = true;
			while(fcost_vec[supplyNum_change - minSup] > fcost_vec[supplyNum - minSup] - stepsize/2*pow(derivative,2))
			{
				stepsize = beta * stepsize;
				if (stepsize < 1)
				{
					return 1;
				}
			}
			return round(stepsize);
		}	

	}
	else
	{
		stepsize = 0;
		changable = false;
	}
	return stepsize;
}

// Linearly finding supply_target
int FindSupply(long supplyNum_prev, long supplyNum_post, long minSup, double *fcost_vec, bool *changed_vec, long **cap, double **cost, long *supply, long **flownet, vector< vector<int> >connectedIndex, long source, long sink, long transStartNum, long transEndNum, vector< vector<int> > Indexassist)
{
	long supply_target, supplyNum, flow;
	bool flowexist;
	double fcost, fcost_target;

	fcost_target = MAX_NUMBER;
	supply_target = MAX_NUMBER;
	for (int supplyNum = supplyNum_prev; supplyNum <= supplyNum_post; supplyNum++)
	{
		if (fcost_vec[supplyNum - minSup] > 0)
		{
			// do nothing
		}
		else
		{
			GetSupply(supply, source, sink, transStartNum, transEndNum, supplyNum);
			flowexist = ShortestPath_MinCostFlow(vertexNum, source, sink, cost, cap, supply, flownet, connectedIndex, Indexassist);
			if(flowexist == true)
			{
				flow = ReverseTransform(supply, flownet, cost, source, sink, transStartNum, transEndNum, fcost);
				fcost_vec[supplyNum - minSup] = fcost;
				changed_vec[supplyNum - minSup] = true;
			}
		}
		if (fcost_vec[supplyNum - minSup] < fcost_target)
		{
			supply_target = supplyNum;
			fcost_target = fcost_vec[supplyNum - minSup];
		}
	}


	return supply_target;
}


// EM algorithm for transcript assembly
Transcript* EM_assembly(string filepath, long **cap, double **cost, long *supply, vector< vector<int> > connectedIndex, long source, long sink, long minSup, long maxSup, long transStartNum, long transEndNum, vector< vector<int> > Indexassist, long **IndexLocator_fragcluster, long **IndexLocator_connectedReadscluter)
{
	string filename;
	long rowIndex, colIndex, index, iLoop, **flownet, flow, supply_target, fragincluded, supplyNum;
	int peak, stepsize, iterCnt;
	double fcost, fcost_target, stepsize_temp;
	bool flowexist, optimal, changable;
	ReadIndicesPath *tmppathList, *del_path;
	Transcript *resultTransList;
	tmppathList = NULL;
	resultTransList = NULL;

	flownet = new long *[vertexNum];
	for (index = 0; index < vertexNum; index++)
	{
		flownet[index] = new long[vertexNum];
	}

	long supplyNum_prev, supplyNum_post, loop, minSupplyNum;
	vector<double> fcost_vec;
	vector<bool> changed_vec;
	supply_target = supplyNum_prev = supplyNum_post = MAX_NUMBER;
	fcost_target = MAX_NUMBER;

// 	while(supply_target == MAX_NUMBER)
// 	{
		for (supplyNum = maxSup; supplyNum >= minSup; supplyNum--)
		{
			GetSupply(supply, source, sink, transStartNum, transEndNum, supplyNum);
			flowexist = ShortestPath_MinCostFlow(vertexNum, source, sink, cost, cap, supply, flownet, connectedIndex, Indexassist);
			if (flowexist == true)
			{
				tmppathList = FindMinPathCover(flownet, source, sink, transStartNum, transEndNum);
				resultTransList = CountCopyNum(tmppathList, transStartNum, transEndNum, source, sink, flownet, IndexLocator_fragcluster, IndexLocator_connectedReadscluter);
				supply_target = supplyNum;
				break;
			}

		}
// 		minSup = maxSup + 1;
// 		maxSup = minSup * 2;
//	}

	fcost_vec.clear();
	changed_vec.clear();

	for (rowIndex = 0; rowIndex < vertexNum; rowIndex++)
	{
		delete [] flownet[rowIndex];
	}
	delete [] flownet;

	return resultTransList;
}

// Merge overlap transcript fragments
Transcript* Reorganize(Transcript *prematureTrans)
{	
	Transcript *TransList, *curTrans, *prevTrans, *newTrans, *Transpointer;
	long transIndex, transIndex_temp;
	bool exist;
	int relation;
	bool *flag; // record whether one transcript fragment can be merged into another fragment;

	TransList = NULL;

	transIndex = 0;
	curTrans = prematureTrans;
	while(curTrans != NULL)
	{
		if (TransList == NULL)
		{
			newTrans = curTrans->clone();
			newTrans->nextList = TransList;
			TransList = newTrans;
		}
		else
		{
			exist = false;
			prevTrans = NULL;
			Transpointer  = TransList;
			while(Transpointer != NULL)
			{
				relation = IdenticalTrans(Transpointer->transcript, curTrans->transcript);
				if (relation == 5 || relation == 1)
				{
					exist = true;
					prevTrans = Transpointer;
					Transpointer = Transpointer->nextList;
					break;
				}
				else if (relation == 2)
				{
					exist = true;
					newTrans = curTrans->clone();
					if (prevTrans == NULL)
					{
						newTrans->nextList = Transpointer->nextList;
						TransList = newTrans;
					}
					else
					{
						newTrans->nextList = Transpointer->nextList;
						prevTrans->nextList = newTrans;
					}
					delete Transpointer;
					Transpointer = newTrans;
					break;
				}
				else
				{
					prevTrans = Transpointer;
					Transpointer = Transpointer->nextList;
				}

			}
			if (exist == false)
			{
				newTrans = curTrans->clone();
				newTrans->nextList = TransList;
				TransList = newTrans;
			}
		}

		transIndex++;
		curTrans = curTrans->nextList;
	}

	// Remove redundant candidate transcript
	// 	Transcript *del_trans, *resultTransList;
	// 	while(prematureTrans != NULL)
	// 	{
	// 		del_trans = prematureTrans;
	// 		prematureTrans = del_trans->nextList;
	// 		delete del_trans;
	// 	}
	// 
	// 	resultTransList = RemoveRedundantTrans(TransList);
	// 	while(TransList != NULL)
	// 	{
	// 		del_trans = TransList;
	// 		TransList = del_trans->nextList;
	// 		delete del_trans;
	// 	}
	// 	return resultTransList;

	return TransList;
}


Transcript* CountCopyNum(ReadIndicesPath *pathList, long transStartNum, long transEndNum, long source, long sink, long **flownet, long **IndexLocator_fragcluster, long **IndexLocator_connectedReadscluter)
{
	ReadIndicesPath *curPath;
	RangeJunctionList *newTrans_path;
	rangeJunction *curfrag, *newfrag, *fragtail;
	ReadIndex *prevVertex, *curVertex, *nextVertex;
	long index, iLoop, transIndex = 0;
	graphpath *connectedPath;
	JuncGraphEdge *graphEdge;
	connectedRead *readcluster;
	int relation;
	long translength, copyCoveredlength;
	Transcript *resultTransList = NULL, *newTrans, *curTrans, *del_trans, *TransList, *copylist = NULL, *curTrans_path;

	int pathIndex = 0;
	curPath = pathList;
	while(curPath != NULL)
	{
		pathIndex++;

		translength = copyCoveredlength = 0;
		// build the transcript based on this molecule
		newTrans_path = new RangeJunctionList;
		newTrans_path->transDirection = curPath->transDirection;
		newTrans_path->nextList = NULL;

		prevVertex = NULL;
		curVertex = curPath->readIndices;
		newTrans_path->rangeLow = ReadList_cluster[curVertex->readIndex]->list->junc->start;
		if (curVertex->next == NULL)
		{
			// contains only one read
			curfrag = ReadList_cluster[curVertex->readIndex]->list;
			while(curfrag != NULL)
			{
				newfrag = new rangeJunction;
				newfrag->junc = curfrag->junc->clone();
				newfrag->next = NULL;


				// insert the newfrag into the 
				if (newTrans_path->list == NULL)
				{
					newTrans_path->list = newfrag;
					fragtail = newfrag;
				}
				else
				{
					fragtail->next = newfrag;
					fragtail = newfrag;
				}

				if (curfrag->junc->type != frag_junction)
				{
					translength += curfrag->junc->end - curfrag->junc->start + 1;
					if (curfrag == ReadList_cluster[curVertex->readIndex]->list)
					{
						copyCoveredlength += curfrag->junc->end - ReadList_cluster[curVertex->readIndex]->start_min + 1;
					}
					else if (curfrag == ReadList_cluster[curVertex->readIndex]->listTail)
					{
						copyCoveredlength += ReadList_cluster[curVertex->readIndex]->end_max - curfrag->junc->start + 1;
					}
					else
					{
						copyCoveredlength += curfrag->junc->end - curfrag->junc->start + 1;
					}
				}
				curfrag = curfrag->next;
			}
			newTrans_path->rangeHigh = ReadList_cluster[curVertex->readIndex]->listTail->junc->end;
		}
		else
		{
			while(curVertex->next != NULL)
			{
				if (prevVertex == NULL || (ReadList_cluster[prevVertex->readIndex]->listTail->junc->start < ReadList_cluster[curVertex->readIndex]->list->junc->start && ReadList_cluster[prevVertex->readIndex]->listTail->junc->end <= ReadList_cluster[curVertex->readIndex]->list->junc->start))
				{
					curfrag = ReadList_cluster[curVertex->readIndex]->list;
				}
				else
				{
					// find the rangeJunction we should start with
					curfrag = ReadList_cluster[curVertex->readIndex]->list;
					while(curfrag != NULL)
					{
						if (curfrag->junc->start == fragtail->junc->start && curfrag->junc->end == fragtail->junc->end)
						{
							break;
						}
						curfrag = curfrag->next;
					}
					if (curfrag != NULL)
					{
						curfrag = curfrag->next;
					}
				}
				while(curfrag != NULL)
				{
					newfrag = new rangeJunction;
					newfrag->junc = curfrag->junc->clone();
					newfrag->next = NULL;

					// insert the newfrag into the 
					if (newTrans_path->list == NULL)
					{
						newTrans_path->list = newfrag;
						fragtail = newfrag;
					}
					else
					{
						fragtail->next = newfrag;
						fragtail = newfrag;
					}

					if (curfrag->junc->type != frag_junction)
					{
						translength += curfrag->junc->end - curfrag->junc->start + 1;
						if (curfrag == ReadList_cluster[curVertex->readIndex]->list)
						{
							copyCoveredlength += curfrag->junc->end - ReadList_cluster[curVertex->readIndex]->start_min + 1;
						}
						else if (curfrag == ReadList_cluster[curVertex->readIndex]->listTail)
						{
							copyCoveredlength += ReadList_cluster[curVertex->readIndex]->end_max - curfrag->junc->start + 1;
						}
						else
						{
							copyCoveredlength += curfrag->junc->end - curfrag->junc->start + 1;
						}
					}

					curfrag = curfrag->next;
				}

				nextVertex = curVertex->next;
				if ((ReadList_cluster[curVertex->readIndex]->listTail->junc->start < ReadList_cluster[nextVertex->readIndex]->list->junc->start) && (ReadList_cluster[curVertex->readIndex]->listTail->junc->end <= ReadList_cluster[nextVertex->readIndex]->list->junc->start))
				{
					// there exist a path connecting two adjacent read vertices

					// add the connecting edges
					connectedPath = FragpathCluster[IndexLocator_fragcluster[curVertex->readIndex][nextVertex->readIndex]]->paths;
					readcluster = ConnectedReadCluster[IndexLocator_connectedReadscluter[curVertex->readIndex][nextVertex->readIndex]];
					while(connectedPath != NULL)
					{
						if (connectedPath->pathindex ==readcluster->pathselectedIndex)
						{
							graphEdge = connectedPath->edgelist;
							while(graphEdge->next != NULL)
							{
								newfrag = new rangeJunction;
								newfrag->junc = graphEdge->linkedVertex->corresJunc->junc->clone();
								newfrag->next = NULL;

								// insert the newfrag into the 
								fragtail->next = newfrag;
								fragtail = newfrag;

								if (newfrag->junc->type != frag_junction)
								{
									translength += newfrag->junc->end - newfrag->junc->start + 1;
								}

								graphEdge = graphEdge->next;
							}
							break;
						}
						connectedPath = connectedPath->next;
					}
				}
				else
				{
					// no connecting edges needed
					// do nothing
				}

				prevVertex = curVertex;
				curVertex = curVertex->next;
			}

			// add the fragments coming from the last vertex
			if (prevVertex == NULL || (ReadList_cluster[prevVertex->readIndex]->listTail->junc->start < ReadList_cluster[curVertex->readIndex]->list->junc->start && ReadList_cluster[prevVertex->readIndex]->listTail->junc->end <= ReadList_cluster[curVertex->readIndex]->list->junc->start))
			{
				curfrag = ReadList_cluster[curVertex->readIndex]->list;
			}
			else
			{
				// find the rangeJunction we should start with
				curfrag = ReadList_cluster[curVertex->readIndex]->list;
				while(curfrag != NULL)
				{
					if (curfrag->junc->start == fragtail->junc->start && curfrag->junc->end == fragtail->junc->end)
					{
						break;
					}
					curfrag = curfrag->next;
				}
				if (curfrag != NULL)
				{
					curfrag = curfrag->next;
				}
			}
			while(curfrag != NULL)
			{
				newfrag = new rangeJunction;
				newfrag->junc = curfrag->junc->clone();
				newfrag->next = NULL;

				// insert the newfrag into the 
				fragtail->next = newfrag;
				fragtail = newfrag;

				if (curfrag->junc->type != frag_junction)
				{
					translength += curfrag->junc->end - curfrag->junc->start + 1;
					if (curfrag == ReadList_cluster[curVertex->readIndex]->list)
					{
						copyCoveredlength += curfrag->junc->end - ReadList_cluster[curVertex->readIndex]->start_min + 1;
					}
					else if (curfrag == ReadList_cluster[curVertex->readIndex]->listTail)
					{
						copyCoveredlength += ReadList_cluster[curVertex->readIndex]->end_max - curfrag->junc->start + 1;
					}
					else
					{
						copyCoveredlength += curfrag->junc->end - curfrag->junc->start + 1;
					}
				}

				curfrag = curfrag->next;
			}
			newTrans_path->rangeHigh = ReadList_cluster[curVertex->readIndex]->listTail->junc->end;
		}

		if (copyCoveredlength > (double)translength *threshold_coverage)
		{
			MergeAdjExons(newTrans_path);

			// back up into a list
			newTrans = new Transcript;
			newTrans->transcript = newTrans_path;
			newTrans->copyNum = curPath->copyNum;
			newTrans->nextList = copylist;
			copylist = newTrans;

			// insert into the transcript list
			if (resultTransList == NULL)
			{
				newTrans = new Transcript;
				newTrans->transcript = newTrans_path->clone();
				newTrans->nextList = resultTransList;
				resultTransList = newTrans;
				transIndex++;
			}
			else
			{
				bool exist = false;
				curTrans = resultTransList;
				while(curTrans != NULL)
				{
					relation = IdenticalTrans(curTrans->transcript, newTrans_path);
					if (relation == 5)
					{
						exist = true;
						break;
					}
					curTrans = curTrans->nextList;
				}
				if (exist == false)
				{
					newTrans = new Transcript;
					newTrans->transcript = newTrans_path->clone();
					newTrans->nextList = resultTransList;
					resultTransList = newTrans;
					transIndex++;
				}
			}
		}

		curPath = curPath->next;
	}

	// sort the candidate transcript fragment set
	long transNum = transIndex;
	void** sortlist_trans = new void* [transNum + 2];
	double* sortkey_trans = new double [transNum + 2];

	curTrans = resultTransList;
	for (long tmpLoop = 1; tmpLoop <= transNum; ++tmpLoop)
	{
		sortlist_trans[tmpLoop] = (void*)curTrans;
		sortkey_trans[tmpLoop] = curTrans->transcript->rangeLow;
		curTrans = curTrans->nextList;
	}
	mergeSort_general(sortlist_trans, sortkey_trans, transNum);

	resultTransList = NULL;
	Transcript *tmpPointer;
	for (long tmpLoop = transNum; tmpLoop >= 1; --tmpLoop)
	{
		tmpPointer = (Transcript *)sortlist_trans[tmpLoop];
		tmpPointer->nextList = resultTransList;
		resultTransList = tmpPointer;
	}

	delete [] sortlist_trans;
	delete [] sortkey_trans;

	// 	ofstream outputfile;
	// 	outputfile.open("yan.txt");
	// 	curTrans = resultTransList;
	// 	while(curTrans != NULL)
	// 	{
	// 		curfrag = curTrans->transcript->list;
	// 		while(curfrag != NULL)
	// 		{
	// 			outputfile << curfrag->junc->start << "-" << curfrag->junc->end << "(";
	// 			if (curfrag->junc->type == frag_exon || curfrag->junc->type == frag_retained_intron)
	// 			{
	// 				outputfile << "1";
	// 			}
	// 			else
	// 				outputfile << "0";
	// 			outputfile << ")";
	// 			curfrag = curfrag->next;
	// 		}
	// 		outputfile << endl;
	// 		curTrans = curTrans->nextList;
	// 	}
	// 	outputfile.close();

	// Merge overlap transcript fragments;
	TransList = Reorganize(resultTransList);
	while(resultTransList != NULL)
	{
		del_trans = resultTransList;
		resultTransList = del_trans->nextList;
		delete del_trans;
	}
	resultTransList = TransList;

	// Count the molecule number of each assembled transcript;
	transIndex = 1;
	curTrans = resultTransList;
	while(curTrans != NULL)
	{
		curTrans->transName = itostr(GeneStart) + "-" + itostr(GeneEnd) + "_trans" + itostr(transIndex);
		curTrans_path = copylist;
		while(curTrans_path != NULL)
		{
			relation = IdenticalTrans(curTrans->transcript, curTrans_path->transcript);
			if (relation == 1 || relation == 5)
			{
				curTrans->copyNum += curTrans_path->copyNum;
			}
			curTrans_path = curTrans_path->nextList;
		}
		transIndex++;
		curTrans = curTrans->nextList;
	}

	while(copylist != NULL)
	{
		curTrans_path = copylist;
		copylist = curTrans_path->nextList;
		delete curTrans_path;
	}

	return resultTransList;
}


// Output GTF format
void outputGTF(Transcript *TransList, string filename, string chromosome, trans_direction direction)
{
	ofstream outputfile;
	long transIndex, readpathIndex, exonNum, junNum, readNum_used;
	Transcript *curTrans;
	rangeJunction *curfrag;
	ReadPath *curPath;
	ReadVertex *curVertex;
	string transname, genename;
	double temp = 0;
	if (direction == antisense)
	{
		genename = itostr(GeneStart) + "-" + itostr(GeneEnd) + "C";
	}
	else
		genename = itostr(GeneStart) + "-" + itostr(GeneEnd) + "W";
	readNum_used = 0;

	outputfile.open(filename.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
	if (TransList != NULL)
	{
		// output GTF file

		curTrans = TransList;
		transIndex = 0;
		while(curTrans != NULL)
		{
			exonNum = 0; 
			junNum = 0;

			transname = genename + "_trans" + itostr(transIndex+1);
			curfrag = curTrans->transcript->list;
			while(curfrag != NULL)
			{
				if (curfrag->junc->type == frag_junction)
				{
					junNum++;
				}
				else
				{
					outputfile << chromosome << '\t' << "Assembly\t" << "exon\t" << curfrag->junc->start << '\t' << curfrag->junc->end << '\t';
					outputfile << temp << '\t';
					if (direction == sense || direction == undetermined)
					{
						outputfile << "+\t";
					}
					else
					{
						outputfile << "-\t";
					}
					outputfile << ".\t";
					outputfile << "gene_id " << "\"" << genename << "\"; ";
					outputfile << "transcript_id " << "\"" << transname << "\";" << " ";
					outputfile << "copyNum " << "\"" << curTrans->copyNum << "\"; " << endl;
					exonNum++;
				}

				curfrag = curfrag->next;
			}

			transIndex++;
			curTrans = curTrans->nextList;
		}
	}

	outputfile.close();

	return;
}

void clearRedundantReads()
{
	while (ReadList.empty() != true)
	{
		delete ReadList[ReadList.size()-1];
		ReadList.pop_back();
	}

	while(ReadList_cluster.empty() != true)
	{
		delete ReadList_cluster[ReadList_cluster.size() - 1];
		ReadList_cluster.pop_back();
	}

	if (ReadVertexIndex.empty() != true)
	{
		ReadVertexIndex.clear();
	}

	return;
}


void clearPathclusters()
{
	while(ReadpathCluster.empty() != true)
	{
		delete ReadpathCluster[ReadpathCluster.size()-1];
		ReadpathCluster.pop_back();
	}

	while(FragpathCluster.empty() != true)
	{
		delete FragpathCluster[FragpathCluster.size()-1];
		FragpathCluster.pop_back();
	}

	while(ConnectedReadCluster.empty() != true)
	{
		delete ConnectedReadCluster[ConnectedReadCluster.size()-1];
		ConnectedReadCluster.pop_back();
	}

	return;
}

// check each retained intron to see whether it is an alternative start/end
bool AlternativeEnd(rangeJunction *frag)
{
	bool isalternative = true, found;
	rangeJunction *curfrag;
	for (long index = 0; index < ReadList.size(); index++)
	{
		found = false;
		curfrag = ReadList[index]->list;
		while(curfrag != NULL)
		{
			if (curfrag->junc->start == frag->junc->start && curfrag->junc->end == frag->junc->end)
			{
				found = true;
				break;
			}
			curfrag = curfrag->next;
		}
		if (found == true)
		{
			if (curfrag != ReadList[index]->list || curfrag != ReadList[index]->listTail)
			{
				isalternative = false;
				break;
			}
		}
	}
	return isalternative;
}

/************************************************************************/
/* Date: 09/19/2013 */
/* Function documents ASM and deterministic fragments information */
/************************************************************************/
vector<int> GetCoveredFrags(GTvertex *genevertex, rangeJunction *list)
{
	vector<int> FragIndexCovered;
	long fragIndex = 1;
	rangeJunction *frag, *readfrag, *firstfrag;
	bool found = false;

	frag = genevertex->junctionInRange->list;
	readfrag = list;
	firstfrag = list;
	while(frag != NULL && readfrag != NULL)
	{
		if (firstfrag->junc->start == frag->junc->start && firstfrag->junc->end == frag->junc->end)
		{
			found = true;
		}
		if (found == true)
		{
			if (readfrag->junc->start == frag->junc->start && readfrag->junc->end == frag->junc->end)
			{
				FragIndexCovered.push_back(fragIndex);
				readfrag = readfrag->next;
			}
		}
		frag = frag->next;
		fragIndex++;
	}

	return FragIndexCovered;
}

bool IdenticialReadInfo(ReadInfo *info1, ReadInfo *info2)
{
	bool identical = true;
	if (info1->FragIndexCovered.size() != info2->FragIndexCovered.size())
	{
		identical = false;
	}
	else
	{
		for (int fragIndex = 0; fragIndex < info1->FragIndexCovered.size(); fragIndex++)
		{
			if (info1->FragIndexCovered[fragIndex] != info2->FragIndexCovered[fragIndex])
			{
				identical = false;
				break;
			}
		}
	}
	return identical;
}

// sort the fragment list
void SortFragments(RangeJunctionList *lists)
{
	long fragNum;
	RangeJunctionList *curlist;
	rangeJunction *curfrag;

	curlist = lists;
	while(curlist != NULL)
	{
		fragNum = 0;
		curfrag = curlist->list;
		while(curfrag != NULL)
		{
			fragNum++;
			curfrag = curfrag->next;
		}

		void** sortlist_frag = new void* [fragNum + 2];
		double* sortkey_frag = new double [fragNum + 2];

		curfrag = curlist->list;
		for (long tmpLoop = 1; tmpLoop <= fragNum; ++tmpLoop)
		{
			sortlist_frag[tmpLoop] = (void*)curfrag;
			sortkey_frag[tmpLoop] = curfrag->junc->start;
			curfrag = curfrag->next;
		}
		mergeSort_general(sortlist_frag, sortkey_frag, fragNum);

		curlist->list = NULL;
		rangeJunction *tmpPointer;
		for (long tmpLoop = fragNum; tmpLoop >= 1; --tmpLoop)
		{
			tmpPointer = (rangeJunction *)sortlist_frag[tmpLoop];
			tmpPointer->next = curlist->list;
			curlist->list = tmpPointer;
		}

		delete [] sortlist_frag;
		delete [] sortkey_frag;

		curlist = curlist->nextList;
	}
}


AdjASMsInfo* DetermineFragInfo(GTvertex *genevertex, long MaxReadLength, long MeanReadLength)
{
	long ASMnumleft, readIndex, prevreadIndex = 0;
	GTvertex *VertexASM, *VertexASMnext, *vertex;
	GTedge *curEdge, *EdgeASM;
	AdjASMsInfo *newinfo, *infoList = NULL;
	ReadInfo *newreadinfo, *curreadinfo;
	TransInfo *newTransinfo;
	bool hasDeterministicRead = false, exist, separable;
	TransNum = 1;
	int prev_ASM_transNum = 0;

	if (genevertex->IndepASMNum > 1)
	{
		EdgeASM = genevertex->child;
		VertexASM = EdgeASM->linkedVertex;
		ASMnumleft = genevertex->IndepASMNum;

		int index = 1;
		rangeJunction *fragtail, *newfrag;

		while(ASMnumleft > 1)
		{
			newinfo = new AdjASMsInfo;
			if (VertexASM->childType > 1)
			{
				// reach an ASM
				hasDeterministicRead = false;
				newinfo->firstASM_start = VertexASM->rangeLow;
				newinfo->Hubstart = VertexASM->rangeHigh;
				VertexASMnext = VertexASM->nextSibling;

				// get all possible transcripts within these two adjacent ASMs
				rangeJunction *newfrag, *curfrag, *tailfrag;
				RangeJunctionList *subList = VertexASM->junctionInRange->clone(), *subpaths, *curpath;
				if (VertexASM->prevSibling != NULL)
				{
					newfrag = new rangeJunction;
					newfrag->junc = VertexASM->prevSibling->junctionInRange->list->junc->clone();
					newfrag->next = subList->list;
					subList->list = newfrag;
					subList->rangeLow = newfrag->junc->start;
				}
				tailfrag = subList->list;
				while(tailfrag->next != NULL)
				{
					tailfrag = tailfrag->next;
				}

				// deal with the connecting hubs
				newfrag = new rangeJunction;
				newfrag->junc = VertexASMnext->junctionInRange->list->junc->clone();
				tailfrag->next = newfrag;
				newfrag->next = NULL;
				tailfrag = newfrag;
				subList->rangeHigh = newfrag->junc->end;
				VertexASMnext =VertexASMnext->nextSibling;
				while(VertexASMnext->childType <= 1)
				{
					newfrag = new rangeJunction;
					newfrag->junc = VertexASMnext->junctionInRange->list->junc->clone();
					tailfrag->next = newfrag;
					newfrag->next = NULL;
					tailfrag = newfrag;
					subList->rangeHigh = newfrag->junc->end;
					VertexASMnext = VertexASMnext->nextSibling;
				}

				// the connecting ASM
				newinfo->secondASM_end = VertexASMnext->rangeHigh;
				newinfo->Hubend = VertexASMnext->rangeLow;
				curfrag = VertexASMnext->junctionInRange->list;
				while(curfrag != NULL)
				{
					newfrag = new rangeJunction;
					newfrag->junc = curfrag->junc->clone();
					tailfrag->next = newfrag;
					tailfrag = newfrag;
					subList->rangeHigh = newfrag->junc->end;
					curfrag = curfrag->next;
				}
				if (VertexASMnext->nextSibling != NULL)
				{
					newfrag = new rangeJunction;
					newfrag->junc = VertexASMnext->nextSibling->junctionInRange->list->junc->clone();
					tailfrag->next = newfrag;
					tailfrag = newfrag;
					subList->rangeHigh = newfrag->junc->end;
				}

				vertex = VertexASM->nextSibling;
				while(vertex != VertexASMnext)
				{
					if(vertex->junctionInRange->list->junc->type != frag_junction)
					{
						newinfo->Hublength += vertex->junctionInRange->list->junc->end - vertex->junctionInRange->list->junc->start + 1;
					}
					vertex = vertex->nextSibling;
				}
				
				subpaths = separateDepPath(subList, separable);
				SortFragments(subpaths);

				curpath = subpaths;
				while(curpath != NULL)
				{
					newTransinfo = new TransInfo;
					newTransinfo->FragIndexCovered = GetCoveredFrags(genevertex, curpath->list);
					curfrag = curpath->list;
					while(curfrag != NULL)
					{
						if (curfrag->junc->type != frag_junction)
						{
							newTransinfo->transLength += curfrag->junc->end - curfrag->junc->start + 1;
						}
						curfrag = curfrag->next;
					}
					newTransinfo->nextTransinfo = newinfo->ASM_Transinfo;
					newinfo->ASM_Transinfo = newTransinfo;
					curpath = curpath->nextList;
				}
				// delete temp result
				int ASM_transNum = 0;
				while(subpaths != NULL)
				{
					ASM_transNum++;
					curpath = subpaths;
					subpaths = curpath->nextList;
					delete curpath;
				}
				if (subList != NULL)
				{
					delete subList;
				}
				
				TransNum = (ASM_transNum/2) * (ASM_transNum - ASM_transNum/2 - prev_ASM_transNum);
				prev_ASM_transNum = ASM_transNum - ASM_transNum/2;

				if (MaxReadLength > newinfo->Hublength)
				{
					// Check for possible deterministic fragments
					for (long readIndex = prevreadIndex; readIndex < ReadList.size(); readIndex++)
					{
						if (ReadList[readIndex]->start >= newinfo->Hubstart)
						{
							prevreadIndex = readIndex;
							break;
						}
						if (ReadList[readIndex]->start < newinfo->Hubstart && ReadList[readIndex]->end > newinfo->Hubend)
						{
							hasDeterministicRead = true;
							newreadinfo = new ReadInfo;
							newreadinfo->FragIndexCovered = GetCoveredFrags(genevertex, ReadList[readIndex]->list);
							exist = false;
							curreadinfo = newinfo->deterministicReadinfo;
							while(curreadinfo != NULL)
							{
								if (IdenticialReadInfo(curreadinfo, newreadinfo) == true)
								{
									exist = true;
									delete newreadinfo;
									break;
								}
								curreadinfo = curreadinfo->nextreadInfo;
							}
							if (exist == false)
							{
								newreadinfo->nextreadInfo = newinfo->deterministicReadinfo;
								newinfo->deterministicReadinfo = newreadinfo;
							}
						}
					}
					if (hasDeterministicRead == false)
					{
						if (MeanReadLength > newinfo->Hublength + 1)
						{
							newinfo->nextinfo = infoList;
							infoList = newinfo;
						}
						else
							delete newinfo;
					}
					else
					{
						newinfo->nextinfo = infoList;
						infoList = newinfo;
					}
					
				}
				else
				{
					delete newinfo;
				}

				VertexASM = VertexASMnext;
				ASMnumleft--;
			}
			else
			{
				VertexASM = VertexASM->nextSibling;
			}
		}

// 		if (hasDeterministicRead == false)
// 		{
// 			AdjASMsInfo *delinfo;
// 			while(infoList != NULL)
// 			{
// 				delinfo = infoList;
// 				infoList = delinfo->nextinfo;
// 				delete delinfo;
// 			}
// 		}
	}

	return infoList;
}


// Filter lowly expressed transcripts
Transcript* FilterLowExpressionTrans(Transcript *TransList)
{
	int totalCopy = 0;
	Transcript *curTrans, *resultTransList = NULL, *newTrans;
	
	curTrans = TransList;
	while (curTrans != NULL)
	{
		totalCopy += curTrans->copyNum;
		curTrans = curTrans->nextList;
	}

	curTrans = TransList;
	while(curTrans != NULL)
	{
		if (curTrans->copyNum < ceil(threshold * totalCopy))
		{
			// this transcript is a minor transcript and can be removed
		}
		else
		{
			newTrans = curTrans->clone();
			newTrans->nextList = resultTransList;
			resultTransList = newTrans;
		}
		curTrans = curTrans->nextList;
	}

	return resultTransList;
}


// Shrink the read fragment cluster to proper size
void ShrinkClusterSize(GTvertex *genevertex, int tolerance)
{
	int index = 1;
	ClusterReadFragments(genevertex, tolerance);
	while (ReadList_cluster.size() > MAXCLUSTERSIZE && index < 5)
	{
		while(ReadList_cluster.empty() != true)
		{
			delete ReadList_cluster[ReadList_cluster.size() - 1];
			ReadList_cluster.pop_back();
		}
		ReadList_cluster.clear();
		tolerance += index * 20;
		ClusterReadFragments(genevertex, tolerance);
		index++;
	}

	return;
}

long Reconstruct_vertex(GTvertex *genevertex, string chromosome, string Destination, string SAMpath, bool &isAssembled, int tolerance, string outputfilename, double scale)
{
	string filename, filepath;
	rangeJunction *transStartfragList = NULL, *transEndfraglist = NULL, *frag1, *frag2, *newfrag;
	Transcript *TransList = NULL, *del_trans = NULL;
	long **cap, *supply, source, sink, minSup, maxSup, totalReadNum, transStartNum = 0, transEndNum = 0, clusterSize = 0, junctionNum, exoniclength, totalCov, **IndexLocator_fragcluster, **IndexLocator_connectedReadscluter;
	double **cost;
	bool altStart, altEnd;
	ReadVertex *newVertex;
	ReadIndex *curindex;
	trans_direction direction = Getdirection(genevertex);
	AdjASMsInfo *adjASMinfoList;
	isAssembled = false;


	GeneStart = genevertex->rangeLow;
	GeneEnd = genevertex->rangeHigh;

	// process SAM file
	RangeJunctionList *backupList;
	backupList = genevertex->junctionInRange->clone();
	spliceGraph = new JuncGraph;
	construct_splice_graph(backupList);
	backupList->list = NULL;

	junctionNum = 0;
	exoniclength = 0;
	frag1 = genevertex->junctionInRange->list;
	while(frag1 != NULL)
	{
		if (frag1->junc->type == frag_junction)
		{
			junctionNum++;
		}
		else
		{
			exoniclength += frag1->junc->end - frag1->junc->start + 1;
		}
		frag1 = frag1->next;
	}

	if (direction == antisense)
	{
		filename = SAMpath + itostr(GeneStart) + "-" + itostr(GeneEnd) + "C.txt";
	}
	else
		filename = SAMpath + itostr(GeneStart) + "-" + itostr(GeneEnd) + "W.txt";
	gapConstraint = 0;
	processReadFile(genevertex, filename);
	ReadNum = ReadList.size();
	totalCov = gapConstraint;
	adjASMinfoList = DetermineFragInfo(genevertex, MaxReadLength, MeanReadLength);

	if (ReadNum > 0)
	{
		if (junctionNum < MAXJUNTOLERENT && genevertex->IndepASMNum > 1)
		{
			// get the lists of transcription start/end fragments
			transStartfragList = transEndfraglist = NULL;
			frag1 = genevertex->junctionInRange->list;
			while(frag1 != NULL)
			{
				if (frag1->junc->type == frag_exon || frag1->junc->type == frag_retained_intron)
				{
					altStart = true;
					frag2 = genevertex->junctionInRange->list;
					while(frag2 != frag1)
					{
						if (abs(frag1->junc->start - frag2->junc->end) <= 1)
						{
							altStart = false;
							if (frag1->junc->type == frag_retained_intron)
							{
								altStart = AlternativeEnd(frag1);
							}
						}
						frag2 = frag2->next;
					}
					if (altStart == true)
					{
						newfrag = new rangeJunction;
						newfrag->junc = frag1->junc->clone();
						newfrag->next = transStartfragList;
						transStartfragList = newfrag;
						transStartNum++;
					}

					altEnd = true;
					frag2 = frag1->next;
					while(frag2 != NULL)
					{
						if (abs(frag2->junc->start - frag1->junc->end) <= 1)
						{
							altEnd = false;
							if(frag1->junc->type == frag_retained_intron)
							{
								altEnd = AlternativeEnd(frag1);
							}
						}
						frag2 = frag2->next;
					}
					if (altEnd == true)
					{
						newfrag = new rangeJunction;
						newfrag->junc = frag1->junc->clone();
						newfrag->next = transEndfraglist;
						transEndfraglist = newfrag;
						transEndNum++;
					}
				}
				frag1 = frag1->next;
			}

			gapConstraint = ceil((double)gapConstraint/ReadList.size());
			gapConstraint = gapConstraint * 3;
			if (gapConstraint > maxgapConstraint)
			{
				gapConstraint = maxgapConstraint;
			}
			totalReadNum = ReadList.size();
			// cluster the read fragments
//			tolerance = 0;
			ShrinkClusterSize(genevertex, tolerance);
			clusterSize = ReadList_cluster.size();
			while (ReadList.empty() != true)
			{
				delete ReadList[ReadList.size()-1];
				ReadList.pop_back();
			}

			if (clusterSize < MAX_READNUM_ALLOWED && clusterSize > 1)
			{
				isAssembled = true;
				// add transcription start/end as read into the ReadList
				frag1 = transStartfragList;
				while(frag1 != NULL)
				{
					newVertex = new ReadVertex;
					newVertex->start = frag1->junc->start;
					newVertex->end = frag1->junc->start;
					newVertex->start_min = newVertex->start;
					newVertex->end_max = newVertex->end;
					newfrag = new rangeJunction;
					newfrag->junc = frag1->junc->clone();
					newVertex->list = newVertex->listTail = newfrag;

					ReadList_cluster.insert(ReadList_cluster.begin(), newVertex);

					frag1 = frag1->next;
				}


				frag2 = transEndfraglist;
				while(frag2 != NULL)
				{
					newVertex = new ReadVertex;
					newVertex->start = frag2->junc->end;
					newVertex->end = frag2->junc->end;
					newVertex->start_min = newVertex->start;
					newVertex->end_max = newVertex->end;
					newfrag = new rangeJunction;
					newfrag->junc = frag2->junc->clone();
					newVertex->list = newVertex->listTail = newfrag;

					if (ReadList_cluster.size() >= ReadList_cluster.capacity())
						ReadList_cluster.reserve(default_dataset_num + ReadList_cluster.capacity());
					ReadList_cluster.push_back(newVertex);

					frag2 = frag2->next;
				}

				for (long index = 0; index < ReadList_cluster.size(); index++)
				{
					ReadList_cluster[index]->Id = index + 1;
				}

				if (ReadList_cluster.size() > transStartNum + transEndNum)
				{
					for (long index = 0; index < ReadList_cluster.size(); index++)
					{
						if (ReadVertexIndex.size() >= ReadVertexIndex.capacity())
							ReadVertexIndex.reserve(default_dataset_num + ReadVertexIndex.capacity());
						ReadVertexIndex.push_back(index);
					}

					// initialize fragment graph
					vertexNum = (ReadVertexIndex.size() - transStartNum - transEndNum) * 2 + 2 + transStartNum + transEndNum; // split each internal vertices and plus source/sink
					cap = new long*[vertexNum];
					cost = new double *[vertexNum];
					supply = new long[vertexNum];
					for (long rowIndex = 0; rowIndex < vertexNum; rowIndex++)
					{
						cap[rowIndex] = new long[vertexNum];
						cost[rowIndex] = new double[vertexNum];
						supply[rowIndex] = 0;
						for (long index = 0; index < vertexNum; index++)
						{
							cap[rowIndex][index] = 0;
							cost[rowIndex][index] = MAX;
						}
					}
					source = 0;
					sink = vertexNum - 1;


					IndexLocator_fragcluster = new long*[ReadList_cluster.size()];
					IndexLocator_connectedReadscluter = new long*[ReadList_cluster.size()];
					for (long rowIndex = 0; rowIndex < ReadList_cluster.size(); rowIndex++)
					{
						IndexLocator_fragcluster[rowIndex] = new long[ReadList_cluster.size()];
						IndexLocator_connectedReadscluter[rowIndex] = new long[ReadList_cluster.size()];
						for (long colIndex = 0; colIndex < ReadList_cluster.size(); colIndex++)
						{
							IndexLocator_fragcluster[rowIndex][colIndex] = -1;
							IndexLocator_connectedReadscluter[rowIndex][colIndex] = -1;
						}
					}

					vector< vector<int> > connectedIndex(vertexNum, vector<int>(0)), Indexassist(vertexNum, vector<int>(0));
					InitializeFlowGraph_Weibull(genevertex, cap, cost, connectedIndex, source, sink, transStartNum, transEndNum, maxSup, totalReadNum, adjASMinfoList, IndexLocator_fragcluster, IndexLocator_connectedReadscluter, scale);
					minSup = sourceDegree;
					if (maxSup > minSup * 2 || maxSup < minSup)
					{
						maxSup = minSup * 2;
					}
					for (long rowIndex = 0; rowIndex < vertexNum; rowIndex++)
					{
						for (long colIndex = 0; colIndex < connectedIndex[rowIndex].size(); colIndex++)
						{
							Indexassist[connectedIndex[rowIndex][colIndex]].push_back(rowIndex);
						}
					}

					filepath = Destination;
					TransList = EM_assembly(filepath, cap, cost, supply, connectedIndex, source, sink, minSup, maxSup, transStartNum, transEndNum, Indexassist, IndexLocator_fragcluster, IndexLocator_connectedReadscluter);
					Transcript *tempTransList = FilterLowExpressionTrans(TransList);
					while(TransList != NULL)
					{
						del_trans = TransList;
						TransList = del_trans->nextList;
						delete del_trans;
					}
					TransList = tempTransList;

					// clear all
					for (long rowIndex = 0; rowIndex < vertexNum; rowIndex++)
					{
						delete [] cap[rowIndex];
						delete [] cost[rowIndex];
						connectedIndex[rowIndex].clear();
						Indexassist[rowIndex].clear();
					}
					for (long rowIndex = 0; rowIndex < ReadList_cluster.size(); rowIndex++)
					{
						delete [] IndexLocator_fragcluster[rowIndex];
						delete [] IndexLocator_connectedReadscluter[rowIndex];
					}
					delete [] cap;
					delete [] cost;
					delete [] supply;
					delete [] IndexLocator_fragcluster;
					delete [] IndexLocator_connectedReadscluter;
					connectedIndex.clear();
					Indexassist.clear();
				}
			}

// 			if (direction == antisense)
// 			{
// 				filename = Destination + itostr(GeneStart) + "-" + itostr(GeneEnd) + "C.gtf";
// 			}
// 			else
// 				filename = Destination + itostr(GeneStart) + "-" + itostr(GeneEnd) + "W.gtf";

			if (TransList == NULL)
			{
				ofstream outputfile;
				string genename, transname;
				double temp = 0;
				int copyNum = round(totalCov/exoniclength);
				if (copyNum == 0)
				{
					copyNum = 1;
				}
				if (direction == antisense)
				{
					genename = itostr(GeneStart) + "-" + itostr(GeneEnd) + "C";
				}
				else
					genename = itostr(GeneStart) + "-" + itostr(GeneEnd) + "W";
				transname = genename + "_trans" + itostr(1);
				outputfile.open(outputfilename.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
				frag1 = genevertex->junctionInRange->list;
				while(frag1 != NULL)
				{
					if (frag1->junc->type != frag_junction)
					{
						outputfile << chromosome << '\t' << "Assembly\t" << "exon\t" << frag1->junc->start << '\t' << frag1->junc->end << '\t';
						outputfile << temp << '\t';
						if (direction == sense || direction == undetermined)
						{
							outputfile << "+\t";
						}
						else
						{
							outputfile << "-\t";
						}
						outputfile << ".\t";
						outputfile << "gene_id " << "\"" << genename << "\"; ";
						outputfile << "transcript_id " << "\"" << transname << "\";" << " ";
						outputfile << "copyNum " << "\"" << copyNum << "\"; " << endl;
					}
					frag1 = frag1->next;
				}

				outputfile.close();
			}
			else
			{
				outputGTF(TransList, outputfilename, chromosome, direction);
			}

			clearRedundantReads();
			clearPathclusters();

			while(TransList != NULL)
			{
				del_trans = TransList;
				TransList = del_trans->nextList;
				delete del_trans;
			}

		}
		else
		{
			ofstream outputfile;
			string genename, transname;
			double temp = 0;
			int copyNum = round(gapConstraint/exoniclength);
			if (copyNum == 0)
			{
				copyNum = 1;
			}

			if (direction == antisense)
			{
				genename = itostr(GeneStart) + "-" + itostr(GeneEnd) + "C";
			}
			else
				genename = itostr(GeneStart) + "-" + itostr(GeneEnd) + "W";
			transname = genename + "_trans" + itostr(1);
//			filename = Destination + itostr(GeneStart) + "-" + itostr(GeneEnd) + ".gtf";
			outputfile.open(outputfilename.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
			frag1 = genevertex->junctionInRange->list;
			while(frag1 != NULL)
			{
				if (frag1->junc->type != frag_junction)
				{
					outputfile << chromosome << '\t' << "Assembly\t" << "exon\t" << frag1->junc->start << '\t' << frag1->junc->end << '\t';
					outputfile << temp << '\t';
					if (direction == sense || direction == undetermined)
					{
						outputfile << "+\t";
					}
					else
					{
						outputfile << "-\t";
					}
					outputfile << ".\t";
					outputfile << "gene_id " << "\"" << genename << "\"; ";
					outputfile << "transcript_id " << "\"" << transname << "\";" << " ";
					outputfile << "copyNum " << "\"" << copyNum << "\"; " << endl;
				}
				frag1 = frag1->next;
			}

			outputfile.close();
		}

	}

	delete spliceGraph;
	delete backupList;

	while(ReadList.empty() != true)
	{
		delete ReadList[ReadList.size() - 1];
		ReadList.pop_back();
	}
	while(ReadList_cluster.empty() != true)
	{
		delete ReadList_cluster[ReadList_cluster.size() - 1];
		ReadList_cluster.pop_back();
	}
	ReadVertexIndex.clear();

	while(transStartfragList != NULL)
	{
		frag1 = transStartfragList;
		transStartfragList = frag1->next;
		delete frag1;
	}
	while(transEndfraglist != NULL)
	{
		frag1 = transEndfraglist;
		transEndfraglist = frag1->next;
		delete frag1;
	}

	return clusterSize;
}


void Reconstruct(GTvertex *rootVertex, string chromosome, string Destination, string SAMPath, int tolerance, string outputfilename, double scale)
{
	long clusterSize;
	GTvertex *curVertex;
	GTedge *curEdge;
	string filename;
	ofstream outputfile;
	string comd;
	int i;
	bool isAssembled;
	trans_direction direction;

	time_t start, stop;

//	filename = Destination + geneName + "_genelist.txt";
//	filename = Destination + itostr(rootVertex->rangeLow) + "-" + itostr(rootVertex->rangeHigh) + "_genelist.txt";
	outputfile.open(filename.c_str());
	if (rootVertex->child == NULL)
	{
		// leaf
		if (rootVertex != NULL)
		{
			time(&start);
			direction = Getdirection(rootVertex);
			if (direction == antisense)
			{
				cout << itostr(rootVertex->rangeLow) + "-" + itostr(rootVertex->rangeHigh) + "C" << endl;
			}
			else
				cout << itostr(rootVertex->rangeLow) + "-" + itostr(rootVertex->rangeHigh) + "W" << endl;
			clusterSize = Reconstruct_vertex(rootVertex, chromosome, Destination, SAMPath, isAssembled, tolerance, outputfilename, scale);

			time(&stop);
			outputfile << difftime(stop, start) << endl;
// 			if (direction == antisense)
// 			{
// 				comd = "rm -r " + SAMPath + itostr(rootVertex->rangeLow) + "-" + itostr(rootVertex->rangeHigh) + "C.txt";
// 			}
// 			else
// 				comd = "rm -r " + SAMPath + itostr(rootVertex->rangeLow) + "-" + itostr(rootVertex->rangeHigh) + "W.txt";
//			i = system(comd.c_str());
		}
	}
	else
	{
		// extend children
		if (rootVertex->level == 0 && rootVertex->child->linkedVertex->level > 0)
		{
			time(&start);
			direction = Getdirection(rootVertex);
			if (direction == antisense)
			{
				cout << itostr(rootVertex->rangeLow) + "-" + itostr(rootVertex->rangeHigh) + "C" << endl;
			}
			else
				cout << itostr(rootVertex->rangeLow) + "-" + itostr(rootVertex->rangeHigh) + "W" << endl;
			clusterSize = Reconstruct_vertex(rootVertex, chromosome, Destination, SAMPath, isAssembled, tolerance, outputfilename, scale);

			time(&stop);
			outputfile << difftime(stop, start) << endl;
// 			if (direction == antisense)
// 			{
// 				comd = "rm -r " + SAMPath + itostr(rootVertex->rangeLow) + "-" + itostr(rootVertex->rangeHigh) + "C.txt";
// 			}
// 			else
// 				comd = "rm -r " + SAMPath + itostr(rootVertex->rangeLow) + "-" + itostr(rootVertex->rangeHigh) + "W.txt";
//			i = system(comd.c_str());
		}
		else
		{
			curEdge = rootVertex->child;
			while(curEdge != NULL)
			{
				time(&start);
				curVertex = curEdge->linkedVertex;
				direction = Getdirection(curVertex);
				if (direction == antisense)
				{
					cout << itostr(curVertex->rangeLow) + "-" + itostr(curVertex->rangeHigh) + "C" << endl;
				}
				else
					cout << itostr(curVertex->rangeLow) + "-" + itostr(curVertex->rangeHigh) + "W" << endl;
				clusterSize = Reconstruct_vertex(curVertex, chromosome, Destination, SAMPath, isAssembled, tolerance, outputfilename, scale);

				time(&stop);
				outputfile << difftime(stop, start) << endl;

// 				if (direction == antisense)
// 				{
// 					comd = "rm -r " + SAMPath + itostr(curVertex->rangeLow) + "-" + itostr(curVertex->rangeHigh) + "C.txt";
// 				}
// 				else
// 					comd = "rm -r " + SAMPath + itostr(curVertex->rangeLow) + "-" + itostr(curVertex->rangeHigh) + "W.txt";
//				i = system(comd.c_str());
				curEdge = curEdge->next;
			}

		}

	}


	outputfile.close();
	return;
}



int main(int argc, char *argv[])
{
	RangeJunctionList *origList;
	string comd, filename, filepath, outputfilename;
	int i, tolerance;
	ifstream inputfile;
	double scale = 200;

	string chromosome, SAMPath, Destination, DistPath, SAMFile_prefix, geneName;
	long chromosome_length = 0;

#ifdef UNIX
	if (argc == 8)
	{
		chromosome = argv[1];
		inputPath = argv[2];
		SAMPath = argv[3];
		SAMFile_prefix = argv[4]; 
		Destination = argv[5];
		outputfilename = argv[6];
		chromosome_length = atol(argv[7]);

		if (chromosome_length <= 0)
			chromosome_length = 300000000;
		else
			chromosome_length += 100;
	}	
	else
	{
		cout << argv[0] << "\t<chromosome>\t<inputPath>\t<SAMPath>\t<Destination>\t<Round>\n";
		exit(1);
	}
#else
	if (argc == 1)
	{
		// debug
		chromosome = "chr17";
		inputPath = "";
		SAMPath = "";
		SAMFile_prefix = "chr17"; 
		Destination = "temp/";
		outputfilename = "chr17";
		chromosome_length = 489429;
	}
#endif

	cout << chromosome << endl;
	initialization();	
	inputData();

	//build tree
	gTree = new GenomeTree;
	origList = buildOrigList();

	constructGTree(origList);

	preCountGTree(gTree->root);


	/************************************************************************/
	/* assembly start                                                       */
	/* access all the genes/ASMs via gTree->root */
	/************************************************************************/

//	preprocessReadFile(gTree->root, Destination, SAMPath, SAMFile_prefix);
	tolerance = 100;
 	Reconstruct(gTree->root, chromosome, Destination, Destination, tolerance, outputfilename, scale);

	cleanAll();

	return 0;
}

