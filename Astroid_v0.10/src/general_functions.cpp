/*    
 *    general_functions.cpp		
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

#include "general_functions.h"

//sort fragment
void merge_general(void** sortlist, double* sortkey, long p, long q, long r, double* mergeSort_Larray, double* mergeSort_Rarray, void** mergeSort_LorderedList, void** mergeSort_RorderedList)
{
	long n1, n2, i, j, k;

	n1 = q - p + 1;
	n2 = r - q;

	for (i = 1; i <= n1; i++)
	{
		mergeSort_Larray[i] = sortkey[p + i - 1];
		mergeSort_LorderedList[i] = sortlist[p + i - 1];
	}
	for (j = 1; j <= n2; j++)
	{
		mergeSort_Rarray[j] = sortkey[q + j];
		mergeSort_RorderedList[j] = sortlist[q + j];
	}

	mergeSort_Larray[n1 + 1] = MAX_NUMBER;
	mergeSort_Rarray[n2 + 1] = MAX_NUMBER;

	i = 1;
	j = 1;

	for (k = p; k <= r; k++)
	{
		if (mergeSort_Larray[i] <= mergeSort_Rarray[j])
		{
			sortkey[k] = mergeSort_Larray[i];
			sortlist[k] = mergeSort_LorderedList[i];

			i++;
		} 
		else
		{
			sortkey[k] = mergeSort_Rarray[j];
			sortlist[k] = mergeSort_RorderedList[j];

			j++;
		}
	}

	return;
}


void mergeSort_general(void** sortlist, double* sortkey, long sortList_size)
{
	if (sortList_size <= 0)
		return;

	//non-recursive merge sort for sorting junctions
	long m, n, i, r;
	m = 1;
	n = sortList_size;

	double* mergeSort_Larray = new double [sortList_size + 2];
	double* mergeSort_Rarray = new double [sortList_size + 2];
	void** mergeSort_LorderedList = new void* [sortList_size + 2];
	void** mergeSort_RorderedList = new void* [sortList_size + 2];

	while (m <= n)
	{
		i = 1;
		while (i <= n - m)
		{
			r = (i + 2 * m - 1) < n ? (i + 2 * m - 1) : n;
			merge_general(sortlist, sortkey, i, i + m - 1, r, mergeSort_Larray, mergeSort_Rarray, mergeSort_LorderedList, mergeSort_RorderedList);
			i = i + 2 * m;
		}

		m = m * 2;
	}

	delete [] mergeSort_Larray;
	delete [] mergeSort_Rarray;
	delete [] mergeSort_LorderedList;
	delete [] mergeSort_RorderedList;

	return;
}

double Pot (long u, long v, double *depth, double *label)
{
	return depth[u] + label[u] - label[v];
}

// Dijkstra's using non-negative edge weights (cost + potential); shortest path algorithm
bool Dijkstra( long VertexNum, long source, long sink, double **cost, long **cap, long **flownet, long *parent, long *degree, double *label, long **adjlist, bool updateflow)
{
	double *depth;
	depth = new double[VertexNum];

	for( long index = 0; index < VertexNum; index++ )
	{
		depth[index] = (MAX_NUMBER/2);
		parent[index] = -1;
	}
	depth[source] = 0;
	parent[source] = - VertexNum - 1;

	while( 1 ) 
	{
		// find u with smallest depth[u]
		long u = -1;
		long bestD = (MAX_NUMBER/2);
		for( long index = 0; index < VertexNum; index++ )
		{
			if( parent[index] < 0 && depth[index] < bestD )
				bestD = depth[u = index];
		}
		if( bestD == (MAX_NUMBER/2) ) 
			break;

		// relax edge (u,index) or (index,u) for all index;
		parent[u] = -parent[u] - 1;
		for( long index = 0; index < degree[u]; index++ )
		{
			// try undoing edge v->u      
			long v = adjlist[u][index];
			if( parent[v] >= 0 ) continue;
			if( flownet[v][u] && depth[v] > Pot(u, v, depth, label) - cost[v][u] ) 
			{
				depth[v] = Pot( u, v, depth, label) - cost[v][u];
				parent[v] = -u-1;
			}

			// try edge u->v
			if( flownet[u][v] < cap[u][v] && depth[v] > Pot(u, v, depth, label) + cost[u][v] ) 
			{
				depth[v] = Pot(u, v, depth, label) + cost[u][v];
				parent[v] = -u - 1;
			}
		}
	}

	if (updateflow == true)
	{
		for( long index = 0; index < VertexNum; index++ )
		{
			if( label[index] < (MAX_NUMBER/2) )
				label[index] += depth[index];
		}
	}

	delete [] depth;

	return parent[sink] >= 0;
}

long Min(long a, long b)
{
	return (a < b ? a : b);
}

long minCostmaxFlow( long VertexNum, long source, long sink, long sourceSup, double &fcost, long **cap, double **cost, long **flownet, long *parent, long **adjlist, long *degree, double *label, bool updateflow)
{
	for (long index = 0; index < VertexNum; index++)
	{
		label[index] = 0;
	}
	long flow = 0, supplyNum = sourceSup;
	fcost = 0;

	// repeatedly, find a cheapest path from source to sink
	while( supplyNum > 0 && Dijkstra( VertexNum, source, sink, cost, cap, flownet, parent, degree, label, adjlist, updateflow )) 
	{
		// get the bottleneck capacity
		long bot = MAX_NUMBER;
		for( long v = sink, u = parent[v]; v != source; u = parent[v = u] )
		{
			bot = Min(bot, flownet[v][u] ? flownet[v][u] : ( cap[u][v] - flownet[u][v] ));
		}
		bot = Min(bot, supplyNum);

		// update the flow network
		for( long v = sink, u = parent[v]; v != source; u = parent[v = u] )
		{
			if( flownet[v][u] ) 
			{ 
				flownet[v][u] -= bot; 
				fcost -= bot * cost[v][u]; 
			}
			else 
			{ 
				flownet[u][v] += bot; 
				fcost += bot * cost[u][v]; 
			}
		}

		flow += bot;
		supplyNum -= bot;
	}

	return flow;
}

void CopyFlownet(long VertexNum, long **flownet, long **flownet_temp)
{
	// copy the contents in flownet to flownet_temp

	for (long rowIndex = 0; rowIndex < VertexNum; rowIndex++)
	{
		for (long colIndex = 0; colIndex < VertexNum; colIndex++)
		{
			flownet_temp[rowIndex][colIndex] = flownet[rowIndex][colIndex];
		}
	}

	return;
}

bool ShortestPath_MinCostFlow(long VertexNum, long source, long sink, double **cost, long **cap, long *supply, long **flownet, vector< vector<int> > connectedIndex, vector< vector<int> > Indexassist)
{
	long flow = 0, ** adjlist, *degree, *parent, **flownet_temp, outvertexIndex, invertexIndex, sourceSup, sourceSup_min, invertexIndex_min;
	double *label, fcost, fcost_min;
	bool found_min, flowexist, foundflow;
	vector<long> supply_plus, supply_minus;

	adjlist = new long*[VertexNum]; // adjacency list
	degree = new long[VertexNum]; // degree
	label = new double[VertexNum]; // labeling function
	parent = new long[VertexNum]; // successor
	flownet_temp = new long*[VertexNum];

	for (long rowIndex = 0; rowIndex < VertexNum; rowIndex++)
	{
		degree[rowIndex] = 0;
		adjlist[rowIndex] = new long[VertexNum];
		flownet_temp[rowIndex] = new long[VertexNum];

		for (long colIndex = 0; colIndex < VertexNum; colIndex++)
		{
			flownet[rowIndex][colIndex] = 0;
			flownet_temp[rowIndex][colIndex] = 0;
			// build the adjacency list
			if( cap[rowIndex][colIndex] || cap[colIndex][rowIndex] ) 
			{
				adjlist[rowIndex][degree[rowIndex]++] = colIndex;
			}
		}
	}

	// find all vertices with supply and demand
	for (long rowIndex = 0; rowIndex < VertexNum; rowIndex++)
	{
		if (supply[rowIndex] > 0)
		{
			if (supply_plus.size() >= supply_plus.capacity())
				supply_plus.reserve(default_dataset_num + supply_plus.capacity());
			supply_plus.push_back(rowIndex);
		}
		if (supply[rowIndex] < 0)
		{
			if (supply_minus.size() >= supply_minus.capacity())
				supply_minus.reserve(default_dataset_num + supply_minus.capacity());
			supply_minus.push_back(rowIndex);
		}
	}

	flowexist = true;
	while (supply_plus.empty() != true)
	{
		outvertexIndex = supply_plus[supply_plus.size() - 1];
		fcost_min = MAX_NUMBER;
		found_min = false;
		foundflow = false;

		for (long index = 0; index < connectedIndex[outvertexIndex].size(); index++)
		{
			invertexIndex = connectedIndex[outvertexIndex][index];
			sourceSup = abs(supply[invertexIndex]);

			for (long index_temp = 0; index_temp < Indexassist[invertexIndex].size(); index_temp++)
			{
				if (Indexassist[invertexIndex][index_temp] != outvertexIndex)
				{
					if (connectedIndex[Indexassist[invertexIndex][index_temp]].size() == 1)
					{
						// invertexIndex should be used for this vertex
						sourceSup -= supply[Indexassist[invertexIndex][index_temp]];
					}
				}
			}
			sourceSup = min(supply[outvertexIndex], sourceSup);
			if (sourceSup > 0)
			{
				CopyFlownet(VertexNum, flownet, flownet_temp);
				flow = minCostmaxFlow( VertexNum, outvertexIndex, invertexIndex, sourceSup, fcost, cap, cost, flownet_temp, parent, adjlist, degree, label, false);
				if (flow > 0)
				{
					if (foundflow == false)
					{
						fcost_min = fcost/sourceSup;
						invertexIndex_min = invertexIndex;
						sourceSup_min = sourceSup;
						found_min = true;
						foundflow = true;
					}
					else
					{
						if (fcost/sourceSup < fcost_min)
						{
							fcost_min = fcost/sourceSup;
							invertexIndex_min = invertexIndex;
							sourceSup_min = sourceSup;
							found_min = true;
						}
					}
				}	
			}
			
		}	
		if (found_min == true)
		{
			sourceSup = sourceSup_min;
			flow = minCostmaxFlow( VertexNum, outvertexIndex, invertexIndex_min, sourceSup, fcost, cap, cost, flownet, parent, adjlist, degree, label, true);
			supply[outvertexIndex] -= flow;
			supply[invertexIndex_min] += flow;

			if (supply[invertexIndex_min] == 0)
			{
				for (long index = 0; index < supply_minus.size(); index++)
				{
					if (supply_minus[index] == invertexIndex_min)
					{
						supply_minus.erase(supply_minus.begin() + index);
						break;
					}
				}			
				for (long index = 0; index < Indexassist[invertexIndex_min].size(); index++)
				{
					for (long index_temp = 0; index_temp < connectedIndex[Indexassist[invertexIndex_min][index]].size(); index_temp++)
					{
						if (connectedIndex[Indexassist[invertexIndex_min][index]][index_temp] == invertexIndex_min)
						{
							connectedIndex[Indexassist[invertexIndex_min][index]].erase(connectedIndex[Indexassist[invertexIndex_min][index]].begin() + index_temp);
							break;
						}
					}
				}
			}
			if (supply[outvertexIndex] == 0)
			{
				supply_plus.pop_back();
			}
		}
		else
		{
			flowexist = false;
// 			cout << outvertexIndex << endl;
// 			for (int index = 0; index < supply_minus.size(); index++)
// 			{
// 				cout << supply_minus[index] << "\t";
// 			}
// 			cout << endl;
			break;
		}

	}

	delete [] degree;
	delete [] label;
	delete [] parent;

	for (long index = 0; index < VertexNum; index++)
	{
		delete [] adjlist[index];
		delete [] flownet_temp[index];
	}
	delete [] adjlist;
	delete [] flownet_temp;

	supply_plus.clear();
	supply_minus.clear();

	return flowexist;
}

