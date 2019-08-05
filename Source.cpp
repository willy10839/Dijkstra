#include <stdio.h>
#include <limits.h>
#include <iostream>
#include <math.h>
#include <stack>

using namespace std;


// A utility function to find the vertex with minimum distance value, from
// the set of vertices not yet included in shortest path tree
int minDistance(double dist[], bool sptSet[],int node)
{
	// Initialize min value
	int min = INT_MAX, min_index;

	for (int v = 0; v < node; v++)
	if (sptSet[v] == false && dist[v] <= min)
		min = dist[v], min_index = v;

	return min_index;
}

// A utility function to print the constructed distance array
void printSolution(double dist[],int dst,int node)
{
	printf("Vertex   Distance from Source\n");
	for (int i = 0; i < node; i++)
	{
		if (i==dst)
			cout << i << "		" << dist[i] << endl;
	}
}

// Funtion that implements Dijkstra's single source shortest path algorithm
// for a graph represented using adjacency matrix representation
int* dijkstra(double** graph, int src,int dst,int node)
{
	double* dist=new double[node];     // The output array.  dist[i] will hold the shortest
	// distance from src to i
	int* parent=new int[node];
	//parent[0] = src;

	bool* sptSet=new bool[node]; // sptSet[i] will true if vertex i is included in shortest
	// path tree or shortest distance from src to i is finalized

	// Initialize all distances as INFINITE and stpSet[] as false
	for (int i = 0; i < node; i++)
		dist[i] = 1000, sptSet[i] = false;


	// Distance of source vertex from itself is always 0
	dist[src] = 0;

	// Find shortest path for all vertices
	for (int count = 0; count < node - 1; count++)
	{
		// Pick the minimum distance vertex from the set of vertices not
		// yet processed. u is always equal to src in first iteration.

		int u = minDistance(dist, sptSet,node);

		// Mark the picked vertex as processed
		sptSet[u] = true;

		// Update dist value of the adjacent vertices of the picked vertex.
		for (int v = 0; v < node; v++)

			// Update dist[v] only if is not in sptSet, there is an edge from 
			// u to v, and total weight of path from src to  v through u is 
			// smaller than current value of dist[v]
		if (!sptSet[v] && graph[u][v] && dist[u] + graph[u][v] < dist[v])
		{
			dist[v] = dist[u] + graph[u][v];
			parent[v] = u;
		}
	}

	// print the constructed distance array
	printSolution(dist,dst,node);

	int pt = dst;
	cout << pt << " ";
	while (parent[pt] != src)
	{
		pt = parent[pt];
		cout << pt << " ";
	}
	pt = parent[pt];
	cout << pt << " "<<endl;
	return parent;
}

void calcost(double** graph,double** cost,double demand,int node)
{
	for (int i = 0; i < node; i++)
	{
		for (int j = 0; j < node; j++)
		{
			if (graph[i][j]!=0)
				cost[i][j] = demand / graph[i][j];
		}
	}
}

int main()
{	
	int nodenumber;
	cout << "Input the number of vertices: ";
	cin >> nodenumber;
	double** graph = new double*[nodenumber];
	for (int i = 0; i < nodenumber; i++)
		graph[i] = new double[nodenumber];

	for (int i = 0; i < nodenumber; i++)
	{
		for (int j = 0; j < nodenumber; j++)
			cin >> graph[i][j];
	}


	double** cost = new double*[nodenumber];
	for (int i = 0; i < nodenumber; i++)
		cost[i] = new double[nodenumber];

	for (int i = 0; i < nodenumber; i++)
	{
		for (int j = 0; j < nodenumber; j++)
			cost[i][j] = 0;
	}


	int x=1;
	int count=0;
	int count_satisfication=0;
	while (x!=0)
	{
		int source, dst;
		double demand;
		cout << "Please input the source: ";
		cin >> source;
		cout << "Please input the destination: ";
		cin >> dst;
		cout << "Please input the demand: ";
		cin >> demand;

		calcost(graph, cost, demand,nodenumber);

		int* parent = dijkstra(cost, source, dst,nodenumber);
		int pt = dst;
		int btk = 10000;
		while (parent[pt] != source)
		{
			if (graph[parent[pt]][pt] < btk)
				btk = graph[parent[pt]][pt];
			pt = parent[pt];
		}
		if (graph[parent[pt]][pt] < btk)
			btk = graph[parent[pt]][pt];
		pt = dst;
		if (demand < btk)
		{
			btk = demand;
			count_satisfication++;
		}
		while (parent[pt] != source)
		{
			graph[pt][parent[pt]] -= btk;
			graph[parent[pt]][pt] -= btk;
			pt = parent[pt];
		}
		graph[pt][parent[pt]] -= btk;
		graph[parent[pt]][pt] -= btk;
		count++;
		cout << "Input 0 to break other continue ";
		cin >> x;

	}
	cout << "the satisfication index is: "<<count_satisfication <<"/"<< count<<endl;

	system("pause");

	return 0;
}