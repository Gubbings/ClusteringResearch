#pragma once
#include <string>
#include <vector>
#include <list>

#define MAX_CLUSTER_FILE_LINE_LENGTH 100000
#define MAX_GRAPH_FILE_LINE_LENGTH 100000
#define MAX_VERTEX_NAME_LENGTH 1000

typedef unsigned int     uint;

using namespace std;

class AdjacencyChecking {
public:
	AdjacencyChecking(string graphFile, string clusterFile, string outputFile);
	~AdjacencyChecking();
	void beginClustering();

private:
	string *vertexNames;
	uint **graph, *adjVertsCount, *clusterSize, *vertexCluster, *vertexPositionInCluster;
	float **vertexWeights;

	uint vertexCount, edgeCount;

	uint **clusters;
	uint clusterCount;

	vector<uint> **clusteredNeighbour;
	float **clusteredVertexDegree;

	list<uint> singletonClusters;

	string outputFile;
	uint swapCount;

	float degree(uint vertex, uint cluster);
	void parseInput(string graphFile, string clusterFileIn);
	void outputClusters();
	void setupClusteredNeighbours();
	void updateNeighbour(uint vertex, uint clusterIn, uint clusterOut);

	template<typename T>
	bool contains(const T *arr, const T key, const int length);
	template<typename T>
	bool contains(const vector<T> vec, const T key);
	template<typename T>
	int find(const T *arr, const T key, const int length);

};