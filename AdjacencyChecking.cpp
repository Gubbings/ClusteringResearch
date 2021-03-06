#include "AdjacencyChecking.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <time.h>
#include <stdlib.h>
#include <cstring>

AdjacencyChecking::AdjacencyChecking(string graphFile, string clusterFile, string outputFile) {
	this->outputFile = outputFile;

	vertexCount = edgeCount = clusterCount = 0;

	parseInput(graphFile, clusterFile);
}

AdjacencyChecking::~AdjacencyChecking() {
	delete[] vertexNames;

	for (uint i = 0; i < vertexCount; i++) {
		delete graph[i];
		delete vertexWeights[i];
		delete clusteredVertexDegree[i];
	}

	delete[] graph;
	delete[] vertexWeights;
	delete[] adjVertsCount;
	delete[] clusteredVertexDegree;


	for (uint i = 0; i < clusterCount; i++) {
		delete clusters[i];
	}

	delete[] clusters;
	delete[] clusterSize;
}

void AdjacencyChecking::parseInput(string graphFile, string clusterFileIn) {
	//loop counter declared outside of the loop for optimization purposes
	uint i;

	list<string> *vertexNamesTemp;
	vertexNamesTemp = new list<string>;

	//maps every vertex name to an integer index that will be used to reference it
	map<string, uint> vertexMap;
	map<string, uint>::iterator vertMapIter;

	ifstream graphInStream(graphFile.c_str());
	string fileInLine;


	string vertexA, vertexB;
	uint indexVertA, indexVertB;
	float weight;

	char tempString[MAX_GRAPH_FILE_LINE_LENGTH];


	//temp lists that will store the neighbours and the weights for the vertices
	list<uint> *mappedVertListA, *mappedVertListB;
	list<float> *weightList;

	mappedVertListA = new list<uint>;
	mappedVertListB = new list<uint>;
	weightList = new list<float>;

	//track the number of neighbours a vertex has so we can transfer its adjacency list from a list to an array
	vector<uint> adjVertCountVec;

	// read the input graph
	if (graphInStream.is_open()) {
		while (!graphInStream.eof()) {

			//read a line from the graph
			graphInStream.getline(tempString, MAX_GRAPH_FILE_LINE_LENGTH, '\n');

			if (strcmp(tempString, "") == 0) {
				break;
			}

			fileInLine = tempString;

			//split the line into vertexA, vertexB and the weight - tab is the delimeter
			istringstream lineStream(tempString);

			lineStream.getline(tempString, MAX_VERTEX_NAME_LENGTH, '\t');
			vertexA = tempString;

			lineStream.getline(tempString, MAX_VERTEX_NAME_LENGTH, '\t');
			vertexB = tempString;

			lineStream >> weight;


			//map the current vertex if it has not already been mapped
			vertMapIter = vertexMap.find(vertexA);
			if (vertMapIter == vertexMap.end()) {
				indexVertA = vertexMap[vertexA] = vertexCount;
				adjVertCountVec.push_back(0);
				vertexNamesTemp->push_back(vertexA);
				vertexCount++;
			}
			else {
				indexVertA = vertMapIter->second;
			}

			vertMapIter = vertexMap.find(vertexB);
			if (vertMapIter == vertexMap.end()) {
				indexVertB = vertexMap[vertexB] = vertexCount;
				adjVertCountVec.push_back(0);
				vertexNamesTemp->push_back(vertexB);
				vertexCount++;
			}
			else {
				indexVertB = vertMapIter->second;
			}

			adjVertCountVec[indexVertA]++;
			adjVertCountVec[indexVertB]++;

			//add the vertices to the tempory lists
			//List A stores column 1 in the input file, List B stores column 2 in the input file
			mappedVertListA->push_back(indexVertA);
			mappedVertListB->push_back(indexVertB);

			//Tempory list to store all edge weights
			weightList->push_back(weight);

			edgeCount++;
		}
	}
	else {
		cout << "Unable to open graph file." << endl;
		return;
	}
	graphInStream.close();

	//convert list of names to array
	list<string>::iterator vertexNameTempItr;
	vertexNames = new string[vertexCount];

	for (i = 0, vertexNameTempItr = vertexNamesTemp->begin(); i < vertexCount; i++, vertexNameTempItr++) {
		vertexNames[i] = *vertexNameTempItr;
	}
	delete vertexNamesTemp;



	//get iterators for the temporary lists
	list<uint>::iterator mappedVertIterA, mappedVertIterB;
	list<float>::iterator weightIter;

	mappedVertIterA = mappedVertListA->begin();
	mappedVertIterB = mappedVertListB->begin();
	weightIter = weightList->begin();


	graph = new uint*[vertexCount];
	vertexWeights = new float*[vertexCount];
	adjVertsCount = new uint[vertexCount];

	for (i = 0; i < vertexCount; i++) {
		graph[i] = new uint[adjVertCountVec[i]];
		vertexWeights[i] = new float[adjVertCountVec[i]];
		adjVertsCount[i] = 0;
	}

	//convert the temp lists into 2D arrays
	for (i = 0; i < edgeCount; i++) {
		indexVertA = *mappedVertIterA;
		indexVertB = *mappedVertIterB;
		weight = *weightIter;

		/*
		The vertices are recorded in order so that:
		ColumnA [0,   3,  4]
		ColumnB [2,   1,  0]
		Weight  [.4, .1, .4]
		Nodes at the same index have an edge between them, ie 0 is connected to 2 with a weight of 0.4
		*/

		graph[indexVertA][adjVertsCount[indexVertA]] = indexVertB;
		graph[indexVertB][adjVertsCount[indexVertB]] = indexVertA;

		vertexWeights[indexVertA][adjVertsCount[indexVertA]] = weight;
		vertexWeights[indexVertB][adjVertsCount[indexVertB]] = weight;


		adjVertsCount[indexVertA]++;
		adjVertsCount[indexVertB]++;

		mappedVertIterA++;
		mappedVertIterB++;
		weightIter++;
	}


	//delete the temporary lists
	delete mappedVertListA;
	delete mappedVertListB;
	delete weightList;



	//-------------------------------------------------------read clusters------------------------------------------------
	ifstream clusterInStream(clusterFileIn.c_str());

	list<uint> *clusterTemp = nullptr;
	list<list<uint>> *clusterListTemp;

	list<list<uint>>::iterator clusterListItr;
	list<uint>::iterator clusterItr;

	clusterListTemp = new list<list<uint>>;

	//track the number of elements in a cluster this vector will be transfered to array later
	vector<uint> clusterSizeVec;


	if (clusterInStream.is_open()) {
		while (!clusterInStream.eof()) {
			//read a line from the cluster
			clusterInStream.getline(tempString, MAX_CLUSTER_FILE_LINE_LENGTH, '\n');

			if (strcmp(tempString, "") == 0) {
				break;
			}

			clusterTemp = new list<uint>;


			fileInLine = tempString;

			//get all vertices from this cluster 
			istringstream lineStream(tempString);

			int currentClusterCount = 0;

			uint vertex;
			while (!lineStream.eof()) {
				lineStream >> tempString;
				vertex = vertexMap.find(tempString)->second;
				currentClusterCount++;
				clusterTemp->push_back(vertex);
			}

			if (currentClusterCount > 1) {
				clusterSizeVec.push_back(currentClusterCount);
				clusterSizeVec[clusterCount];
				clusterListTemp->push_back(*clusterTemp);

				clusterCount++;
			}
			else {
				singletonClusters.push_back(vertex);
			}
		}
	}
	else {
		cout << "Unable to open cluster file." << endl;
		return;
	}
	clusterInStream.close();

	clusters = new uint*[clusterCount];
	clusterSize = new uint[clusterCount];
	vertexCluster = new uint[vertexCount];
	vertexPositionInCluster = new uint[vertexCount];

	for (uint i = 0; i < vertexCount; i++) {
		vertexCluster[i] = -1;
	}

	//delete the temporary lists used that store the clusters and transfer to arrays
	for (i = 0, clusterListItr = clusterListTemp->begin(); i < clusterCount; i++, clusterListItr++) {
		size_t size = clusterListItr->size();
		clusters[i] = new uint[size];

		list<uint> cluster = *clusterListItr;
		clusterItr = cluster.begin();

		for (int j = 0; j < size; j++, clusterItr++) {
			clusters[i][j] = *clusterItr;
			vertexCluster[*clusterItr] = i;
			vertexPositionInCluster[*clusterItr] = j;
		}

		clusterSize[i] = clusterSizeVec[i];
	}





	delete clusterListTemp;
	delete clusterTemp;
}


void AdjacencyChecking::beginClustering() {
	setupClusteredNeighbours();

	clock_t time = clock();

	int swapCount = 0;
	bool swapped = false;

	int checkCount = 0;

	do {
		swapped = false;

		for (uint incidentVertex = 0; incidentVertex < vertexCount; incidentVertex++) {
			for (uint adjIndex = 0; adjIndex < adjVertsCount[incidentVertex]; adjIndex++) {

				checkCount++;

				uint terminalVertex = graph[incidentVertex][adjIndex];

				uint clusterA = vertexCluster[incidentVertex];
				uint clusterB = vertexCluster[terminalVertex];

				if (terminalVertex < incidentVertex || vertexCluster[incidentVertex] == vertexCluster[terminalVertex] || clusterA == -1 || clusterB == -1) {
					continue;
				}

				float degreeAB = 0;

				int other = find(graph[incidentVertex], terminalVertex, adjVertsCount[incidentVertex]);
				if (other != -1) {
					degreeAB = vertexWeights[incidentVertex][other];
				}



				float gain = degree(incidentVertex, clusterB) + degree(terminalVertex, clusterA) - 2 * degreeAB;

				if (gain > degree(incidentVertex, clusterA) + degree(terminalVertex, clusterB)) {
					
					uint posA = vertexPositionInCluster[incidentVertex];
					uint posB = vertexPositionInCluster[terminalVertex];
					
					clusters[clusterA][posA] = terminalVertex;
					clusters[clusterB][posB] = incidentVertex;


					updateNeighbour(incidentVertex, clusterB, clusterA);
					updateNeighbour(terminalVertex, clusterA, clusterB);

					vertexCluster[incidentVertex] = clusterB;
					vertexCluster[terminalVertex] = clusterA;

					vertexPositionInCluster[incidentVertex] = posB;
					vertexPositionInCluster[terminalVertex] = posA;
					
					swapCount++;				
					swapped = true;
				}
			}
		}
		
	} while (swapped);

	cout << "Number of Swaps: " << swapCount << endl;
	cout << "Number of checks: " << checkCount << endl;

	outputClusters();
}

void AdjacencyChecking::outputClusters() {
	ofstream outStream(outputFile.c_str());

	if (outStream.is_open()) {
		for (uint i = 0; i < clusterCount; i++) {
			for (uint j = 0; j < clusterSize[i]; j++) {
				outStream << vertexNames[clusters[i][j]] << '\t';
			}
			outStream << '\n';
		}

		for (list<uint>::iterator i = singletonClusters.begin(); i != singletonClusters.end(); i++) {
			outStream << vertexNames[*i] << '\n';
		}
	}
	else {
		cout << "Unable to open output file for writing" << endl;
		return;
	}
}

/*
clusteredVertexDegree[VERTEX][CLUSTER] = support value of vertex in cluster
*/
void AdjacencyChecking::setupClusteredNeighbours() {
	clusteredVertexDegree = new float*[vertexCount];

	for (uint i = 0; i < vertexCount; i++) {
		clusteredVertexDegree[i] = new float[clusterCount] {};
	}

	for (uint i = 0; i < clusterCount; i++) {
		//check nodes that are in the current cluster 
		for (uint j = 0; j < clusterSize[i]; j++) {
			uint nodeU = clusters[i][j];

			for (uint k = 0; k < adjVertsCount[nodeU]; k++) {
				clusteredVertexDegree[graph[nodeU][k]][i] += vertexWeights[nodeU][k];
			}
		}
	}
}

void AdjacencyChecking::updateNeighbour(uint vertex, uint clusterIn, uint clusterOut) {
	for (uint i = 0; i < adjVertsCount[vertex]; i++) {
		clusteredVertexDegree[graph[vertex][i]][clusterIn] += vertexWeights[vertex][i];
		clusteredVertexDegree[graph[vertex][i]][clusterOut] -= vertexWeights[vertex][i];
	}
}

float AdjacencyChecking::degree(uint vertex, uint cluster) {
	return clusteredVertexDegree[vertex][cluster];
}


template<typename T>
int AdjacencyChecking::find(const T *arr, T key, const int length) {
	for (int i = 0; i < length; i++) {
		if (arr[i] == key) {
			return i;
		}
	}

	return -1;
}

template<typename T>
bool AdjacencyChecking::contains(const T *arr, const T key, const int length) {
	for (int i = 0; i < length; i++) {
		if (arr[i] == key) {
			return true;
		}
	}

	return false;
}


template<typename T>
bool AdjacencyChecking::contains(const vector<T> vec, const T key) {
	for (uint i = 0; i < vec.size(); i++) {
		if (vec[i] == key) {
			return true;
		}
	}

	return false;
}