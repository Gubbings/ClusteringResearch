#include "VariableDepthLSA.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <time.h>
#include <stdlib.h>
#include <cstring>
#include <algorithm>

using namespace std;

struct candidateSwap {	
	candidateSwap() {
		vertexIndexA = -1;
		vertexIndexB = -1;
		gain = 0;
	}
		
	uint vertexIndexA;
	uint vertexIndexB;
	float gain;
};

struct concreteSwap {
	concreteSwap() {
		vertexA = -1;
		vertexB = -1;
		gain = 0;
	}

	uint vertexA;
	uint vertexB;
	float gain;
};

struct vertex {
	uint vertexID;
	uint clusterIndex;
};

bool compareSwaps(candidateSwap a, candidateSwap b) {
	return a.gain > b.gain;
}

VariableDepthLSA::VariableDepthLSA(string graphFile, string clusterFile, string outputFile) {
	this->outputFile = outputFile;

	vertexCount = edgeCount = clusterCount = 0;

	parseInput(graphFile, clusterFile);
}

VariableDepthLSA::~VariableDepthLSA() {
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

	delete[] vertexCluster;
	delete[] clusters;
	delete[] clusterSize;
}

void VariableDepthLSA::parseInput(string graphFile, string clusterFileIn) {
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
	vertexClusterIndex = new uint[vertexCount];

	//delete the temporary lists used that store the clusters and transfer to arrays
	for (i = 0, clusterListItr = clusterListTemp->begin(); i < clusterCount; i++, clusterListItr++) {
		size_t size = clusterListItr->size();
		clusters[i] = new uint[size];

		list<uint> cluster = *clusterListItr;
		clusterItr = cluster.begin();

		for (int j = 0; j < size; j++, clusterItr++) {
			clusters[i][j] = *clusterItr;
			vertexCluster[*clusterItr] = i;
			vertexClusterIndex[*clusterItr] = j;
		}

		clusterSize[i] = clusterSizeVec[i];
	}





	delete clusterListTemp;
	delete clusterTemp;
}


void VariableDepthLSA::beginClustering() {
	clock_t time = clock();
	int swapCount = 0;

	//initialize support values
	clusteredVertexDegree = new float*[vertexCount];

	for (uint i = 0; i < vertexCount; i++) {
		clusteredVertexDegree[i] = new float[clusterCount] {};
	}

/*
	int currentIterationSwapCount = 0;

	do {
		currentIterationSwapCount = 0;
*/

		//arbitrary ordering and selection of pairs of clusters	
		for (uint clusterA = 0; clusterA < clusterCount - 1; clusterA++) {

			for (uint clusterB = clusterA + 1; clusterB < clusterCount; clusterB++) {

				list<concreteSwap> gainSequence;

				//make copies of clusterA and clusterB where we will swap to find the sequence			
				uint *clusterACopy = new uint[clusterSize[clusterA]];
				memcpy(clusterACopy, clusters[clusterA], clusterSize[clusterA] * sizeof(uint));

				uint *clusterBCopy = new uint[clusterSize[clusterB]];
				memcpy(clusterBCopy, clusters[clusterB], clusterSize[clusterB] * sizeof(uint));

				//calculate support for all vertices in the cluster pair
				setupClusteredNeighbours(clusterA);
				setupClusteredNeighbours(clusterB);

				//create a copy of the support matrix to use while computing the gain sequence
				float *clusterASupportCopy = new float[vertexCount];
				float *clusterBSupportCopy = new float[vertexCount];

				//copy the clusterA and clusterB columns of the clusteredVertexDegree array
				for (uint i = 0; i < vertexCount; i++) {
					clusterASupportCopy[i] = clusteredVertexDegree[i][clusterA];
					clusterBSupportCopy[i] = clusteredVertexDegree[i][clusterB];
				}

				int availableA = clusterSize[clusterA];
				int availableB = clusterSize[clusterB];
				int candidateSwapCount = availableA * availableB;

				candidateSwap *swaps = new candidateSwap[candidateSwapCount];

				/*
				-------------------------------------- Search for swaps until no more are possible ---------------------------------------------------
					* for a pair of uniform clusters this would mean net gain is 0 this is not always true for non-uniform clusters
				*/
				while (availableA > 0 && availableB > 0) {

					for (uint i = 0; i < availableA; i++) {
						uint vertexA = clusterACopy[i];

						for (uint j = 0; j < availableB; j++) {
							uint vertexB = clusterBCopy[j];

							float degreeAB = 0;
							int other = find(graph[vertexA], vertexB, adjVertsCount[vertexA]);
							if (other != -1) {
								degreeAB = vertexWeights[vertexA][other];
							}

							float costA = clusterBSupportCopy[vertexA] - clusterASupportCopy[vertexA];
							float costB = clusterASupportCopy[vertexB] - clusterBSupportCopy[vertexB];

							swaps[i * availableB + j].vertexIndexA = i;
							swaps[i * availableB + j].vertexIndexB = j;
							swaps[i * availableB + j].gain = costA + costB - 2 * degreeAB;
						}
					}

					//select the best swap - the swap array will be partiall overriden in each iteration of the above loops up to at most availableA * availableB
					sort(swaps, swaps + availableA * availableB, compareSwaps);

					concreteSwap bestSwap;
					bestSwap.vertexA = clusterACopy[swaps[0].vertexIndexA];
					bestSwap.vertexB = clusterBCopy[swaps[0].vertexIndexB];
					bestSwap.gain = swaps[0].gain;

					gainSequence.push_back(bestSwap);

					/*
					-------------------------------------- Swap the pair with best gain ---------------------------------------------------
					*/
					//perform the swap
					uint temp = clusterACopy[swaps[0].vertexIndexA];
					clusterACopy[swaps[0].vertexIndexA] = clusterBCopy[swaps[0].vertexIndexB];
					clusterBCopy[swaps[0].vertexIndexB] = temp;

					//update cluster index for these vertices
					vertexCluster[swaps[0].vertexIndexA] = clusterB;
					vertexCluster[swaps[0].vertexIndexB] = clusterA;


					//move non-available elements to the back of the cluster where we wont access them
					temp = clusterACopy[availableA - 1];
					clusterACopy[availableA - 1] = clusterACopy[swaps[0].vertexIndexA];
					clusterACopy[swaps[0].vertexIndexA] = temp;
					availableA--;

					temp = clusterBCopy[availableB - 1];
					clusterBCopy[availableB - 1] = clusterBCopy[swaps[0].vertexIndexB];
					clusterBCopy[swaps[0].vertexIndexB] = temp;
					availableB--;


					/*
						-------------------------------------- Update costs ---------------------------------------------------
					*/

					//update costs of vertices adjacent to vertexA that were already in clusterA - their cost increases because A is leaving
					for (int i = 0; i < adjVertsCount[swaps[0].vertexIndexA]; i++) {
						//consider only adjacent vertices that are present in clusterA
						if (vertexCluster[graph[swaps[0].vertexIndexA][i]] == clusterA) {
							clusterASupportCopy[graph[swaps[0].vertexIndexA][i]] += 2 * vertexWeights[swaps[0].vertexIndexA][i];
						}
					}

					//update costs of vertices adjacent to vertexB that were already in clusterA - their cost decreses because B is entering
					for (int i = 0; i < adjVertsCount[swaps[0].vertexIndexB]; i++) {
						//consider only adjacent vertices that are present in clusterA
						if (vertexCluster[graph[swaps[0].vertexIndexB][i]] == clusterA) {
							clusterASupportCopy[graph[swaps[0].vertexIndexB][i]] -= 2 * vertexWeights[swaps[0].vertexIndexB][i];
						}
					}

					//update costs of vertices adjacent to vertexB that were already in clusterB - their cost increases because B is leaving
					for (int i = 0; i < adjVertsCount[swaps[0].vertexIndexB]; i++) {
						//consider only adjacent vertices that are present in clusterB
						if (vertexCluster[graph[swaps[0].vertexIndexB][i]] == clusterB) {
							clusterASupportCopy[graph[swaps[0].vertexIndexB][i]] += 2 * vertexWeights[swaps[0].vertexIndexB][i];
						}
					}

					//update costs of vertices adjacent to vertexA that were already in clusterB - their cost decreses because A is entering
					for (int i = 0; i < adjVertsCount[swaps[0].vertexIndexA]; i++) {
						//consider only adjacent vertices that are present in clusterB
						if (vertexCluster[graph[swaps[0].vertexIndexA][i]] == clusterB) {
							clusterASupportCopy[graph[swaps[0].vertexIndexA][i]] -= 2 * vertexWeights[swaps[0].vertexIndexA][i];
						}
					}
				}

				delete[] swaps;
				delete[] clusterACopy;
				delete[] clusterBCopy;
				delete[] clusterASupportCopy;
				delete[] clusterBSupportCopy;

				//search for best net gain
				float totalGain = 0;
				float bestGain = 0;
				int bestGainIndex = -1;

				int k = 0;
				for (list<concreteSwap>::iterator i = gainSequence.begin(); i != gainSequence.end(); i++, k++) {
					totalGain += (*i).gain;
					if (totalGain > bestGain) {
						bestGain = totalGain;
						bestGainIndex = k;
					}
				}
				
				//debugging
				if (bestGainIndex > 0) {
					cout << "Swap Sequence Size: " << bestGainIndex << endl;
				}

				list<concreteSwap>::iterator i = gainSequence.begin();
				for (k = 0; k <= bestGainIndex; k++, i++) {
					//perform the swap				
					int indexA = vertexClusterIndex[(*i).vertexA];
					int indexB = vertexClusterIndex[(*i).vertexB];

					uint temp = clusters[clusterA][indexA];
					clusters[clusterA][indexA] = clusters[clusterB][indexB];
					clusters[clusterB][indexB] = temp;

					//update cluster for these vertices
					vertexCluster[indexA] = clusterB;
					vertexCluster[indexB] = clusterA;

					//update cluster index for these vertices
					vertexClusterIndex[(*i).vertexA] = indexB;
					vertexClusterIndex[(*i).vertexB] = indexA;

					//currentIterationSwapCount++;
					swapCount++;
				}
			}
		}
/*
		cout << "Current Iteration Swap Count = " << currentIterationSwapCount << endl;
	} while (currentIterationSwapCount > 0);
*/

	time = clock() - time;
	cout << "Clock ticks taken = " << time << endl;
	cout << "Swap Count = " << swapCount << endl;

	outputClusters();
}

void VariableDepthLSA::outputClusters() {
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
void VariableDepthLSA::setupClusteredNeighbours(uint clusterIndex) {

	for (uint j = 0; j < clusterSize[clusterIndex]; j++) {
		uint nodeU = clusters[clusterIndex][j];

		for (uint k = 0; k < adjVertsCount[nodeU]; k++) {
			clusteredVertexDegree[graph[nodeU][k]][clusterIndex] += vertexWeights[nodeU][k];
		}
	}

}

void VariableDepthLSA::updateNeighbour(uint vertex, uint clusterIn, uint clusterOut) {
	for (uint i = 0; i < adjVertsCount[vertex]; i++) {
		clusteredVertexDegree[graph[vertex][i]][clusterIn] += vertexWeights[vertex][i];
		clusteredVertexDegree[graph[vertex][i]][clusterOut] -= vertexWeights[vertex][i];
	}
}


/*
float VariableDepthLSA::externalCost(uint vertex, uint clusterIndexB) {
	return clusteredVertexDegree[vertex][clusterIndexB];
}

float VariableDepthLSA::internalCost(uint vertex, uint clusterIndexA) {
	return clusteredVertexDegree[vertex][clusterIndexA];
}

//TODO
float VariableDepthLSA::cost(uint vertex, uint clusterIndexA, uint clusterIndexB) {
	//convert this to not use functions later this is just for semantic description
	return externalCost(vertex, clusterIndexB) - internalCost(vertex, clusterIndexA);
}
*/

template<typename T>
int VariableDepthLSA::find(const T *arr, T key, const int length) {
	for (int i = 0; i < length; i++) {
		if (arr[i] == key) {
			return i;
		}
	}

	return -1;
}

template<typename T>
bool VariableDepthLSA::contains(const T *arr, const T key, const int length) {
	for (int i = 0; i < length; i++) {
		if (arr[i] == key) {
			return true;
		}
	}

	return false;
}


template<typename T>
bool VariableDepthLSA::contains(const vector<T> vec, const T key) {
	for (uint i = 0; i < vec.size(); i++) {
		if (vec[i] == key) {
			return true;
		}
	}

	return false;
}