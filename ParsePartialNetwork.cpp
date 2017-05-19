//used to find example of swap for conference presentation
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

int main(){
	ifstream graphInStream("in.txt");
	ifstream clusterIn("in.cluster");
	ofstream out("out-graph.txt");
	char tempString[1000];

	vector<string> nodes;
	if (clusterIn.is_open()) {
		while (!clusterIn.eof()) {
			
			clusterIn.getline(tempString, 1000, '\n');
			istringstream lineStream(tempString);
			
			while(!lineStream.eof()){
				lineStream.getline(tempString, 1000, '\t');
				nodes.push_back(tempString);
			}
		}
	}

	for (int i = 0; i < nodes.size(); i++) {
		cout << nodes[i] << endl;
	}

	vector<string> graphSubset;
	if (graphInStream.is_open()) {
		while (!graphInStream.eof()) {

			//read a line from the graph
			graphInStream.getline(tempString, 100000, '\n');

			if (strcmp(tempString, "") == 0) {
				break;
			}
		
			

			string line, vertexA, vertexB;
			line = tempString;

			//split the line into vertexA, vertexB and the weight - tab is the delimeter
			istringstream lineStream(tempString);

			lineStream.getline(tempString, 1000, '\t');
			vertexA = tempString;

			bool containsA = false;
			for (int i = 0; i < nodes.size(); i++) {
				if (strcmp(nodes[i].c_str(), vertexA.c_str()) == 0) {
					containsA = true;
					break;
				}
			}


			lineStream.getline(tempString, 1000, '\t');
			vertexB = tempString;

			bool containsB = false;
			for (int i = 0; i < nodes.size(); i++) {
				if (strcmp(nodes[i].c_str(), vertexB.c_str()) == 0) {
					containsB = true;
					break;
				}
			}

			

			int weight;
			lineStream >> weight;

			if (containsA && containsB) {
				graphSubset.push_back(line);
			}
		}
	}

	for (int i = 0; i < graphSubset.size(); i++) {
		out << graphSubset[i] << endl;
	}
}