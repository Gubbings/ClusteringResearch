#include <iostream>
#include <string>
#include <stdlib.h>
#include "LocalSearchClustering.h"
#include "LocalSearchClustering2.h"
#include "VariableDepthLSA.h"
#include "EjectionChainMethod.h"
#include "AdjacencyChecking.h"

using namespace std;

int main(int argc, char *argv[]) {

	string in = "in.txt";
	string clusters = "in.cluster";
	string out = "out.txt";
	int mode = 3;

	if (argc != 5) {
		cout << "Arguments not properly specified using defaults. " << endl;
		cout << "Expected arguments 1: graph file, 2: cluster file, 3: output file" << endl;
		cout << "Default settings: \ngraph file = in.txt, \nclusters file = in.clsuter, \noutput file = out.txt" << endl;
	}
	else {
		in = argv[1];
		clusters = argv[2];
		out = argv[3];		
	}
	
	if (mode == 0) {
		LocalSearchClustering localSearchAlg(in, clusters, out);
		localSearchAlg.beginClustering();
	}
	else if(mode == 1){
		VariableDepthLSA variableDepth(in, clusters, out);
		variableDepth.beginClustering();		
	}
	else if (mode == 2) {
		EjectionChainMethod ejectionChainMethod(in, clusters, out);
		ejectionChainMethod.beginClustering();
	}
	else if (mode == 3) {
		AdjacencyChecking adjChecking(in, clusters, out);
		adjChecking.beginClustering();
	}
	

	return 0;
}

