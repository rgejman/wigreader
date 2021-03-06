#include <iostream>
#include "wigreader.h"
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <cstdlib>
#include "functions.h"
#include <cstring>
#include <map>


using namespace std;

void *read(void *ptr );


void *read( void *ptr) {
	WigReader *reader;
	reader = (WigReader *) ptr;
	(*reader).Read();
}

int main(int argc, char* argv[]) {
	if(argc != 3 ){
		cerr << "Purpose: To convert un-normalized wiggle files into normalized (to RPKM) wiggle files." << endl;
		cerr << "Usage: normalizeWig 0.001203 example1.wig => example1.normalized.wig" << endl;
		cerr << "The normalization factor should be calculated as 1,000,000/num_mapped_reads." << endl;
		exit(1);
	}
	double	normalization = atof(argv[1]); // e.g. 1000
	WigReader reader(argv[2]);
	
	ostringstream output_filename_ss;
	output_filename_ss << argv[2] << ".normalized.wig";
	string output_file = output_filename_ss.str();
	cout << "Output file: " << output_file << endl;
	reader.Read();
	
	ofstream output_stream;
	output_stream.open (output_file.c_str());
	//output_stream.precision(10);
	//output_stream.setf(0,ios::floatfield); /* Not fixed. */
	
	/* Output data */
	
	map<string, vector<DATAPOINT> >::const_iterator itr;
	std::vector<string> chroms;
	for (itr = reader.data.begin(); itr != reader.data.end(); ++itr) {
		/* chr header */
		output_stream << "variableStep chrom=" << (*itr).first << " " << "span=" << reader.getStep() << endl;
		vector<DATAPOINT> v = (*itr).second;
		for (vector<DATAPOINT>::iterator dp = v.begin(); dp!=v.end(); ++dp) {
			if((*dp).rpkm > 0) {
				output_stream << (*dp).start << "\t" << ((*dp).rpkm * normalization) << endl;
			}
		}
	}
	output_stream.close();
	return 0;
}