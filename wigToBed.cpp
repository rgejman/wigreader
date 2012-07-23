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

string outputfileName(string wig_file) {
	string output_filename = file_complete_prefix(basename(wig_file));
	output_filename.append(".bed");
	cout << "Output: " << output_filename << endl;
	return output_filename;
}


int main(int argc, char* argv[]) {

	if(argc != 2) {
		cout << endl;
		cout << "Usage: wigToBed sample.wig => sample.bed" << cout;
		cout << endl;
		exit(1);
	}

	string wig_file			= argv[1];
	string output_file		= outputfileName(wig_file);

	WigReader reader(wig_file);
	reader.Read();

	ofstream output_stream;
	output_stream.open (output_file.c_str());
	std::vector<string> chromosomes = reader.chromosomes();
	int width = reader.getStep();

	/* Iterate over the chromosomes */
	int i = 0;
	for (std::vector<string>::iterator chr_it = chromosomes.begin(); chr_it != chromosomes.end(); ++chr_it) {
		/* Iterate over the chromosome length in steps according to $width */
		int s = reader.startForChromosome((*chr_it));
		int e = reader.endForChromosome((*chr_it));
		for(int pos = s; (pos+width) < e; pos += width) {
			i++;
			output_stream << *chr_it << "\t" << pos << "\t" << pos + width << "\tbin-" << i<< "\t" << reader.rpkm(*chr_it,pos,pos+width) << endl;
		}
	}
	output_stream.close();

	return 0;
}