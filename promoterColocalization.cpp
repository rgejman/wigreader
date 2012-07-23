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
#include <pthread.h>


using namespace std;

void *read(void *ptr );


void *read( void *ptr) {
	WigReader *reader;
	reader = (WigReader *) ptr;
	(*reader).Read();
}

int main(int argc, char* argv[]) {
	if(argc <= 6 ){
		cerr << "Usage: promoterColocalization 4000 1000 genes.bed output_basename example1.wig example2.wig => promoter_colocalization_output_basename.txt" << endl;
		exit(1);
	}
	int	promoter_width 	= atoi(argv[1]); // e.g. 1000
	int	block_width 	= atoi(argv[2]); // e.g. 1000
	
	string gene_bed_file	= argv[3];
	string output_basename	= argv[4];
	std::vector<WigReader *> wig_readers;
	std::vector<string> file_names;
	std::vector<pthread_t *> reading_threads;
	
	map<string, BEDLINE> genes;
	readBedFile(gene_bed_file,genes);
	
	ostringstream output_filename_ss;
	output_filename_ss << "promoter_colocalization_" << output_basename << ".txt";
	string output_file = output_filename_ss.str();
	
	for(int i = 5; i < argc; i++) {
		cout << "Preparing " << argv[i] << " for reading." << endl;
		string file_name_prefix = file_prefix(basename(argv[i]));
		file_names.push_back(file_name_prefix);
		WigReader *wr = new WigReader(argv[i]);
		pthread_t *thr = new pthread_t;
		pthread_create(thr, NULL, read, wr);
		reading_threads.push_back(thr);
		wig_readers.push_back(wr);
	}
	for (vector<pthread_t *>::iterator thr = reading_threads.begin(); thr!=reading_threads.end(); ++thr) {
		pthread_join(*(*thr),NULL);
	}
	
	ofstream output_stream;
	output_stream.open (output_file.c_str());
	//output_stream.precision(3);
	//output_stream << fixed;
	/* output headers */
	
	output_stream << "chr\tpos\tgene";
	for (vector<string>::iterator it = file_names.begin(); it!=file_names.end(); ++it) {
		output_stream << "\t" << *it;
	}
	output_stream << endl;
	
	/* Output data */
	
	std::vector<string> chromosomes = ((WigReader)(*wig_readers.front())).chromosomes();
	
	/* check to make sure all the step sizes are the same and
	    get the latest start position and earliest end positions for each chr */
	int step = ((WigReader)(*wig_readers.front())).getStep();
	
	for (std::vector<WigReader *>::iterator reader_it = wig_readers.begin(); reader_it!=wig_readers.end(); ++reader_it) {
		if((*(*reader_it)).getStep() != step) {
			cerr << "Step sizes were not the same!" << endl;
			exit(1);
		}
	}
	
	/* Iterate over the chromosomes */
	map<string, BEDLINE >::const_iterator itr;
	for(itr = genes.begin(); itr != genes.end(); ++itr){
		BEDLINE gene = (*itr).second;
		int s;
		int e;
		if(gene.strand == '+') {
			s = gene.start - (promoter_width / 2);
			e = s + promoter_width;
		}
		else {
			s = gene.end - (promoter_width / 2);
			e = s + promoter_width;
		}
		for(int pos = s; pos < e; pos += block_width) {
			ostringstream tmp;
			tmp << gene.chr << "\t" << pos << "\t" << gene.gene;
			/* Output a RPKM value for each wig file */
			try {
				for (std::vector<WigReader *>::iterator reader_it = wig_readers.begin(); reader_it!=wig_readers.end(); ++reader_it) {
					tmp << "\t" << (*(*reader_it)).rpkm(gene.chr,pos,pos+block_width);
				}
				
				output_stream << tmp.str() << endl;
			}
			catch(...) {
				// Doesn't matter. we simply skip this one.
			}
		}
	}
	output_stream.close();
	return 0;
}