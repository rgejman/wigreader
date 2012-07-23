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

/* 
	1. Takes an arbitrary # of wig files.
	2. Look at wig files, starting from lowest position in each chr (across all files)
	3. Take $width blocks and compare RPKM at that location.
	4. Output those.
*/

void *read(void *ptr );


void *read( void *ptr) {
	WigReader *reader;
	reader = (WigReader *) ptr;
	(*reader).Read();
}

int main(int argc, char* argv[]) {
	if(argc <=4 ){
		cerr << "Usage: genomicColocalization 1000 output_basename example1.wig example2.wig => genomic_colocalization_output_basename.txt" << endl;
		exit(1);
	}
	int	width = atoi(argv[1]); // e.g. 1000
	string output_basename = argv[2];
	std::vector<WigReader *> wig_readers;
	std::vector<string> file_names;
	std::vector<pthread_t *> reading_threads;
	
	ostringstream output_filename_ss;
	output_filename_ss << "genomic_colocalization_" << output_basename << ".txt";
	string output_file = output_filename_ss.str();
	
	for(int i = 3; i < argc; i++) {
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
	/* output headers */
	
	output_stream << "chr\tpos";
	for (vector<string>::iterator it = file_names.begin(); it!=file_names.end(); ++it) {
		output_stream << "\t" << *it;
	}
	output_stream << endl;
	
	/* Output data */
	
	std::vector<string> chromosomes = ((WigReader)(*wig_readers.front())).chromosomes();
	
	/* check to make sure all the step sizes are the same and
	    get the latest start position and earliest end positions for each chr */
	int step = ((WigReader)(*wig_readers.front())).getStep();
	std::map<string,int> latestStarts;
	std::map<string,int> earliestEnds;
	
	for (std::vector<WigReader *>::iterator reader_it = wig_readers.begin(); reader_it!=wig_readers.end(); ++reader_it) {
		if((*(*reader_it)).getStep() != step) {
			cerr << "Step sizes were not the same!" << endl;
			exit(1);
		}
		;
		for (std::vector<string>::iterator chr_it = chromosomes.begin(); chr_it != chromosomes.end(); ++chr_it) {
			int s = (*(*reader_it)).startForChromosome((*chr_it));
			int e = (*(*reader_it)).endForChromosome((*chr_it));
			
			if(s == 0 || e == 0) {
				cerr << "Chromosome " << *chr_it << " does not exist in all wigfiles." << endl;
				exit(1);
			}
			if(latestStarts[*chr_it] == 0) {
				latestStarts[*chr_it] = s;
			}
			else if(latestStarts[*chr_it] < s) {
				latestStarts[*chr_it] = s;
			}
			
			if(earliestEnds[*chr_it] == 0) {
				earliestEnds[*chr_it] = e;
			}
			else if(earliestEnds[*chr_it] > e) {
				earliestEnds[*chr_it] = e;
			}
		}
	}
	
	/* Iterate over the chromosomes */
	for (std::vector<string>::iterator chr_it = chromosomes.begin(); chr_it != chromosomes.end(); ++chr_it) {
		/* Iterate over the chromosome length in steps according to $width */
		int s = latestStarts[*chr_it];
		int e = earliestEnds[*chr_it];
		for(int pos = s; (pos+width) < e; pos += width) {
			output_stream << *chr_it << "\t" << pos;
			/* Output a RPKM value for each wig file */
			for (std::vector<WigReader *>::iterator reader_it = wig_readers.begin(); reader_it!=wig_readers.end(); ++reader_it) {
				output_stream << "\t" << (*(*reader_it)).rpkm(*chr_it,pos,pos+width);
			}
			output_stream << endl;
		}
	}
	output_stream.close();
	return 0;
}