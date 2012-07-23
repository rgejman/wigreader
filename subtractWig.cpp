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
	if(argc != 4 ){
		cerr << "Usage: subractWig output_basename chip.wig input.wig => output_basename.wig" << endl;
		exit(1);
	}
	string output_basename = argv[1];
	std::vector<WigReader *> wig_readers;
	std::vector<string> file_names;
	std::vector<pthread_t *> reading_threads;
	
	ostringstream output_filename_ss;
	output_filename_ss << output_basename << ".wig";
	string output_file = output_filename_ss.str();
	
	for(int i = 2; i < argc; i++) {
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
	
	WigReader *chip = wig_readers.front();
	WigReader *input = wig_readers.back();
	
	/* Iterate over the chromosomes */
	for (std::vector<string>::iterator chr_it = chromosomes.begin(); chr_it != chromosomes.end(); ++chr_it) {
		/* Iterate over the chromosome length in steps according to $width */
		output_stream << "variableStep chrom=" << *chr_it << " " << "span=" << step << endl;
		int s = latestStarts[*chr_it];
		int e = earliestEnds[*chr_it];
		float v;
		for(int pos = s; (pos+step) < e; pos += step) {
			v = (*chip).rpkm(*chr_it, pos,pos+step) - (*input).rpkm(*chr_it, pos,pos+step);
			output_stream << pos << "\t" << v << endl;
		}
	}
	output_stream.close();
	return 0;
}