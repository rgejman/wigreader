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
#include <algorithm>


using namespace std;


string profileOutputFile(string wig_file) {
	string output_filename = file_prefix(basename(wig_file));
	output_filename.append(".intergenic_profile.txt");
	cout << output_filename << endl;
	return output_filename;
}


void profileForRegion(BEDLINE &gene, WigReader &reader, ofstream &output_stream) {
	ostringstream outbuff;
	//outbuff.precision(3);
	//outbuff << fixed;
	outbuff << gene.gene << "\t" << gene.chr << "\t" << gene.start << "\t" << gene.end << "\t" << gene.strand;
	try {
		float v = reader.rpkm(gene.chr, gene.start, gene.end);
		outbuff << "\t" << v;
		outbuff << endl;
		output_stream << outbuff.str();
	}
	catch(OutOfWigRangeException outofwigrangeexception) {
		return;
	}
}

int main(int argc, char* argv[]) {
	
	if(argc != 4) {
		cout << endl;
		cout << "Usage: intergenicProfile genes.bed example.wig dist_away_from_genes" << cout;
		cout << endl;
		exit(1);
	}
	
	string gene_bed_file	= argv[1]; //"/home/tarakhovsky/genomics/useful_bed_files/mm9.tss.2kb.bed"
	string wig_file			= argv[2];
	string header_postpend	= file_prefix(basename(wig_file));
	int	buffer_dist			= atoi(argv[3]);
	const int MIN_SIZE		= 1000;
	
	string output_file		= profileOutputFile(wig_file);
	
	
	vector<BEDLINE> genes;
	map<string,vector<BEDLINE> > genesMapByChr;
	
	readBedFile(gene_bed_file,genes);
	bedlinesVectorToMapByChr(genes,genesMapByChr);
	
	WigReader reader(wig_file);
	reader.Read();
	
	ofstream output_stream;
	output_stream.open (output_file.c_str());
	//output_stream.precision(3);
	//output_stream << fixed;
	output_stream << "chr\tstart\tend\trpkm" << endl;
	map<string, vector<BEDLINE>  >::const_iterator it1;
	for(it1 = genesMapByChr.begin(); it1 != genesMapByChr.end(); ++it1){
		string chr = (*it1).first;
		vector<BEDLINE> genesInThisChr = (*it1).second;
		/* sort the genes in this Chr by start and then end. */
		sort(genesInThisChr.begin(),genesInThisChr.end());
		
		/* Plan:
					1. Find the beginning of the recorded data for this chr
					2. Calculate mean(RPKM) from beginning until beginning of first gene record (minus buffer).
					3. Find end of last used gene record + buffer and record as new "beginning"
		*/
		
		int currentPos = reader.startForChromosome(chr);
		int endChr = reader.endForChromosome(chr);
		for(int currentGeneIndex=0; currentGeneIndex < genesInThisChr.size(); currentGeneIndex +=1) {
			BEDLINE g = genesInThisChr[currentGeneIndex];
			
			/* If we are past the end of the chr */
			if(currentPos >= endChr) {
				break;
			}
			
			// If our currentPos is beyond our current gene, we need to hop to the next gene.
			// This should not happen.
			if(currentPos >= (g.end + buffer_dist)) {
				currentGeneIndex += 1;
				continue;
			}
			// If currentPos is inside the current gene+buffer_dist, we need to skip to the end of this
			//   gene + buffer_dist and to the next gene.
			if(currentPos >= (g.start - buffer_dist)) {
				currentPos = g.end + buffer_dist;
				currentGeneIndex += 1;
				continue;
			}
			// Now we know we are somewhere before the next gene's buffer_dist and after the last gene's buffer_dist
			int endPos = g.start - buffer_dist; // This is where we are going until
			// If the region is smaller than the MIN_SIZE, then we skip it and move to the next gene.
			if((endPos - currentPos) < MIN_SIZE) {
				currentPos = g.end + buffer_dist;
				currentGeneIndex += 1;
				continue;
			}
			
			// Finally, we do due diligance.
			// Are we past the limit of the chr?
			// If so, let's shorten our limit
			if(endPos > endChr) {
				endPos = endChr;
				if(endPos < currentPos) {
					cerr << "ERROR: currentPos got ahead of the end of the chr." << endl;
					exit(1);
				}
			}
			float rpkm = reader.rpkm(chr, currentPos, endPos);
			// We are gt $MIN_SIZE away from the next gene's buffer_dist. Let's output the RPKM here.
			output_stream << chr << "\t" << currentPos << "\t" << endPos << "\t" << rpkm << endl;
			
			currentPos = endPos + 1;
		}
		
		
	}
	output_stream.close();
	
	return 0;
}