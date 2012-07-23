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

string profileOutputFile(string wig_file,string gene_file, string analysis_type="profile") {
	string output_filename = file_prefix(basename(wig_file));
	string bed_prefix_name = bed_prefix(basename(gene_file));
	output_filename.append(".");
	output_filename.append(bed_prefix_name);
	output_filename.append(".");
	output_filename.append(analysis_type);
	output_filename.append(".wholegene.txt");
	cout << output_filename << endl;
	return output_filename;
}


string profileOutputFile(string wig_file, string gene_file ,int distance_around_tss, int blocks) {
	string output_filename = file_prefix(basename(wig_file));
	string bed_prefix_name = bed_prefix(basename(gene_file));
	output_filename.append(".");
	output_filename.append(bed_prefix_name);
	output_filename.append(".profile." + intToString(distance_around_tss) + "." + intToString(blocks) + ".txt");
	cout << output_filename << endl;
	return output_filename;
}

void profileForGene(string analysis_type, BEDLINE &gene, WigReader &reader, ofstream &output_stream) {
	ostringstream outbuff;
	//outbuff.precision(3);
	//outbuff << fixed;
	outbuff << gene.gene << "\t" << gene.chr << "\t" << gene.start << "\t" << gene.end << "\t" << gene.strand;
	try {
		float v;
		if(analysis_type == "profile") {
			v = reader.rpkm(gene.chr, gene.start, gene.end);
		}
		else if(analysis_type == "max") {
			v = reader.max_rpkm(gene.chr, gene.start, gene.end);
		}
		else if(analysis_type == "area") {
			v = reader.area(gene.chr, gene.start, gene.end);
		}
		outbuff << "\t" << v;
		outbuff << endl;
		output_stream << outbuff.str();
	}
	catch(OutOfWigRangeException outofwigrangeexception) {
		return;
	}
}


void profileForGene(BEDLINE &gene, WigReader &reader, int distance_around_tss, int blocks, int block_len, ofstream &output_stream) {
	int d 			= gene.end - gene.start;
	int to_add 		= (distance_around_tss - d)/2;
	int start		= gene.start - to_add;
	int end			= gene.end + to_add;
	ostringstream outbuff;
	//outbuff.precision(3);
	//outbuff << fixed;
	outbuff << gene.gene << "\t" << gene.chr << "\t" << start << "\t" << end << "\t" << gene.strand;
	try {
		int s;
		int e;
		int offset;
		for(int i = 0; i < blocks; i++) {
			offset = i * block_len;
			if(gene.strand == '+') {
				s = start + offset;
				e = s + block_len - 1;
			}
			else {
				e = end - offset;
				s = e - block_len + 1;
			}
			float v = reader.rpkm(gene.chr, s, e);
			outbuff << "\t" << v;
		}
		outbuff << endl;
		output_stream << outbuff.str();
	}
	catch(OutOfWigRangeException outofwigrangeexception) {
		return;
	}
}

int main(int argc, char* argv[]) {
	
	if(argc < 3 || argc > 5) {
		cout << endl;
		cout << "Usage: geneProfile area|max|profile genes.bed example.wig => Profile for entire gene" << endl;
		cout << "Usage: geneProfile tss.bed example.wig 6000 80 => Profile of TSSs +/-3kb in 80 blocks" << endl;
		cout << endl;
		exit(1);
	}
	string gene_bed_file;
	string wig_file;
	string header_postpend;
	int distance_around_tss;
	int blocks;
	int block_len;
	string output_file;
	string analysis_type;
	
	bool wholegene;
	
	if(argc > 4) {
		gene_bed_file		= argv[1]; //"/home/tarakhovsky/genomics/useful_bed_files/mm9.tss.2kb.bed"
		wig_file			= argv[2];
		header_postpend		= file_prefix(basename(wig_file));
		distance_around_tss = atoi(argv[3]);
		blocks				= atoi(argv[4]);
		output_file			= profileOutputFile(wig_file, gene_bed_file,distance_around_tss, blocks);
		if(distance_around_tss % blocks != 0) {
			cerr << "Distance around TSS must be evenly divisible by blocks." << endl;
			exit(1);
		}
		block_len			= distance_around_tss / blocks;
		wholegene			= false;
	}
	else {
		analysis_type		= argv[1];
		gene_bed_file		= argv[2]; //"/home/tarakhovsky/genomics/useful_bed_files/mm9.tss.2kb.bed"
		wig_file			= argv[3];
		header_postpend		= file_prefix(basename(wig_file));
		output_file			= profileOutputFile(wig_file, gene_bed_file,analysis_type);
		wholegene			= true;
	}
	cout << "Analysis type: " << analysis_type << endl;

	vector<BEDLINE> genes;
	readBedFile(gene_bed_file,genes);
	WigReader reader(wig_file);
	reader.Read();
	
	ofstream output_stream;
	output_stream.open (output_file.c_str());
	//output_stream.precision(3);
	//output_stream << fixed;
	output_stream << "symbol\tchr\tstart\tend\tstrand";
	if(!wholegene) {
		for(int i = 0; i < blocks; i++) {
			output_stream << "\tb" + intToString(i) + "_" + header_postpend;
		}
	}
	else {
		output_stream << "\t" + header_postpend;
	}
	output_stream << endl;
	vector<BEDLINE>::const_iterator itr;
	for(itr = genes.begin(); itr != genes.end(); ++itr){
		
		BEDLINE gene = *itr;
		if(wholegene) {
			profileForGene(analysis_type, gene, reader, output_stream);
		}
		else {
			profileForGene(gene, reader, distance_around_tss, blocks, block_len, output_stream);
		}
	}
	output_stream.close();
	
	return 0;
}