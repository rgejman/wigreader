#include <iostream>
#include "wigreader.h"
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <cstdlib>
#include "functions.h"
#include <cstring>
#include <pthread.h>


using namespace std;
void *read(void *ptr );


void *read( void *ptr) {
	WigReader *reader;
	reader = (WigReader *) ptr;
	(*reader).Read();
}

void readExons(void *ptr) {
	
}

void readIntrons(void *ptr) {
	
}

void make_profile(vector<MODEL> &models, WigReader &reader, const string model_kind, bool start, int length, int blocks) {
	if(length % blocks != 0) {
		cerr << "Length must be evenly divisible by blocks" << endl;
		exit(1);
	}
	int block_len = length / blocks;
	string output_filepath = model_kind + "." + (start ? "start" : "end") + ".profile." + intToString(length) + "." + intToString(blocks) + ".txt";
	cout << "Output file: " << output_filepath << endl;
	ofstream output_file;
	output_file.open (output_filepath.c_str());
	//output_file.precision(3);
	//output_file << fixed;
	output_file << "symbol\tchr\tstart\tend\tstrand";
	for(int i=0; i < blocks; i++) {
		output_file << "\tb" + intToString(i) + "_" + model_kind;
	}
	output_file << endl;
	for(vector<MODEL>::iterator it = models.begin(); it != models.end(); ++it) {
		int block_starts[blocks];
		if((*it).strand == '+') {
			int s = start ? (*it).start : (*it).end - length;
			for(int i=0; i < blocks; i++) {
				block_starts[i] = s + (i *(block_len-1));
			}
		}
		else {
			int e = start ? (*it).end : (*it).start + length;
			for(int i=1; i <= blocks; i++) {
				block_starts[i-1] = (e+1) - i*(block_len-1);
			}
		}
		output_file << (*it).gene_id << "|" << (*it).identifier << "\t" << (*it).chr << "\t" << block_starts[0] << "\t" << block_starts[blocks-1] + block_len << "\t" << (*it).strand;
		for(int i=0; i < blocks; i++) {
			int s = block_starts[i];
			float v;
			try {
				v = reader.rpkm((*it).chr, s, s+block_len);
			}
			catch(OutOfWigRangeException outofwigrangeexception) {
				v = 0.0;
			}
			output_file << "\t" << v;
		}
		output_file << endl;
		/* std::cout << *it; ... */
	}
	output_file << endl;
	output_file.close();

}
/*
int main(int argc, char* argv[]) {
	if(argc != 5 ){
		cerr << "Usage: profileFromEnds genes.bed example1.wig 6000 300" << endl;
		exit(1);
	}
	cerr.precision(3);
	cerr << fixed;
	
	string genes_file 		= argv[1];
	string wig_file 		= argv[2];
	int total_length		= atoi(argv[3]);
	int breaks				= atoi(argv[4]);
	string prefix			= file_prefix(basename(wig_file));
	
	int one_side_length = total_length / 2;
		
	vector<MODEL> genes;
	ReadBedFileIntoModelsVector(genes_file, genes,6000);
	
	WigReader reader(wig_file.c_str());	
	reader.Read();
	
	make_profile(genes, reader, prefix, true, one_side_length, breaks);
	make_profile(genes, reader, prefix, false, one_side_length, breaks);
	
	
	return 0;
}
*/

int main(int argc, char* argv[]) {
	
	cerr.precision(3);
	cerr << fixed;
	
	const string WIG_BASE    = "/media/bigdisk/sequencing/wig";
	const string MODELS_BASE = "/home/tarakhovsky/genomics/useful_bed_files";

	const string Eugene_WT_H3K9me2_file  	= WIG_BASE + "/Eugene/Eugene_WT_H3K9me2_CD4_3x.sorted.25.wig";
	const string Eugene_KR_H3K9me2_file  	= WIG_BASE + "/Eugene/Eugene_KR_H3K9me2_CD4_3x.sorted.25.wig";
	const string Terry_WT_H3K9me2_file  	= WIG_BASE + "/Terry/Terry_WT_H3K9me2_5x.sorted.25.wig";
	const string Terry_KO_H3K9me2_file  	= WIG_BASE + "/Terry/Terry_KO_H3K9me2_5x.sorted.25.wig";
	

	const string EXONS_FILE              = MODELS_BASE + "/mm9.ensembl.exons.uniq_constitutive_or_longer.txt";
	const string INTRONS_FILE            = MODELS_BASE + "/mm9.ensembl.introns_interpolated.txt";
	//const string GENES_FILE              = MODELS_BASE + "/mm9.ensembl_with_symbols.genes.prot_coding.bed";
	const string GENES_FILE              = MODELS_BASE + "/mm9.ensembl_with_symbols.genes.prot_coding.bed";

	const string HETEROCHROMATIC_GENES_FILE  = "/media/bigdisk/projects/eugene/heterochromatic_gene_mappings.txt";
	
	cout << "Will parse HT Genes file" << endl;
	
	vector<string> heterochromatic_genes;
	readAndParseHCGenesFile(HETEROCHROMATIC_GENES_FILE, heterochromatic_genes);
	cout << "Done parsing HT Genes file" << endl;
	
	vector<MODEL> exons;
	vector<MODEL> introns;
	ReadModelFile(EXONS_FILE, exons, 200);
	removeChrs(exons, "X");
	removeIfInGeneList(exons,heterochromatic_genes);
	
	ReadModelFile(INTRONS_FILE, introns, 1000);
	removeChrs(introns, "X");
	removeIfInGeneList(introns,heterochromatic_genes);
	
	vector<MODEL> genes;
	ReadBedFileIntoModelsVector(GENES_FILE, genes,6000);
	removeChrs(genes, "X");
	vector<MODEL> HC_genes = genes;
	removeIfInGeneList(genes,heterochromatic_genes);
	removeIfNotInGeneList(HC_genes,heterochromatic_genes);
	
	pthread_t eugene_wt_thread, eugene_kr_thread, terry_wt_thread, terry_ko_thread;
	
	WigReader eugene_wt_reader(Eugene_WT_H3K9me2_file.c_str());	
	WigReader eugene_kr_reader(Eugene_KR_H3K9me2_file.c_str());
	
	WigReader terry_wt_reader(Terry_WT_H3K9me2_file.c_str());	
	WigReader terry_ko_reader(Terry_KO_H3K9me2_file.c_str());
	
	
	pthread_create( &eugene_wt_thread, NULL, read, &eugene_wt_reader);
	pthread_create( &eugene_kr_thread, NULL, read, &eugene_kr_reader);
	pthread_create( &terry_wt_thread, NULL, read, &terry_wt_reader);
	pthread_create( &terry_ko_thread, NULL, read, &terry_ko_reader);
	
	pthread_join(eugene_wt_thread, NULL);
	pthread_join(eugene_kr_thread, NULL);
	pthread_join(terry_wt_thread, NULL);
	pthread_join(terry_ko_thread, NULL);
	
	
	
	make_profile(exons, eugene_wt_reader, "Eugene_H3K9me2_WT_exons", true, 100, 10);
	make_profile(exons, eugene_wt_reader, "Eugene_H3K9me2_WT_exons", false, 100, 10);
	make_profile(exons, eugene_kr_reader, "Eugene_H3K9me2_KR_exons", true, 100, 10);
	make_profile(exons, eugene_kr_reader, "Eugene_H3K9me2_KR_exons", false, 100, 10);

	make_profile(introns, eugene_wt_reader, "Eugene_H3K9me2_WT_introns", true, 500, 50);
	make_profile(introns, eugene_wt_reader, "Eugene_H3K9me2_WT_introns", false, 500, 50);
	make_profile(introns, eugene_kr_reader, "Eugene_H3K9me2_KR_introns", true, 500, 50);
	make_profile(introns, eugene_kr_reader, "Eugene_H3K9me2_KR_introns", false, 500, 50);
	
	make_profile(genes, eugene_wt_reader, "Eugene_H3K9me2_WT_genes", true, 3000, 300);
	make_profile(genes, eugene_wt_reader, "Eugene_H3K9me2_WT_genes", false, 3000, 300);
	make_profile(genes, eugene_kr_reader, "Eugene_H3K9me2_KR_genes", true, 3000, 300);
	make_profile(genes, eugene_kr_reader, "Eugene_H3K9me2_KR_genes", false, 3000, 300);
	
	make_profile(HC_genes, eugene_wt_reader, "Eugene_H3K9me2_WT_HC_genes", true, 3000, 300);
	make_profile(HC_genes, eugene_wt_reader, "Eugene_H3K9me2_WT_HC_genes", false, 3000, 300);
	make_profile(HC_genes, eugene_kr_reader, "Eugene_H3K9me2_KR_HC_genes", true, 3000, 300);
	make_profile(HC_genes, eugene_kr_reader, "Eugene_H3K9me2_KR_HC_genes", false, 3000, 300);
	
	
	make_profile(exons, terry_wt_reader, "Terry_H3K9me2_WT_exons", true, 100, 10);
	make_profile(exons, terry_wt_reader, "Terry_H3K9me2_WT_exons", false, 100, 10);
	make_profile(exons, terry_ko_reader, "Terry_H3K9me2_KO_exons", true, 100, 10);
	make_profile(exons, terry_ko_reader, "Terry_H3K9me2_KO_exons", false, 100, 10);

	make_profile(introns, terry_wt_reader, "Terry_H3K9me2_WT_introns", true, 500, 50);
	make_profile(introns, terry_wt_reader, "Terry_H3K9me2_WT_introns", false, 500, 50);
	make_profile(introns, terry_ko_reader, "Terry_H3K9me2_KO_introns", true, 500, 50);
	make_profile(introns, terry_ko_reader, "Terry_H3K9me2_KO_introns", false, 500, 50);
	
	make_profile(genes, terry_wt_reader, "Terry_H3K9me2_WT_genes", true, 3000, 300);
	make_profile(genes, terry_wt_reader, "Terry_H3K9me2_WT_genes", false, 3000, 300);
	make_profile(genes, terry_ko_reader, "Terry_H3K9me2_KO_genes", true, 3000, 300);
	make_profile(genes, terry_ko_reader, "Terry_H3K9me2_KO_genes", false, 3000, 300);
	
	
	return 0;
}
