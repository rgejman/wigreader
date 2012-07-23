using namespace std;
#include "functions.h"
#include <string>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "wigreader.h"

string revcomp(const string &dna, bool reverse) {
	string ans = dna;
	for (int i=0; i<dna.length(); i++) {
		// find number for character
		char base;
		if (dna[i] == 'A') { base = 'T'; }
		else if (dna[i] == 'C') { base = 'G'; }
		else if (dna[i] == 'G') { base = 'C'; }
		else if (dna[i] == 'T') { base = 'A'; }
		else if (dna[i] == 'N') { base = 'N'; }
		if (reverse) {
			ans[dna.length() - 1 - i] = base;
		} else {
			ans[i] = base;
		}

	}
	return ans;
}

int stringToInt(string &str) {
	return atoi(str.c_str());
}

int stringToLongInt(string &str) {
	return atol(str.c_str());
}

string intToString(int number)
{
	stringstream ss;//create a stringstream
	ss << number;//add number to the stream
	return ss.str();//return a string with the contents of the stream
}


void FindAndReplace( std::string& source, const char* find, const char* replace )
{
	size_t findLen = strlen(find);
	size_t replaceLen = strlen(replace);
	size_t pos = 0;

//search for the next occurrence of find within source
	while ((pos = source.find( find, pos)) != std::string::npos)
	{
//replace the found string with the replacement
		source.replace( pos, findLen, replace );

//the next line keeps you from searching your replace string, 
//so your could replace "hello" with "hello world" 
//and not have it blow chunks.
		pos += replaceLen; 
	}
}

string basename(string file) {
	vector<string> tokens = split(file, "/");
	return tokens.back();
}

string enclosing_directory(string file){
	vector<string> tokens = split(file, "/");
	tokens.erase(tokens.end());
	std::stringstream ss;
	for(size_t i = 0; i < tokens.size(); ++i)
	{
		if(i == 0 && file[0] == '/') {
			ss << "/";
		}
		ss << tokens[i];
		ss << "/";
	}
	ss << "/";
	return ss.str();
}

/* removes .bed */
string bed_prefix(string file) {
	vector<string> tokens = split(file, ".bed");
	return tokens.front();
}

string file_prefix(string file) {
	vector<string> tokens = split(file, ".");
	return tokens.front();
}

/* Returns all but final extension */
/* e.g. file.5000.help.txt => file.5000.help */
string file_complete_prefix(string file) {
	vector<string> tokens = split(file, ".");
	tokens.erase(tokens.end());
	std::stringstream ss;
	for(size_t i = 0; i < tokens.size(); ++i)
	{
		if(i != 0)
			ss << ".";
		ss << tokens[i];
	}
	return ss.str();
}


string file_postfix(string file) {
	vector<string> tokens = split(file, ".");
	tokens.erase(tokens.begin());
	std::stringstream ss;
	for(size_t i = 0; i < tokens.size(); ++i)
	{
		if(i != 0)
			ss << ".";
		ss << tokens[i];
	}
	return ss.str();
}

vector<string> split(string line, const char *delim) {
	FindAndReplace(line, delim, " ");
	istringstream stm(line);
	vector<string> tokens;
	for (;;) {
		string word;
		if (!(stm >> word)) break;
		tokens.push_back(word);
	}
	return tokens;
}

void ReadBedFileIntoModelsVector(string file, vector<MODEL> &models, int min_size) {
	cout << "Reading and Parsing bed file into models file: " << file << endl;
	vector<BEDLINE> bedlines;	
	readBedFile(file, bedlines);
	for(vector<BEDLINE>::iterator it = bedlines.begin(); it != bedlines.end(); ++it) {
		MODEL m;
		vector<string> tokens = split((*it).gene, "|");
		m.gene_id = tokens[0];
		m.identifier = tokens[1];
		m.chr = (*it).chr;
		m.strand = (*it).strand;
		m.start = (*it).start;
		m.end = (*it).end;
		if((m.end - m.start) < min_size){
			continue;
		}
		models.push_back(m);
	}
	cout << "Done reading and Parsing bed file into models file: " << file << endl;

}

void ReadModelFile(string file, vector<MODEL> &models, int min_size) {
// e.g. ENSMUSG00000072821      ENSMUSE00000694008      5       +       94505411        94505464        0
	cout << "Reading and Parsing model file: " << file << endl;

	string line;
	string word;
	ifstream modelfilestream (file.c_str());
	if (modelfilestream.is_open()) {
		getline (modelfilestream,line); // header
		while ( modelfilestream.good() )
		{
			getline (modelfilestream,line);
			if(line.compare("") == 0) {
				continue;
			}
			MODEL model;
			FindAndReplace(line, "\t", " ");
			istringstream stm(line);
			stm >> model.gene_id;
			stm >> model.identifier;
			stm >> model.chr;
			stm >> word;
			model.strand		= word.c_str()[0]; // Converts "+" or "-" to '+' or '-'
			stm >> word;
			model.start			= atoi(word.c_str());
			stm >> word;
			model.end			= atoi(word.c_str());
			if((model.end - model.start) < min_size){
				continue;
			}
			models.push_back(model);
		}
		modelfilestream.close();
	}
	else {
		cerr << "Unable to open file: " << file << endl;
		exit(1);
	}
	cout << "Done Reading and Parsing model file: " << file << endl;
}


vector<string> tokenizeHTSymbolMapGenesFile(std::string l) 
{
	FindAndReplace(l,"\t"," ");
	FindAndReplace(l,","," ");

//std::replace(l.begin(), l.end(), '\t', ' ');
//std::replace(l.begin(), l.end(), ',', ' ');
	istringstream stm(l);
	vector<string> tokens;
	for (;;) {
		string word;
		if (!(stm >> word)) break;
		tokens.push_back(word);
	}
	return tokens;
}

void readAndParseHCGenesFile(string file, vector<string> &hcg) {
	cout << "Reading and Parsing HT Genes file" << endl;

	string line;
	ifstream genefilestream (file.c_str());
	if (genefilestream.is_open()) {
		while ( genefilestream.good() )
		{
			getline (genefilestream,line);
			vector<string> tokens = tokenizeHTSymbolMapGenesFile(line);
			tokens.erase (tokens.begin());
			hcg.insert(hcg.end(),tokens.begin(), tokens.end());
		}
		genefilestream.close();
	}
	else {
		cerr << "Unable to open file: " << file << endl;
		exit(1);
	}
}

void removeChrs(vector<MODEL> &models, string chr) {
	for (int i=0; i<models.size(); i++) {
		if(models[i].chr == chr) {
		//cout << models[i].chr << ' ' << chr << endl;
			models.erase(models.begin() + i);
			i--;
		}
	}
}

void removeLtLength(vector<MODEL> &models, int length) {
	for (int i=0; i<models.size(); i++) {
		if((models[i].end - models[i].start) < length) {
		//cout << "Removing: " << models[i].end - models[i].start << endl;
			models.erase(models.begin() + i);
			i--;
		}
	}
}

void removeIfNotInGeneList(vector<MODEL> &models, const vector<string>&gene_list) {
	for (int i=0; i<models.size(); i++) {
		bool del = true;
		for (int j=0; j<gene_list.size(); j++) {
			if(models[i].gene_id == gene_list[j]) {
				del = false;
				break;
			}
		}
		if(del) {
			models.erase(models.begin() + i);
			i--;
		}
	}
}

void removeIfInGeneList(vector<MODEL> &models, const vector<string>&gene_list) {
	for (int i=0; i<models.size(); i++) {
		for (int j=0; j<gene_list.size(); j++) {
			if(models[i].gene_id == gene_list[j]) {
			//cout << "Removing: " << s << endl;
				models.erase(models.begin() + i);
				i--;
				break;
			}
		}
	}
}

void bedlinesVectorToMapByChr(vector<BEDLINE> &bedlines, map<string,vector<BEDLINE> > &bedlinesByChr) {
	for(vector<BEDLINE>::iterator it = bedlines.begin(); it != bedlines.end(); ++it) {
		bedlinesByChr[(*it).chr].push_back(*it);
	}
}

void readBedFile(string file, map<string, BEDLINE> &bedlines) {
	cout << "Reading Bed File (map)" << endl;

	string line;
	ifstream bedfilestream (file.c_str());
	if (bedfilestream.is_open()) {
		while ( bedfilestream.good() )
		{
			getline (bedfilestream,line);
			if(line == "") {
				continue;
			}
			vector<string> tokens = split(line,"\t");
			BEDLINE bedline;
			bedline.chr		= tokens[0];
			bedline.start	= stringToLongInt(tokens[1]);
			bedline.end		= stringToLongInt(tokens[2]);
			bedline.gene	= tokens[3];
			// then there is a "value" field that we skip.
			if(tokens.size() >= 6) {
				bedline.strand	= tokens[5].c_str()[0];
			}
			bedlines[bedline.gene] = bedline;
		}
		bedfilestream.close();
	}
	else {
		cerr << "Unable to open file: " << file << endl;
		exit(1);
	}
	cout << "Done reading Bed File" << endl;
}

void readBedFile(string file, vector<BEDLINE> &bedlines) {
	cout << "Reading Bed File (vector)" << endl;

	string line;
	ifstream bedfilestream (file.c_str());
	if (bedfilestream.is_open()) {
		while ( bedfilestream.good() )
		{
			getline (bedfilestream,line);
			if(line == "") {
				continue;
			}
			vector<string> tokens = split(line,"\t");
			BEDLINE bedline;
			bedline.chr		= tokens[0];
			bedline.start	= stringToLongInt(tokens[1]);
			bedline.end		= stringToLongInt(tokens[2]);
			bedline.gene	= tokens[3];
			// then there is a "value" field that we skip.
			if(tokens.size() >= 6) {
				bedline.strand	= tokens[5].c_str()[0];
			}
			bedlines.push_back(bedline);
		}
		bedfilestream.close();
	}
	else {
		cerr << "Unable to open file: " << file << endl;
		exit(1);
	}
	cout << "Done reading Bed File" << endl;
}