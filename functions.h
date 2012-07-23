#include <vector>
#include <iostream>
#include <map>
class WigReader;



struct MODEL {
	string gene_id;
	string identifier;
	string chr;
	char strand;
	long int start;
	long int end;
};

struct BEDLINE {
	string chr;
	long int start;
	long int end;
	string gene;
	char strand;
	bool operator < (const BEDLINE &s) const {
		if(chr == s.chr) {
			if(start == s.start) {
				return(end < s.end);
			}
			return(start < s.start);
		}
		return(chr < s.chr);
	}
};
string bed_prefix(string file);
string file_complete_prefix(string file);
string revcomp(const string &dna, bool reverse);
void FindAndReplace( std::string& source, const char* find, const char* replace );
void ReadModelFile(string file, vector<MODEL> &models, int min_size);
void ReadBedFileIntoModelsVector(string file, vector<MODEL> &models, int min_size);
vector<string> split(string line, const char *delim);
vector<string> tokenizeHTSymbolMapGenesFile(std::string l);
void readAndParseHCGenesFile(string file, vector<string> &hcg);
void removeChrs(vector<MODEL> &models, string chr);
void removeLtLength(vector<MODEL> &models, int length);
void removeIfInGeneList(vector<MODEL> &models, const vector<string>&gene_list);
void removeIfNotInGeneList(vector<MODEL> &models, const vector<string>&gene_list);
string basename(string file);
string file_prefix(string file);
string file_postfix(string file);
void bedlinesVectorToMapByChr(vector<BEDLINE> &bedlines, map<string,vector <BEDLINE> > &bedlinesByChr);
void readBedFile(string file, map<string, BEDLINE> &bedlines);
void readBedFile(string file, vector<BEDLINE> &bedlines);
int stringToInt(string &str);
string intToString(int number);
string enclosing_directory(string file);
