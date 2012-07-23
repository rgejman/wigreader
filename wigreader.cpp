using namespace std;
#include <cstdlib>
#include <iostream>
#include "wigreader.h"
#include <fstream>
#include <string>
#include <sstream>
#include "functions.h"
#include <math.h>


bool data_point_compare (DATAPOINT a, DATAPOINT b)
{
	return a.start < b.start;
}

std::vector<std::string> parseHeader(std::string l) 
{
	FindAndReplace(l,"="," ");
	std::istringstream stm(l);
	std::vector<std::string> tokens;
	for (;;) {
		std::string word;
		if (!(stm >> word)) break;
		tokens.push_back(word);
	}
	return tokens;
}

WigReader::WigReader(string file) {
	init(file.c_str());
}


WigReader::WigReader(const char* file) {
	init(file);
}

void WigReader::init(const char* file) {
	wigfile 				= file;
	step 					= -1;
	last_pos 				= -1;
	skip_until_next_header 	= false;
}

void WigReader::Read() {
	string line;
	ifstream wigfilestream (wigfile);
	cout << "Will open wigreader on file: " << wigfile << endl;
	if (wigfilestream.is_open()) {
		cout << "Opened WigReader on: " << wigfile << endl;
		while ( wigfilestream.good() )
		{
			getline (wigfilestream,line);
			ProcessLine(line);
		}
		wigfilestream.close();
	}
	else {
		cerr << "Unable to open file: " << wigfile << endl;
		exit(1);
	}
	map<string, vector<DATAPOINT> >::const_iterator itr;

	for(itr = data.begin(); itr != data.end(); ++itr){
		string key = (*itr).first;
		DATAPOINT last_dp = data[key].back();		
		ends[key] = last_dp.start;
	}
}

void WigReader::ProcessLine(const string &line) {
	//cout << "ProcessLine: " << line << endl;
	if(line[0] == 'v'){
		ProcessHeader(line);
	}
	else if(!skip_until_next_header) {
		ProcessDataLine(line);	
	}
}

void WigReader::ProcessHeader(const string &header){
	// e.g. variableStep chrom=chr1 step=25
	skip_until_next_header = false;
	vector<string> tokens = parseHeader(header);
	tokens.erase(tokens.begin()); // remove "variableStep"
	chr = tokens[1];
	if(data.find(chr) != data.end()) {
		cerr << "Skipping some lines from " << chr << "." << endl;
		skip_until_next_header = true;
		return;
	}

	int step_l = atoi(tokens[4].c_str());
	if(step == -1) {
		step = step_l;
		cout << "Step: " << step << endl;
	}
	else if(step != step_l) {
		cerr <<  "Error: Step was " << step << " until " << chr << " which had step=" << step_l << endl;
		cerr <<  "Skipping " << chr << endl;
		skip_until_next_header = true;
		return;
		//exit(1);
	}
	last_pos = -1;
	current_data_points_vector = &data[chr]; // initializes the vector for this chr and assigns it to the current_data_points_vector var.
	cout << "Reading " << chr << " with step size " << step << endl;

}

void WigReader::ProcessDataLine(const string &line) {
	// e.g. 6748876	2.0
	DATAPOINT dp;
	if(sscanf(line.c_str(), "%ld\t%f", &dp.start, &dp.rpkm)  == EOF ) {
		if(line.compare("") != 0) {
			cerr << "Error reading data line: " << line << endl;
			exit(1);
		}
		else {
			return;
		}
	}
	/* Set the starting position for this chr */
	if(last_pos == -1) {
		starts[chr] = dp.start;
	}
	while(last_pos >= 0 && last_pos < (dp.start - step)) {
		last_pos += step;
		DATAPOINT fillin_dp;
		fillin_dp.start = last_pos;
		fillin_dp.rpkm = 0.0;
		(*current_data_points_vector).push_back(fillin_dp);
		//cout << "Fillin: " << ": " << fillin_dp.start << ' ' << fillin_dp.rpkm << endl;
	}
	(*current_data_points_vector).push_back(dp);
	last_pos = dp.start;
	//cout << "DataLine:" << line << ": " << dp.start << ' ' << dp.rpkm << endl;
}

string WigReader::locate_chr_key(const string chr) {
	string key = chr;
	map<string,vector<DATAPOINT> >::iterator chr_data = data.find(key);
	if(chr_data == data.end()) {
		string new_key;
		if(key.compare("MT") == 0) {
			new_key = "chrM";
		}
		else {
			new_key = chr;
			new_key.insert(0,"chr");
		}
		chr_data = data.find(new_key);
		if(chr_data == data.end()) {
			cerr << "Neither " << chr << " nor " << new_key << " are found in the wig file." << endl;
			throw OutOfWigRangeException();
		}
		else {
			key = new_key;
		}
	}
	return(key);
}

float WigReader::max_rpkm(const string chr, long int a, long int b){
	string key = locate_chr_key(chr);
	map<string,vector<DATAPOINT> >::iterator chr_data = data.find(key);
	if(a == b) {
		return(0);
	}
	long int start 	= starts[key];
	long int end	= ends[key];
	if(a < b) {
		if(a < start) {
			return(0);
		}
		if(b > end) {
			return(0);
		}
		return calculate_max_rpkm(chr_data,a,b,start);
	}
	else {
		if(b < start) {
			return(0);
		}
		if(a > end) {
			return(0);
		}
		try {
			return calculate_max_rpkm(chr_data,b,a,start);
		}
		catch(...) {
			cerr << "Error occurred in calculate_rpkm: " << a << " " << " " << b << " " << start << endl;
		}
	}
	
}
float WigReader::area(const string chr, long int a, long int b) {
	string key = locate_chr_key(chr);
	map<string,vector<DATAPOINT> >::iterator chr_data = data.find(key);
	if(a == b) {
		return(0);
	}
	long int start 	= starts[key];
	long int end	= ends[key];
	if(a < b) {
		if(a < start) {
			return(0);
		}
		if(b > end) {
			return(0);
		}
		return calculate_rpkm_area(chr_data,a,b,start);
	}
	else {
		if(b < start) {
			return(0);
		}
		if(a > end) {
			return(0);
		}
		try {
			return calculate_rpkm_area(chr_data,b,a,start);
		}
		catch(...) {
			cerr << "Error occurred in calculate_rpkm: " << a << " " << " " << b << " " << start << endl;
		}
	}
}

float WigReader::rpkm(const string chr, long int a, long int b) {
	string key = locate_chr_key(chr);
	map<string,vector<DATAPOINT> >::iterator chr_data = data.find(key);
	if(a == b) {
		return(0);
	}
	long int start 	= starts[key];
	long int end	= ends[key];
	if(a < b) {
		if(a < start) {
			return(0);
		}
		if(b > end) {
			return(0);
		}
		return calculate_rpkm(chr_data,a,b,start);
	}
	else {
		if(b < start) {
			return(0);
		}
		if(a > end) {
			return(0);
		}
		try {
			return calculate_rpkm(chr_data,b,a,start);
		}
		catch(...) {
			cerr << "Error occurred in calculate_rpkm: " << a << " " << " " << b << " " << start << endl;
		}
	}
}

float WigReader::calculate_rpkm(map<string,vector<DATAPOINT> >::iterator chr_data, long int a, long int b, long int start, bool debug) {
	if(debug){
		cerr << "Debugging WigReader::calculate_rpkm()" << endl;
		cerr << "a:\t\t" << a << endl;
		cerr << "b:\t\t" << b << endl;
		cerr << "step:\t\t" << step << endl;		
		cerr << "start:\t\t" << start << endl;
	}
	
	float numerator   = 0.0;
	float denominator = 0.0;
    
    // The "pos" variable holds the current position along the chromosome.
    // We calculate the initial position by finding the first wiggle block containing "a"
	if(debug) {
		cerr << "Pos Calculation:" << endl;
		long double a_start = (long double)a - (long double)start;
		cerr << "(double)a - (double)start: " << a_start << endl;
		a_start = a_start/(long double)step;
		cerr << "((((double)a - (double)start))/(double)step): " << a_start << endl;
		long int pos = (long int)a_start;
		cerr << "(long int)((((double)a - (double)start))/(double)step): " << pos << endl;
		pos = pos * step + start;
		cerr << "(long int)((((double)a - (double)start))/(double)step) * step + start: " << pos << endl;
		
	}
	long int pos = (long int)((((double)a - (double)start))/(double)step); // (int) will act as floor.
	pos = pos * step + start;
	if(debug) {
		cerr << "Initial Pos:\t" << pos << endl;;
	}
    // If the first base of the wiggle block is not equal to the first base we are interested in
    // we must take a fraction of that wiggle block (equal to the # of bases overlapping with the region of interest)
    if(pos < a) {
		if(debug) {
			cerr << "pos < a" << endl;
		}
		int bases = a - pos + 1;
		if(bases < 0) {
			cerr << "pos<a: # of bases (" << bases << ") is <0!: " << step << ' ' << start << ' ' << pos;
			exit(1);
		}
		int offset = (pos - start)/step;
		float v = (chr_data->second)[offset].rpkm;
		numerator += v * bases / step;
		denominator += (float)bases / step;
		pos += step;
	}
	while(pos <= (b-step)) {
		if(debug) {
			cerr << "pos <= (b-step)" << endl;
			cerr << pos << " <= " << b-step << endl;
		}
		int offset = (pos - start)/step;
		numerator += (chr_data->second)[offset].rpkm;;
		denominator += (float)1;
		pos += step;
	}
	if(pos <= b) {
		if(debug) {
			cerr << "pos <= b" << endl;
			cerr << pos << " <= " << b << endl;
		}
		int bases = b - pos + 1;
		int offset = (pos - start)/step;
		float v = (chr_data->second)[offset].rpkm;
		numerator += v * bases / step;
		denominator += (float)bases / step;
	}
	float v = numerator / denominator;
	if(isnan(v) and !debug) {
		cerr << "Going into debug mode because v = " << v << endl;
		calculate_rpkm(chr_data, a, b, start, true);
		exit(1);
	}
	return(v);
}

float WigReader::calculate_rpkm_area(map<string,vector<DATAPOINT> >::iterator chr_data, long int a, long int b, long int start) {
    // The "pos" variable holds the current position along the chromosome.
    // We calculate the initial position by finding the first wiggle block containing "a"
	long int pos = (long int)((((double)a - (double)start))/(double)step); // (int) will act as floor.
	pos = pos * step + start;
	
	// 'area' is the accumulator variable for the area
	float area = 0.0;
	float v;
	
    // If the first base of the wiggle block is not equal to the first base we are interested in
    // we must take a fraction of that wiggle block (equal to the # of bases overlapping with the region of interest)
    if(pos < a) {
		int bases = a - pos + 1;
		if(bases < 0) {
			cerr << "pos<a: # of bases (" << bases << ") is <0!: " << step << ' ' << start << ' ' << pos;
			exit(1);
		}
		int offset = (pos - start)/step;
		v = (chr_data->second)[offset].rpkm;
		area += v * bases;
		pos += step;
	}
	while(pos <= (b-step)) {
		int offset = (pos - start)/step;
		area += (chr_data->second)[offset].rpkm;;
		pos += step;
	}
	if(pos <= b) {
		int bases = b - pos + 1;
		int offset = (pos - start)/step;
		float v = (chr_data->second)[offset].rpkm;
		area += v * bases;
	}
	return(area);
}

float WigReader::calculate_max_rpkm(map<string,vector<DATAPOINT> >::iterator chr_data, long int a, long int b, long int start) {
	float max = 0.0;
	long int pos = (long int)((((double)a - (double)start))/(double)step); // (int) will act as floor.
	pos = pos * step + start;

    // If the first base of the wiggle block is not equal to the first base we are interested in
	if(pos < a) {
		int offset = (pos - start)/step;
		float v = (chr_data->second)[offset].rpkm;
		if(v > max) {
			max = v;
		}
		pos += step;
	}
	while(pos <= (b-step)) {
		int offset = (pos - start)/step;
		float v = (chr_data->second)[offset].rpkm;
		if(v > max) {
			max = v;
		}
		pos += step;
	}
	if(pos <= b) {
		int bases = b - pos + 1;
		int offset = (pos - start)/step;
		float v = (chr_data->second)[offset].rpkm;
		if(v > max) {
			max = v;
		}
	}
	return(max);
}

std::vector<string> WigReader::chromosomes() {
	map<string, vector<DATAPOINT> >::const_iterator itr;
	std::vector<string> chroms;
	for (itr = data.begin(); itr != data.end(); ++itr) {
		chroms.push_back(itr->first);
	}
	return(chroms);
}

int	WigReader::startForChromosome(const string chr) {
	return(starts[chr]);
}

int	WigReader::endForChromosome(const string chr) {
	return(ends[chr]);
}

int WigReader::getStep() {
	return(step);
}