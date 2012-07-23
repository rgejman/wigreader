#include <iostream>
#include <map>
#include <vector>
#include <stdexcept>
class functions;

using std::runtime_error;


using namespace std;

struct DATAPOINT {
	long int start;
	float rpkm;
};

bool data_point_compare (DATAPOINT a, DATAPOINT b);
typedef std::map < std::string, vector<DATAPOINT> > DataMapType;
typedef std::map < std::string, long int > EndPositionsVector;

class OutOfWigRangeException : public runtime_error {
	public:
   		OutOfWigRangeException(): runtime_error( "position outside bounds of wig file" ) {}
};

class WigReader {
	public:
		WigReader(string file);
		WigReader(const char* file);
		void Read();
		float rpkm(const string chr, long int a, long int b);
		float max_rpkm(const string chr, long int a, long int b);
		float area(const string chr, long int a, long int b);
		std::vector<string> chromosomes();
		int	startForChromosome(const string chr);
		int	endForChromosome(const string chr);
		int getStep();
		DataMapType data;
		string locate_chr_key(const string chr);

	private:
		void init(const char* file);
		float calculate_rpkm(map<string,vector<DATAPOINT> >::iterator chr_data, long int a, long int b, long int start, bool debug=false);
		float calculate_rpkm_area(map<string,vector<DATAPOINT> >::iterator chr_data, long int a, long int b, long int start);
		float calculate_max_rpkm(map<string,vector<DATAPOINT> >::iterator chr_data, long int a, long int b, long int start);

		EndPositionsVector starts;
		EndPositionsVector ends;
		int step;
		const char* wigfile;
		string chr;
		long int last_pos;
		vector<DATAPOINT>* current_data_points_vector;
		bool skip_until_next_header;
		void ProcessLine(const string &line);
		void ProcessHeader(const string &header);
		void ProcessDataLine(const string &line);

};