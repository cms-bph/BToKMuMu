// -----------------------------------------------
//       Author: Xin Shi <Xin.Shi@cern.ch> 
//       Created:   [2013-08-20 Tue 13:54] 
// -----------------------------------------------
#include <TChain.h>

void set_root_style(int stat=1110, int grid=0); 
TChain* add_chain(TString datatype, TString label, TString cut, int verbose=1); 
char* get_option(char ** begin, char ** end, const std::string & option); 
bool option_exists(char** begin, char** end, const std::string& option); 
int get_number_of_lines(TString, TString, TString &); 
double calc_scale_factor(TString datatype, TString energy); 

