// Iterative HFold files
#include "hotspot.hh"
#include "Result.hh"
#include "cmdline.hh"
#include "hfold_validation.h"
#include "HFold.hh"
#include "W_final.h"
//simfold files
#include "s_specific_functions.h"
#include "simfold.h"
#include "externs.h"
#include "h_globals.h"
#include "constants.h"
#include "params.h"
#include "common.h"
// a simple driver for the HFold
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stack>
#include <sys/stat.h>
#include <string>
#include <getopt.h>


// this causes less compilation warnings than the #defines
char HFOLD[8] =						"./HFold";
char SIMFOLD[10] =                  "./simfold";

bool exists (const std::string path) {
  struct stat buffer;   
  return (stat (path.c_str(), &buffer) == 0); 
}

void get_input(std::string file, std::string &sequence, std::string &structure){
	if(!exists(file)){
		std::cout << "Input file does not exist" << std::endl;
		exit(EXIT_FAILURE);
	}
	std::ifstream in(file.c_str());
	std::string str;
	int i = 0;
	while(getline(in,str)){
		if(str[0] == '>') continue;
		if(i==0) sequence = str;
		if(i==1) structure = str;
		++i;
	}
	in.close();
}

//check length and if any characters other than ._()
void validateStructure(std::string sequence, std::string structure){
	if(structure.length() != sequence.length()){
		std::cout << " The length of the sequence and corresponding structure must have the same length" << std::endl;
		exit(EXIT_FAILURE);
	}

	//check if any characters are not ._()
	for(char c : structure) {
		if (!(c == '.' || c == '_' || c == '(' || c == ')')){
			std::cout << "Structure must only contain ._(): " << c << std::endl;
			exit(EXIT_FAILURE);
		}
	}
}

//check if sequence is valid with regular expression
//check length and if any characters other than GCAUT
void validateSequence(std::string sequence){

	if(sequence.length() == 0){
		std::cout << "sequence1 or sequence2 is missing" << std::endl;
		exit(EXIT_FAILURE);
	}
  // return false if any characters other than GCAUT -- future implement check based on type
  for(char c : sequence) {
    if (!(c == 'G' || c == 'C' || c == 'A' || c == 'U' || c == 'T')) {
		std::cout << "Sequence contains character " << c << " that is not G,C,A,U, or T." << std::endl;
		exit(EXIT_FAILURE);
    }
  }
}

double hfold(char *sequence, char *restricted, char *structure){
	W_final *min_fold = new W_final (sequence, restricted);
	if (min_fold == NULL) giveup ("Cannot allocate memory", "HFold");
	double energy = min_fold->hfold();
    min_fold->return_structure (structure);
    delete min_fold;
    return energy;
}


int main (int argc, char *argv[])
{
    args_info args_info;

	// get options (call getopt command line parser)
	if (cmdline_parser (argc, argv, &args_info) != 0) {
	exit(1);
	}

	std::string seq;
	if (args_info.inputs_num>0) {
	seq=args_info.inputs[0];
	} else {
		if(!args_info.input_file_given) std::getline(std::cin,seq);
	}
	int n = seq.length();

	validateSequence(seq);

	std::string restricted;
    args_info.input_structure_given ? restricted = input_struct : restricted = "";

	if(restricted != "") validateStructure(seq,restricted);

	std::string fileI;
    args_info.input_file_given ? fileI = input_file : fileI = "";

	std::string fileO;
    args_info.output_file_given ? fileO = output_file : fileO = "";	

	int number_of_suboptimal_structure = args_info.subopt_given ? subopt : 1;

	bool pk_free = args_info.pk_free_given;

	if(fileI != ""){
		
		if(exists(fileI)){
			get_input(fileI,seq,restricted);
		}
		if(seq == ""){
			std::cout << "sequence is missing from file" << std::endl; 
		}
		
	}

    // configuration file, the path should be relative to the location of this executable
   	char config_file[400];
	strcpy (config_file, SIMFOLD_HOME "/params/multirnafold.conf");

   //what to fold: RNA or DNA
	int dna_or_rna= RNA;
	// represents degrees Celsius
	double temperature = 37.0;
	// initialize the thermodynamic parameters
	init_data ("./simfold", config_file, dna_or_rna, temperature);

	fill_data_structures_with_new_parameters ( SIMFOLD_HOME "/params/turner_parameters_fm363_constrdangles.txt");

	// in HotKnots and ComputeEnergy package the most up-to-date parameters set is DP09.txt
	// so we add it here
	fill_data_structures_with_new_parameters ( SIMFOLD_HOME "/params/parameters_DP09.txt");

	std::vector<Hotspot> hotspot_list;

	char sequence[n+1];
	strcpy(sequence,seq.c_str());
	// Hotspots

	
	if(restricted != ""){
		Hotspot hotspot(0,restricted.length()-1,restricted.length());
		hotspot.set_structure(restricted);
		hotspot_list.push_back(hotspot);
	}
	else {
		get_hotspots(sequence, hotspot_list,number_of_suboptimal_structure);
	}

	// Data structure for holding the output
	std::vector<Result> result_list;

    //double min_energy;
	char final_structure[MAXSLEN];
	// Iterate through all hotspots or the single given input structure
	for(int i = 0;i<hotspot_list.size();++i){
		char structure[n+1];
		double energy;
		std::string struc = hotspot_list[i].get_structure();
		strcpy(structure,struc.c_str());

		char final_structure[MAXSLEN];
		if(!pk_free){
			energy = hfold(sequence, structure, final_structure);
		}
		else{
			call_simfold(SIMFOLD, sequence, structure, final_structure, &energy);
		}

		std::string final(final_structure);
		Result result(seq,hotspot_list[i].get_structure(),hotspot_list[i].get_energy(),final,energy);
		result_list.push_back(result);
	}

    

Result::Result_comp result_comp;
	std::sort(result_list.begin(), result_list.end(),result_comp );

	int number_of_output = 1;

	if(number_of_suboptimal_structure != 1){
			number_of_output = std::min( (int) result_list.size(),number_of_suboptimal_structure);
	}

	//output to file
	if(fileO != ""){
		std::ofstream out(fileO);
		out << sequence << std::endl;
		for (int i=0; i < number_of_output; i++) {
			out << "Restricted_" << i << ": " << result_list[i].get_restricted() << std::endl;;
			out << "Result_" << i << ":     " << result_list[i].get_final_structure() << " (" << result_list[i].get_final_energy() << ")" << std::endl;	
		}

	}else{
		//kevin: june 22 2017
		//Mateo: Sept 13 2023
		//changed format for ouptut to stdout
		std::cout << seq << std::endl;
		if(result_list.size() == 1){
			// std::cout << "Restricted_" << 0 << ": " << result_list[0].get_restricted() << std::endl;;
			std::cout << result_list[0].get_final_structure() << " (" << result_list[0].get_final_energy() << ")" << std::endl;
		}
		else{
			for (int i=0; i < number_of_output; i++) {
				if(result_list[i].get_final_structure() == result_list[i-1].get_final_structure()) continue;
				std::cout << "Restricted_" << i << ": " << result_list[i].get_restricted() << std::endl;;
				std::cout << "Result_" << i << ":     " << result_list[i].get_final_structure() << " (" << result_list[i].get_final_energy() << ")" << std::endl;
			}
		}
	}
	cmdline_parser_free(&args_info);

    return 0;
}

bool call_simfold (char *programPath, char *input_sequence, char *input_structure, char *output_structure, double *output_energy) {
        std::string result = "";

		char config_file[400];
		strcpy (config_file, SIMFOLD_HOME "/params/multirnafold.conf");

		double temperature;
		temperature = 37;
		init_data ("./simfold", config_file, RNA, temperature);

        fill_data_structures_with_new_parameters (SIMFOLD_HOME "/params/turner_parameters_fm363_constrdangles.txt");
	// when I fill the structures with DP09 parameters, I get a segmentation fault for 108 base sequence!!!!
	// So I chopped the parameter set to only hold the exact number as the turner_parameters_fm363_constrdangles.txt,
	// but still getting seg fault!
	fill_data_structures_with_new_parameters (SIMFOLD_HOME "/params/parameters_DP09_chopped.txt");

	*output_energy = simfold_restricted (input_sequence, input_structure, output_structure);
//	printf ("Call_Simfold_RES( can be called by different methods): %s  %.2lf\n", output_structure, output_energy);
	if(is_invalid_restriction(input_structure,output_structure)){
			fprintf(stderr,"ERROR!!! There is something wrong with the structure, doesn't match restricted\n");
			fprintf(stderr,"  %s\n  %s\n  %s\t%.2lf\n", input_sequence, input_structure, output_structure, *output_energy);
			fprintf(stderr,"ERROR!!! There is something wrong with the structure, doesn't match restricted\n");
			exit(11);
	}
	return true;
}

//---------------------------------------this function is suppose to be the same as the one in Hfold_interacting, if any changes are made, please change that one too--------------------
//kevin 30 Aug 2017
//check if the computed structure matches the restricted structure
int is_invalid_restriction(char* restricted_structure, char* current_structure){
	std::string openBracketArray ("({[");
	std::string closeBracketArray (")}]");

	for (int i=0; i < strlen(restricted_structure); i++){
        if(restricted_structure[i] != '_' && restricted_structure[i] != current_structure[i]){
			if( (openBracketArray.find_first_of(restricted_structure[i]) != -1) && ((openBracketArray.find_first_of(current_structure[i]) != -1)) ){
				continue;
			}else if ( (closeBracketArray.find_first_of(restricted_structure[i]) != -1) && ((closeBracketArray.find_first_of(current_structure[i]) != -1)) ){
				continue;
			}else{
				return 1;
			}
		}

    }
	return 0;
}
