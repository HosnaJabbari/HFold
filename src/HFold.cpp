// Iterative HFold files
#include "hotspot.hh"
#include "Result.hh"
#include "cmdline.hh"
#include "W_final.hh"
#include "h_globals.hh"
// a simple driver for the HFold
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdio.h>
#include <sys/stat.h>
#include <string>
#include <getopt.h>

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
		if (!(c == '.' || c == '_' || c == '(' || c == ')' || c == 'x')){
			std::cout << "Structure must only contain ._()x: " << c << std::endl;
			exit(EXIT_FAILURE);
		}
	}
}

//check if sequence is valid with regular expression
//check length and if any characters other than GCAUT
void validateSequence(std::string sequence){

	if(sequence.length() == 0){
		std::cout << "sequence is missing" << std::endl;
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

std::string hfold(std::string seq,std::string res, double &energy, sparse_tree &tree, bool pk_free, bool pk_only, int dangle){
	W_final min_fold(seq,res, pk_free, pk_only, dangle);
	energy = min_fold.hfold(tree);
    std::string structure = min_fold.structure;
    return structure;
}

void seqtoRNA(std::string &sequence){
    for (char &c : sequence) {
      	if (c == 'T') c = 'U';
    }
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

	std::string restricted;
    args_info.input_structure_given ? restricted = input_struct : restricted = "";

	std::string fileI;
    args_info.input_file_given ? fileI = input_file : fileI = "";

	std::string fileO;
    args_info.output_file_given ? fileO = output_file : fileO = "";	

	int number_of_suboptimal_structure = args_info.subopt_given ? subopt : 1;

	bool pk_free = args_info.pk_free_given;

	bool pk_only = args_info.pk_only_given;

	int dangles = args_info.dangles_given ? dangle_model : 2;

	if(fileI != ""){
		
		if(exists(fileI)){
			get_input(fileI,seq,restricted);
		}
		if(seq == ""){
			std::cout << "sequence is missing from file" << std::endl; 
		}
		
	}
	int n = seq.length();
	std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
	if(!args_info.noConv_given) seqtoRNA(seq);

	validateSequence(seq);
	if(restricted != "") validateStructure(seq,restricted);

	std::string file= args_info.paramFile_given ? parameter_file : "params/rna_Turner2004.par";
	if(exists(file)){
		vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
	}
	else if (seq.find('T') != std::string::npos){
		vrna_params_load_DNA_Mathews2004();
	}

	std::vector<Hotspot> hotspot_list;

	// Hotspots

	vrna_param_s *params;
	params = scale_parameters();
	if(restricted != ""){
		Hotspot hotspot(1,restricted.length(),restricted.length()+1);
		hotspot.set_structure(restricted);
		hotspot_list.push_back(hotspot);
	}
	if((number_of_suboptimal_structure-hotspot_list.size())>0) {
		get_hotspots(seq, hotspot_list,number_of_suboptimal_structure,params);
	}
	free(params);
	// Data structure for holding the output
	std::vector<Result> result_list;

    //double min_energy;
	// Iterate through all hotspots or the single given input structure
	cand_pos_t size = hotspot_list.size();
	for(int i = 0;i<size;++i){
		double energy;
		std::string structure = hotspot_list[i].get_structure();

		sparse_tree tree(structure,n);
		std::string final_structure = hfold(seq,structure, energy,tree,pk_free,pk_only, dangles);
		
		Result result(seq,hotspot_list[i].get_structure(),hotspot_list[i].get_energy(),final_structure,energy);
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
		out << seq << std::endl;
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
			std::cout << "Restricted_" << 0 << ": " << result_list[0].get_restricted() << std::endl;;
			std::cout << "Result_" << 0 << ":     " << result_list[0].get_final_structure() << " (" << result_list[0].get_final_energy() << ")" << std::endl;
			for (int i=1; i < number_of_output; i++) {
				if(result_list[i].get_final_structure() == result_list[i-1].get_final_structure()) continue;
				std::cout << "Restricted_" << i << ": " << result_list[i].get_restricted() << std::endl;;
				std::cout << "Result_" << i << ":     " << result_list[i].get_final_structure() << " (" << result_list[i].get_final_energy() << ")" << std::endl;
			}
		}
	}
	cmdline_parser_free(&args_info);

    return 0;
}
