
// a simple driver for the HFold

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

// include the simfold header files
#include "simfold.h"
#include "externs.h"
#include "h_globals.h"
//#include "h_externs.h"
#include "constants.h"
#include "params.h"

#include <iostream>


// Hosna June 20th, 2007
//#include "W_final.h"
#include "hfold.h"

//kevin June 23 2017
#include "hfold_validation.h"
#include <getopt.h>
#include <unistd.h>

// Ian Wark October 18 2017
#include "shape_data.h"

//kevin 27 oct

#include "s_specific_functions.h"
#include "Hotspot.h"
#include "h_common.h"



void printUsage();

int main (int argc, char *argv[])
{
    char sequence[MAXSLEN];
    char structure[MAXSLEN];
    char restricted[MAXSLEN];
    double energy;

/*
    if (argc != 3)
    {
        printf ("Usage: %s <sequence> <restricted_structure>\n", argv[0]);
        printf ("Example: %s \"GCAACGAUGACAUACAUCGCUAGUCGACGC\" \"(____________________________)\"\n", argv[0]);
        printf ("\tRestricted structure symbols:\n");
        printf ("\t\t() restricted base pair\n");
        printf ("\t\t_ no restriction\n");
        return 1;
    }
    */

    //kevin: june 22 2017
	//validation for command line argument
    char* inputPath;
	inputPath = (char*) malloc(sizeof(char) * 1000);

	char* outputPath;
	outputPath = (char*) malloc(sizeof(char) * 1000);

	bool sequenceFound = false;
	bool restrictedFound = false;
	bool inputPathFound = false;
	bool outputPathFound = false;
	bool errorFound = false;
	int option;
	int number_of_suboptimal_structure = 0;

	//kevin: june 23 2017 https://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html
        static struct option long_options[] = 
                {
                        {"s", required_argument, 0, 's'},
                        {"r", required_argument, 0, 'r'},
                        {"i", required_argument, 0, 'i'},
                        {"o", required_argument, 0, 'o'},
                        {"shape", required_argument, 0, 'g'},
                        {"b", required_argument, 0, 'h'},
                        {"m", required_argument, 0, 'm'},
						{"n", required_argument, 0, 'n'},
                        {0, 0, 0, 0}
                };


	while (1){
		// getopt_long stores the option index here.
                int option_index = 0;

		option = getopt_long (argc, argv, "s:r:i:o:", long_options, &option_index);

		// Detect the end of the options
		if (option == -1)
			break;

		switch (option)
		{
		case 's':
			if(sequenceFound){
				fprintf(stderr, "-s is duplicated\n");
				errorFound = true;
				break;
			}
			if(inputPathFound){
				printf("Cannot combine -i with -s/-r \n");
				errorFound = true;
				break;
			}
			strcpy(sequence, optarg);
			//printf("seq: %s\n",sequence);
			sequenceFound = true;
			break;
		case 'r':
			if(restrictedFound || inputPathFound){
				fprintf(stderr, "-r is duplicated\n");
				errorFound = true;
				break;
			}
			if(inputPathFound){
				fprintf(stderr, "Cannot combine -i with -s/-r \n");
				errorFound = true;
				break;
			}
			if(number_of_suboptimal_structure != 0){
				fprintf(stderr, "Cannot combine --r with --n \n");
				errorFound = true;
				break;
			}
			strcpy(restricted, optarg);
			restrictedFound = true;
			break;
		case 'i':
			if(restrictedFound || sequenceFound){
				fprintf(stderr, "Cannot combine -i with -s/-r \n");
				errorFound = true;
				break;
			}
			strcpy(inputPath,optarg);
			//printf("file: %s %d\n", file,access(file, F_OK|R_OK));
			if(access(inputPath, F_OK) == -1) { //if file does not exist
				fprintf(stderr, "Input file not exist\n");
				exit(4);
			}
			if (!validateHFOLDInputFile(inputPath, sequence, restricted, &sequenceFound, &restrictedFound)) {
				fprintf(stderr, "Input file is invalid\n");
				errorFound = true;
				break;
			}
			inputPathFound = true;
			break;
		case 'o':
			strcpy(outputPath, optarg);
			//printf("access: %d\n",access(output_path, F_OK));
			if(access(outputPath, F_OK) != -1) { //if file already exist
				addTimestamp(&outputPath);
			}
			outputPathFound = true;
			break;
		case 'g': //--shape (shape file path)
			if(!sequenceFound){
				fprintf(stderr, "Must define sequence before shape file\n");
				errorFound = true;
				break;
			}
				// important that this is before set_shape_file
				shape.set_sequence_length(strlen(sequence));
				shape.set_shape_file(std::string(optarg));
				break;
		case 'h': //--b (shape intercept)
				if (shape.is_number(optarg))
						shape.set_b(atof(optarg));
				else
						errorFound = true;
				break;
		case 'm': //--m (shape slope)
				if (shape.is_number(optarg))
						shape.set_m(atof(optarg));
				else
						errorFound = true;
				break;
		case 'n':
			number_of_suboptimal_structure = atoi(optarg);
			if(number_of_suboptimal_structure <= 0){
				fprintf(stderr, "number must be > 0\n");
				errorFound = true;
				break;
			}
			if(restrictedFound){
				fprintf(stderr, "Cannot combine --r with --n \n");
				errorFound = true;
				break;
			}
			break;
		default:
			errorFound = true;
			break;
		}
		//clean up when error
		if(errorFound){
			printUsage();
			exit(1);
		}
	}

	if(!(sequenceFound)){
		fprintf(stderr, "--s is missing\n");
		printUsage();
		free(inputPath);
    	free(outputPath);
		exit(1);
	}

	if(!validateSequence(sequence)){
		fprintf(stderr, "-s is invalid\n");
		//printUsage();
		free(inputPath);
    	free(outputPath);
		exit(1);
	}

	
	//kevin: june 22 2017 if we have output path and input path, try to combine both
	if(outputPathFound && inputPathFound){
		addPath(&outputPath, inputPath);
		//printf("out path: %s\n",outputPath);
	}
	//kevin: june 22 2017
	//end of validation for command line arguments

    //kevin; june 23 2017 took this out because variable is set during validation
    //strcpy (sequence, argv[1]);
    //strcpy (restricted, argv[2]);

    // Before calling any function in the library, you have to initialize config_file, dna_or_rna, temperature
    //     and to call the function init_data, which loads the thermodynamic parameters into memory

    // configuration file, the path should be relative to the location of this executable
    char config_file[200];
    strcpy (config_file, SIMFOLD_HOME "/params/multirnafold.conf");

    // what to fold: RNA or DNA
    int dna_or_rna;
    dna_or_rna = RNA;

    // temperature: any integer or real number between 0 and 100
    // represents degrees Celsius
    double temperature = 37.0;

    // initialize the thermodynamic parameters
    // call init_data only once for the same dna_or_rna and same temperature
    // if one of them changes, call init_data again
    init_data (argv[0], config_file, dna_or_rna, temperature);

	// Hosna, July 18, 2012
	// In simfold we have the following for RNA && temp=37
	fill_data_structures_with_new_parameters (SIMFOLD_HOME "/params/turner_parameters_fm363_constrdangles.txt");

	// Hosna, July 25, 2012
	// in HotKnots and ComputeEnergy package the most up-to-date parameters set is DP09.txt
	// so we add it here
	fill_data_structures_with_new_parameters (SIMFOLD_HOME "/params/parameters_DP09.txt");

	//kevin 27 oct
	std::vector<Hotspot*> hotspot_list;
	if(restrictedFound){
		if(!validateStructure(restricted, sequence)){
			fprintf(stderr, "-r is invalid\n");
			//printUsage();
			exit(1);
		}else{
			replaceBrackets(restricted);
		}
	}else{
		//printf("getting hotspot\n");
		get_hotspots(sequence, &hotspot_list);
		//printf("done hotspot\n");
		//exit(999);
	}


	Result* result;
	std::vector<Result*> result_list;

	if(restrictedFound){
		energy = hfold(sequence, restricted, structure);
		result = new Result(sequence,restricted,structure,energy);
		result_list.push_back(result);
	}else{

		//printf("number of hotspots: %d\n",hotspot_list.size());
		for (int i=0; i < hotspot_list.size(); i++){
			//printf("hotspot substructure #%d: %s\n",i,hotspot_list[i]->get_structure());
			energy = hfold(sequence,hotspot_list[i]->get_structure(),structure);
			result = new Result(sequence,hotspot_list[i]->get_structure(),structure,energy);
			//printf("%s\n%s\n%s\n%lf%d\n",result->get_sequence(),result->get_restricted(),result->get_final_structure(),result->get_energy(),result->get_method_chosen());
			result_list.push_back(result);
		}

		std::sort(result_list.begin(), result_list.end(),compare_result_ptr);
	}

/*
    // check if restricted is included in structure
	// Hosna March 7, 2012
	// for optimality we get the value once, use it many times
	int seqLen = strlen (sequence);
    for (int i=0; i < seqLen; i++)
    {
        if ((restricted[i] == '(' || restricted[i] == ')' || restricted[i] == '.') &&
            (restricted[i] != structure[i]))
        {
            fprintf(stderr, "There is something wrong with the structure, doesn't match restricted\n");
			fprintf(stderr, "  %s\n  %s\n  %s\t%.2lf\n", sequence, restricted, structure, energy);
			exit(1);
        }
    }
*/

/*
    //kevin 22 June 2017
	//different ways of outputing
	if(outputPathFound){
		FILE* fp;
		fp = fopen(outputPath,"w");
		if(fp){
			fprintf(fp,"Sequence: %s\n",sequence);
			fprintf(fp,"Input_structure: %s\n",restricted);
			fprintf(fp,"Output_structure: %s\n",structure);
			fprintf(fp,"Energy: %.2lf\n",energy);
			fclose(fp);
		}
	}else{
		printf ("Seq: %s\n", sequence);
        printf ("RES: %s  %.2lf\n", structure, energy);
	}
*/

	//kevin 5 oct 2017
	int number_of_output;
	//printf("number_of_suboptimal_structure: %d\n",number_of_suboptimal_structure);
	if(number_of_suboptimal_structure != 0){
		number_of_output = MIN(result_list.size(),number_of_suboptimal_structure);
	}else{
		number_of_output = 1;
	}

	//kevin: june 22 2017
	//output to file
	if(outputPathFound){
		bool write_success = write_output_file(outputPath, number_of_output, result_list);
		if(!write_success){
			fprintf(stderr, "write to file fail\n");
			exit(4);
		}
	}else{
		//kevin 5 oct 2017
		printf("Seq: %s\n",sequence);
		for (int i=0; i < number_of_output; i++) {
			printf("Restricted_%d: %s\n",i, result_list[i]->get_restricted());
			printf("Result_%d: %s\nEnergy_%d: %lf\n",i, result_list[i]->get_final_structure(),i,result_list[i]->get_final_energy());
		}
	}

    free(inputPath);
    free(outputPath);
	for (int i=0; i < result_list.size(); i++) {
		delete result_list[i];
	}

	for(int i =0; i<hotspot_list.size(); i++){
		delete hotspot_list[i];
	}
    return 0;
}

void printUsage(){
	/*
	printf ("\nUsage: HFold_iterative <sequence> <structure>\n");
	printf ("Example: ./HFold_iterative \"GCAACGAUGACAUACAUCGCUAGUCGACGC\" \"(____________________________)\" \n");
	printf ("Or \nUsage: HFold_iterative <path to input file> <path to output file>\n");
	printf ("Example: ./HFold_iterative \"/home/username/Desktop/inputFile.txt\" -o \"/home/username/Desktop/outFile.txt\" \n");
	printf ("\tRestricted structure symbols:\n");
	printf ("\t\t() restricted base pair\n");
	printf ("\t\t _ no restriction\n");
*/
	printf("Usage ./HFold --s <sequence> --r <structure> [--o </path/to/file>]\n");
	printf("or\n");
	printf("Usage ./HFold --i </path/to/file> [--o </path/to/file>]\n");
	printf ("  Restricted structure symbols:\n");
	printf ("    () restricted base pair\n");
	printf ("    _ no restriction\n");
	printf("Example:\n");
	printf("./HFold --s \"GCAACGAUGACAUACAUCGCUAGUCGACGC\" --r \"(____________________________)\"\n");
	printf("./HFold --i \"/home/username/Desktop/myinputfile.txt\" --o \"/home/username/Desktop/some_folder/outputfile.txt\"\n\n");

	printf("You can also include SHAPE data to be used.\n");
        printf("The SHAPE data must be in a file with 1 number per line.\n");
	printf("The number corresponds with each nucleotide in order, and the file must be exactly the same length as the sequence.\n");
        printf("--shape (\"filename\") to specify a file for shape data\n");
        printf("--b (number) to specify an intercept for the shape data (default is -0.600000)\n");
        printf("--m (number) to specify a slope for the shape data (default is 1.800000)\n\n");

	printf("Example:\n");
        printf("./HFold --s \"GCAACGAUGACAUACAUCGCUAGUCGACGC\" --r \"(____________________________)\" --shape \"shapefile\" --b -0.4 --m 1.3\n\n");

	printf("Please read README for more details\n");
}

//kevin 27 oct
bool write_output_file(char* path_to_file, int num_of_output, std::vector<Result*> result_list){
	FILE* fp = fopen(path_to_file,"w");
	if (fp == NULL) {
		return false;
	}
	fprintf(fp,"Seq: %s\n",result_list[0]->get_sequence());
	for (int i=0; i < num_of_output; i++) {
		fprintf(fp,"Restricted_%d: %s\n",i, result_list[i]->get_restricted());
		fprintf(fp,"Result_%d: %s\nEnergy_%d: %lf\n",i, result_list[i]->get_final_structure(),i,result_list[i]->get_final_energy());
	}
	fclose(fp);
	return true;
}