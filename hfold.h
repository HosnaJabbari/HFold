#ifndef HFOLD_H_
#define HFOLD_H_

#include <vector>
#include "Result.h"

double hfold(char *sequence, char *restricted, char *structure);
double hfold_pkonly(char *sequence, char *restricted, char *structure); // April 3, 2012
bool write_output_file(char* path_to_file, int num_of_output, std::vector<Result*> result_list);
#endif /*HFOLD_H_*/
