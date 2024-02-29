#ifndef HFOLD_H_
#define HFOLD_H_


bool call_simfold (char *programPath, char *input_sequence, char *input_structure, char *output_structure, double *output_energy);
int is_invalid_restriction(char* restricted_structure, char* current_structure);

#endif /*HFOLD_H_*/
