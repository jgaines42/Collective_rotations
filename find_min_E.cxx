#include "AminoAcid.h"
#include "E_recursive.h" 
#include <iostream>
#include <fstream>
#include <ostream>
#include <cstdlib>
#include <math.h>
#include <string.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>  

using namespace std;



int main(int argc, char **argv){

/*******************************
 * Process initial arguments
 *******************************/
	// Arguments:
	// pdb_name: 4 characters
	// variant: int 
	// database: 1 character
	// group: int
	if (argc < 5){
			cerr << " Invalid input " << endl ;
			cerr << endl;
			return 0;
		
	}
	
	// Get pdb file abbrev.
	string file_name;
	file_name= argv[1];
	
	// Get variant
	int variant = atoi(argv[2]);
	int size_char = 0;
	if (variant < 10)
		size_char = 2;
	else if (variant < 100)
		size_char = 3;
	else
		size_char = 4;
	char *charV = argv[2];
	char char_var[size_char];
	strcpy(char_var, charV);
	
	string folder_data;
	folder_data = argv[4];
	string folder_pdb;
	folder_pdb = argv[5];
	
	string group;
	group = argv[3];
	int return_value = 0;
/*******************************
 * Load parameter file
 *******************************/	
	string folder_name;
	string param_file;
	param_file = folder_data + file_name + "/" + file_name + "_parameters_" + char_var +"_g" +  group + ".txt";
	folder_name = folder_data + file_name + "/";
	
	string line;
	ifstream myfile;
	cerr << param_file << endl;

	myfile.open(param_file.c_str());

	//ifstream myfile (param_file);
	if (myfile.is_open())
	  {
		// Parameter data
		char name[4];
		int max_depth;
		float min_E;
		myfile >> name >> max_depth >> min_E;
		cerr << max_depth << endl;
		string *all_res1  = new string[max_depth+1];
		int all_id[max_depth+1];

		vector<float> all_min;
		float temp_min = 0;
		int num_atoms[max_depth+1];
		int coord_width[max_depth+1];
		int counter = 0;
		while (!myfile.eof() )
		{
		  myfile >> all_res1[counter] >> all_id[counter] >> temp_min >> num_atoms[counter] >> coord_width[counter] ;
		  all_min.push_back(temp_min);
		  cerr << all_res1[counter] << endl;
		  counter++;
		}
		myfile.close();
 		/*******************************
		 * Set up Residues
		 *******************************/	
 		vector<AminoAcid>  all_residues;
 		 for (int i = 0; i < max_depth; i ++){
 		 	AminoAcid x  =  AminoAcid(variant, file_name, all_id[i], all_res1[i], coord_width[i]/3, num_atoms[i], folder_name);
 		 	all_residues.push_back(x);
 		 	}
 		 /*****************************
 		 * Set up recursion
 		 ******************************/	
 		 string save_name;
		save_name = folder_data +  file_name+ "_loc_out_" + char_var + "_g" + group + ".txt";
		vector<int> loc_array;
		for (int i = 0; i < max_depth; i++){
			loc_array.push_back(0);
			}
		vector<float> prior_coordX;
		vector<float> prior_coordY;
		vector<float> prior_coordZ;
		vector<float> prior_Size;

			max_depth = max_depth -1; //Change to 0 indexing

		//Clear output file
		ofstream myfile1(save_name.c_str(), ios::out);

		if (myfile1.is_open()){
			myfile1 << min_E << endl;
			//myfile1.close();
		}

		cerr << "got to recursion" << endl;
		float E = E_recursive(save_name, loc_array, min_E, prior_coordX, prior_coordY,  prior_coordZ, prior_Size, 0, max_depth, all_min, all_residues,0, myfile1);
		cerr << "finished recursion" << endl;

		if (E < min_E)
		{
			return_value = 1;
			}
		

		 /*******************************
		 * Delete pointer arrays
		 *******************************/	
 		AminoAcid this_res;
 		for (int i = 0; i<=max_depth; i++){
 		 	cerr << "i " << i << endl;

 			this_res = all_residues[i];
 			cerr << "i " << i << endl;
 			this_res.delete_file(variant, file_name, all_id[i], all_res1[i], coord_width[i]/3, num_atoms[i], folder_name);
 			}
 		
 		delete [] all_res1;
 		
 		
 	}
	else 
	{
		cerr << "Unable to open parameter file" << endl; // Parameter file didn't exist
		return_value = 1;
 	}

 cerr << return_value << endl;
 cout << return_value << endl;
	

}