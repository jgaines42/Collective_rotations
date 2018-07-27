#include "E_recursive.h"  

using namespace std;

float  E_recursive(string save_name, vector<int> loc_array, float min_sofar, vector<float> prior_coordX, vector<float> prior_coordY,  vector<float> prior_coordZ, vector<float> prior_Size, int depth, int max_depth, vector<float> mins, vector<AminoAcid> all_residues, float current_E, ofstream& myfile){

AminoAcid this_res = all_residues[depth];
int num_toLoop = this_res.get_numDih();
float num_rest = 0;
if (depth < max_depth){ // if there will be more residues to go
	for (int i = depth+1; i <= max_depth; i ++){ // Loop over the rest of the residues
		num_rest = num_rest + mins[i];
	}
}


float this_E;
float new_E;

AminoAcid test_res = all_residues[1];

for (int loop_thisRes = 0; loop_thisRes < num_toLoop; loop_thisRes++){ //num_toLoop
	if (depth == 0)
		cerr << loop_thisRes << endl;
	new_E = 0;
	this_E = this_res.get_E(loop_thisRes);

	if (this_E + num_rest + current_E <= min_sofar){
		vector <float> new_priorX = prior_coordX;
		vector <float> new_priorY = prior_coordY;
		vector <float> new_priorZ = prior_coordZ;
		vector <float> new_priorSize = prior_Size;
		vector<float> new_X = this_res.get_x(loop_thisRes);

        vector<float> new_Y = this_res.get_y(loop_thisRes);
        vector<float> new_Z = this_res.get_z(loop_thisRes);

        vector<float> new_size = this_res.get_sizes();
		new_E = this_E + current_E;
			
        if (depth > 0){

            new_E = new_E + get_energy(prior_coordX, prior_coordY, prior_coordZ, new_X, new_Y, new_Z, prior_Size, new_size);	
        }

        if (new_E + num_rest <=min_sofar){

        	if (depth < max_depth){
        		for (int atoms = 0;atoms < this_res.get_numAtoms(); atoms++){

        			new_priorX.push_back(new_X[atoms]);
        			new_priorY.push_back(new_Y[atoms]);
        			new_priorZ.push_back(new_Z[atoms]);
					new_priorSize.push_back(new_size[atoms]);
				}

				loc_array[depth] = loop_thisRes;
				min_sofar = E_recursive(save_name, loc_array, min_sofar,new_priorX, new_priorY,  new_priorZ,  new_priorSize, depth+1,  max_depth,  mins, all_residues, new_E, myfile);
			} //if depth < max_depth
			else{ // depth == max_depth
				loc_array[depth] = loop_thisRes;

				if (new_E == min_sofar){

					if (loc_array[max_depth] > 0) {
						if (this_res.get_group(loc_array[max_depth]) > this_res.get_group(loc_array[max_depth]-1)) {
							//ofstream myfile(save_name.c_str(), ios::out | ios::app); 
							//myfile.open
							if (myfile.is_open()){
						
								for (int members = 0; members <= max_depth; members ++){
									myfile << loc_array[members] << " ";
								}
								myfile<<endl;
								//myfile << min_sofar << endl;
								//myfile.close();
							}
							else
								cerr << "Error opening output file " << endl;
						} //If new group
					}// If not first in list
					else { // First in group, so write
						
						if (myfile.is_open()){
						
							for (int members = 0; members <= max_depth; members ++){
								myfile << loc_array[members] << " ";
							}
							myfile<<endl;
							
						}
						else
							cerr << "Error opening output file " << endl;
					} //End else first in group
				}
				else {//new_E < min_sofar

					min_sofar = new_E;
					myfile.close();
					myfile.open(save_name.c_str(), ios::out); 

					if (myfile.is_open()){
						myfile << min_sofar << endl;

						for (int members = 0; members <= max_depth; members ++){
							myfile << loc_array[members] << " ";
						}
						myfile << endl;
					}
					else{
						cerr << "Error opening output file " << endl;
						}	
					}
			}
		}		
	}
}
 
   
return min_sofar;
} //End of E_recursive function

	
float get_energy(vector<float> prior_resX, vector<float> prior_resY, vector<float> prior_resZ,  vector<float> new_resX, vector<float> new_resY, vector<float> new_resZ, vector<float> As_prior, vector<float> As_new){
	int size_prior = prior_resX.size();
	int size_new = new_resX.size();
	float E = 0;
	float D1 = 0;
	float D2 = 0;

	for (int i = 0;i < size_prior;i++){
		for (int j = 0; j < size_new; j++){
			D1 = pow(prior_resX[i]-new_resX[j],2)+pow(prior_resY[i]-new_resY[j],2)+pow(prior_resZ[i]-new_resZ[j],2);
			D2 = sqrt(D1);
			if (D2 < (As_prior[i]+As_new[j])){
				 E = E+ pow((1 - pow((As_prior[i]+As_new[j])/D2,6)),2);

				 }
		
		}
	}
	return E;
} //End of get_energy function

	
	
	
