
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [] = find_interactions(pdb,res_to_run,folder_name, save_folder)
%
% Finds interaction networks of a set of residues in a protein.
%
% Input:
%   pdb: 4 letter PDB code.
%   res_to_run: array of residue ids to cluster
%   folder_name: folder where PDB.mat files are stores
%   save_folder: folder to save results to
%
% Output:
%   PDB_all_core_paired_data.mat which contains 'paired', an nxn
%       array with 1s for each interaction pair. All amino acids in
%       res_to_run are included
%   PDB_all_core_binned.mat
%       Column 1: Amino acid ids
%       Column 2: Which group the amino acid is in. All pairs that are not
%       ILMTSWYV are removed and put in individual groups. To change this, 
%       comment out as indicated on line 74
%
% Notes:
% - Interactions are defined as 2 residues being able to overlap without
%  first overlapping with the backbone
% - Set up PDB file using download_preprocess_pdb.py
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = find_interactions(pdb,res_to_run,folder_name, save_folder)

res_to_run = res_to_run;
PDB_name = strcat(folder_name, pdb, '.mat');

load(PDB_name);
tempModel2(:,1) = num2cell([1:size(tempModel2,1)]);
tempModel2 = add_sizes_protein(tempModel2,9);

ind1 = strcmp(tempModel2(:,3),'B');
x = tempModel2(~ind1,:);
tempModel2 = tempModel2(~ind1,:);
tempModel2(:,1) = num2cell([1:size(tempModel2,1)]);
ind1 = strcmp(tempModel2(:,5), 'A');
tempModel2 = tempModel2(ind1,:);

%Create NxN array of interactions
paired = zeros(size(res_to_run,1), size(res_to_run,1));

%First check for clashes in initial positions
for i = 1:size(res_to_run,1)-1
    for j = i+1:size(res_to_run,1)
        paired(i,j) = residue_interaction(tempModel2,i,j, 0, res_to_run);
    end
end

%Then do rotations of residues to see if they can interact
for i = 1:size(res_to_run,1)-1
    for j = i+1:size(res_to_run,1)
        if paired(i,j) == 0 %If an interaction was already made, don't do it again
            paired(i,j) = residue_interaction(tempModel2,i,j,1,res_to_run);
        else
            
        end
    end
end

save(strcat(save_folder, pdb, '_all_paired_data.mat'), 'paired')

% Also make list of the groups
ind0 = ismember(cell2mat(tempModel2(:,6)), res_to_run);
ind1 = ismember(tempModel2(:,2), {'N'});
these_res = tempModel2(ind0&ind1,:);

%Remove ala,gly and cys from the pairing
%% Remove the following lines to keep all amino acids to the line marked
% 'end remove'
ind0 = ismember(these_res(:,4),{'ALA', 'GLY', 'CYS'});
ind1 = find(ind0 == 1);
paired(ind1,:) = 0;
paired(:,ind1) = 0;

for i = 1:size(paired,1)
    for j = i+1:size(paired,1)
        if paired(i,j) == 1
            paired(j,i) =1;
        end
    end
end

cdata = zeros(size(paired,1),1);
ind0 = ismember(these_res(:,4), {'ILE','LEU','PHE','VAL','MET', 'THR', 'SER', 'TYR', 'TRP'});
cdata(ind0) = 15;
cdata(~ind0) = 5;
%% End remove

g = graph(paired);

bins = conncomp(g);
new_moving = [cell2mat(these_res(:,6)),bins'];
save(strcat(save_folder, pdb, '_all_core_binned.mat'), 'new_moving');

end