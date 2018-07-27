%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function []  = set_up_recursive_run( pdb,folder_data, folder_code,folder_pdb)
%
% pdb: 4 letter code for PDB
% folder_data: folder containing *all_core_binned_2.mat file
% folder_code: folder containing recursive code
% folder_pdb: folder containing XXXX.mat file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function []  = set_up_recursive_run( pdb,folder_data, folder_code,folder_pdb)


load(strcat(folder_data, pdb, '_all_core_binned_2.mat'));
groups = unique(new_moving(:,2));
f = fopen(strcat('run_', pdb, '_combined.sh'), 'w');

for i = 1:size(groups,1)
    this_group = new_moving(new_moving(:,2)==groups(i),:);
    if size(this_group,1) > 1
        for j = 1:50
             fprintf(f, 'cd %s ; bash run_HQ54_bash.sh %s %d %d %s %s ;\n', folder_code,pdb, j, groups(i), folder_data, folder_pdb);
        end
    end
end
fclose(f);

if ~(exist(strcat(folder_data, pdb, '/'))==7)
    a = strcat('mkdir', {' '}, folder_data, pdb, '/');
    system(a{:});
end
end