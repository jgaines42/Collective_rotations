%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [] = process_C_recursive(pdb_names, variant, group, folder_data, folder_pdb)
% 
% Processes the output of the recursive C++ code. Takes all the output
% files and condenses the output into matlab files (removes cases where the
% only difference is CH3 rotations)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = process_C_recursive(pdb_names, variant, group, folder_data, folder_pdb, database)
folder_name = folder_data;
folder1 = strcat(folder_data, pdb_names, '/');

load(strcat(folder_name, pdb_names, '_all_core_binned_2.mat'));
ind0 = new_moving(:,2) == group;
all_moving = double(new_moving(ind0,1));

load(strcat(folder_pdb, pdb_names, '.mat'));

if size(all_moving,1) >1 && size(all_moving,1) < 11
     res_ids = cell2mat(tempModel2(:,6));
    for i = 2:size(tempModel2,1)
        if res_ids(i,1) < res_ids(i-1,1)
            res_ids(i,1) = res_ids(i,1) + 1000;
        end
    end
    tempModel2(:,6) = num2cell(res_ids);
    Residue_id = [];
    Residue_names = {};
    all_DOF = zeros(size(all_moving,1),1);
    for i = 1:size(all_moving,1)
        ind0 = find(res_ids == all_moving(i));
        resiName = tempModel2{ind0(1),4}; %Get name of interface resiude
        resiName(2:3) = lower(resiName(2:3));
        if database == 'h'
            Residue_id = [Residue_id,(all_moving(i))];
            Residue_names = [Residue_names, resiName];
        elseif ismember(resiName, { 'Ala','Ile', 'Leu', 'Phe', 'Met', 'Val'})
            Residue_id = [Residue_id,(all_moving(i))];
            Residue_names = [Residue_names, resiName];
        end
        if  resiName == 'Leu' | resiName == 'Ile' | resiName == 'Phe' | resiName == 'His' | resiName == 'Tyr' | resiName == 'Trp'
            all_DOF(i) = 2;
        elseif resiName == 'Met'
            all_DOF(i) = 3;
        else
            all_DOF(i) = 1;
        end
    end
    
    
    f = fopen(strcat(folder1, pdb_names,  '_loc_out_', num2str(variant), '_g', num2str(group), '.txt'));
    counter = 1;
    all_Energy = [];
    unique_dihedrals = [];
    all_Energy = textscan(f, '%f',1);
    Energy = all_Energy{:};
    old_size  = 1000;
    while ~feof(f)
        data = textscan(f,strcat(repmat('%d ',1,size(all_moving,1))),1000);
       % Energy = data{size(all_moving,1)+1};
        all_dih = zeros(size(data{1},1), sum(all_DOF));
        for i = 1:size(all_moving,1)
            load(strcat(folder1, pdb_names, '_', Residue_names{i},num2str( Residue_id(i)), '_dihedrals_fast2_', num2str(variant), '.mat'))
            this_res_ind = data{i}+1; %Cause of zero indexing in C++
            this_res_dih = dihedrals(this_res_ind,:);
            all_dih(:,sum(all_DOF(1:i-1))+1:sum(all_DOF(1:i))) =[this_res_dih(:,1:all_DOF(i))];
        end
        unique_dihedrals = [unique_dihedrals;unique(all_dih,'rows')];
        if mod(counter,200)==0 && size(unique_dihedrals,1)-old_size > 50000;
            unique_dihedrals = unique(unique_dihedrals,'rows');
            old_size = size(unique_dihedrals,1);
          % fprintf('%f %f \n', counter, old_size)
        end
        counter = counter + 1;
        
      %  all_Energy = unique([all_Energy;Energy]);
    end  
    fclose(f);
    Energy = all_Energy;
   

    unique_dihedrals = unique(unique_dihedrals, 'rows');

    save(strcat(folder1, pdb_names, '_unique_dihedrals_fast2_', num2str(variant), '_g', num2str(group), '.mat'), 'unique_dihedrals');
    save(strcat(folder1, pdb_names, '_Energy_fast2_', num2str(variant), '_g', num2str(group), '.mat'), 'Energy');
    if variant > 1
        e = strcat('rm ', {' '}, folder1, pdb_names,  '_loc_out_', num2str(variant), '_g', num2str(group), '.txt');
        system(e{:});
    end
end
end