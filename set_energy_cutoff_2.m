%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function min_energy = set_energy_cutoff_2(coor_m1,coor_m2,coor_m3,coor_m4,coor_m5,coor_m6,coor_m7, coor_m8, coor_m9,coor_m10, num_res)
%
% Calculates the global energy due to lowest energy state of each residue
%
% Input:
%   coord_*: the coordinates of the lowest energy state for each residue.
%       The energy state value is stored in the last row
%   num_res: The number of residues
%
% Output:
%   min_energy: The energy due to atomic overlaps of the residue in this
%       configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function min_energy = set_energy_cutoff_2(coor_m1,coor_m2,coor_m3,coor_m4,coor_m5,coor_m6,coor_m7, coor_m8, coor_m9,coor_m10, num_res)
min_energy = 0;

%Calculate the global energy by doing pairwise interactions of each residue
for loop_cores = 1:num_res
    current_coord = eval(strcat('coor_m', num2str(loop_cores)));
    %Get the energy due to overlaps of the residue with non-moving atoms
    %from the last row of the coordinates
    min_energy = min_energy + current_coord(size(current_coord,1),1);
    
    %Loop over all other residues
    for rest_res  = loop_cores+1:num_res
        next_coord = eval(strcat('coor_m', num2str(rest_res)));
        new_energy = check_clash(current_coord(1:size(current_coord,1)-1,:), next_coord(1:size(next_coord,1)-1,:));
        min_energy = min_energy + new_energy;
    end
end


end
