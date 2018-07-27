%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [overall_min,best_loc] = set_energy_cutoff_3(num_res, overall_min)
%
% Randomly samples combinations of the residues to try to find a lower
% energy
%
% Input: 
%   num_res: Number of residues
%   overall_min: The lowest energy state found so far
%
% Output:
%   overall_min: The lowest energy state foun
%   best_loc: The index of the dihedral angles of the lowest energy state
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [overall_min,best_loc] = set_energy_cutoff_3(num_res, overall_min)
global coord_1 coord_2 coord_3 coord_4 coord_5 coord_6 coord_7 coord_8 coord_9 coord_10;

best_loc = zeros(1,num_res);
loc = zeros(1,num_res);
num_rand = round(100000/num_res); %The number of random states to sample

% Calculate the number of possible permutations
total_sample = 1;
for residues = 1:num_res
    current_coor = eval(strcat('coord_', num2str(residues)));
    total_sample = total_sample*(size(current_coor,2)-1)/3;
end

% If the sqrt of the total permutations is < num_rand, then reduce num_rand
% to sqrt(total_sample), this is to keep from sampling excessively
if sqrt(total_sample) < num_rand
    num_rand = round(sqrt(total_sample));
end

%Generate a bunch of random numbers
for residues = 1:num_res
    switch residues
        case 1
            r1 = randi([1 (size(coord_1,2)-1)/3],1,num_rand);
        case 2
            r2 = randi([1 (size(coord_2,2)-1)/3],1,num_rand);
        case 3
            r3 = randi([1 (size(coord_3,2)-1)/3],1,num_rand);
        case 4
            r4 = randi([1 (size(coord_4,2)-1)/3],1,num_rand);
        case 5
            r5 = randi([1 (size(coord_5,2)-1)/3],1,num_rand);
        case 6
            r6 = randi([1 (size(coord_6,2)-1)/3],1,num_rand);
        case 7
            r7 = randi([1 (size(coord_7,2)-1)/3],1,num_rand);
        case 8
            r8 = randi([1 (size(coord_8,2)-1)/3],1,num_rand);
        case 9
            r9 = randi([1 (size(coord_9,2)-1)/3],1,num_rand);
        case 10
            r10 = randi([1 (size(coord_10,2)-1)/3],1,num_rand);
    end
end

% Sample random dihedral angles
for random_dihedrals = 1:num_rand
    for residues = 1:num_res
        
        %Get the coordinates of this randomly chosen state
        switch residues
            case 1
                loc(1,residues) = r1(random_dihedrals);
                coor_m1 = coord_1(:,[(loc(1,residues)-1)*3+1:loc(1,residues)*3,size(coord_1,2)]);
            case 2
                loc(1,residues) = r2(random_dihedrals);
                coor_m2 = coord_2(:,[(loc(1,residues)-1)*3+1:loc(1,residues)*3,size(coord_2,2)]);
            case 3
                loc(1,residues) = r3(random_dihedrals);
                coor_m3 = coord_3(:,[(loc(1,residues)-1)*3+1:loc(1,residues)*3,size(coord_3,2)]);
            case 4
                loc(1,residues) = r4(random_dihedrals);
                coor_m4 = coord_4(:,[(loc(1,residues)-1)*3+1:loc(1,residues)*3,size(coord_4,2)]);
            case 5
                loc(1,residues) = r5(random_dihedrals);
                coor_m5 = coord_5(:,[(loc(1,residues)-1)*3+1:loc(1,residues)*3,size(coord_5,2)]);
            case 6
                loc(1,residues) = r6(random_dihedrals);
                coor_m6 = coord_6(:,[(loc(1,residues)-1)*3+1:loc(1,residues)*3,size(coord_6,2)]);
            case 7
                loc(1,residues) = r7(random_dihedrals);
                coor_m7 = coord_7(:,[(loc(1,residues)-1)*3+1:loc(1,residues)*3,size(coord_7,2)]);
            case 8
                loc(1,residues) = r8(random_dihedrals);
                coor_m8 = coord_8(:,[(loc(1,residues)-1)*3+1:loc(1,residues)*3,size(coord_8,2)]);
            case 9
                loc(1,residues) = r9(random_dihedrals);
                coor_m9 = coord_9(:,[(loc(1,residues)-1)*3+1:loc(1,residues)*3,size(coord_9,2)]);
            case 10
                loc(1,residues) = r10(random_dihedrals);
                coor_m10 = coord_10(:,[(loc(1,residues)-1)*3+1:loc(1,residues)*3,size(coord_10,2)]);
        end
        
    end
    
    % Do pairwise energy calculations for all residues in the randomly
    % chosen state
    min_energy = 0;
    loop_all = 1;
    for loop_cores = 1:num_res
        current_coord = eval(strcat('coor_m', num2str(loop_cores)));
        min_energy = min_energy + current_coord(size(current_coord,1),1);
        if min_energy < overall_min
            for rest_res  = loop_cores+1:num_res
                next_coord = eval(strcat('coor_m', num2str(rest_res)));
                new_energy = check_clash(current_coord(1:size(current_coord,1)-1,:), next_coord(1:size(next_coord,1)-1,:));
                min_energy = min_energy + new_energy;
            end
        else
            loop_all = 0;
        end
    end
    % If the lowest energy is less than the overall_min, save the value
    if loop_all == 1 & min_energy < overall_min
        overall_min = min_energy;
        best_loc = loc;
    end
    
end
end