load(strcat(folder_pdb,pdb_names, '_all_core_paired_data_2.mat'));
ind0 = new_moving(:,2) == group;
sub_paired = paired(ind0,ind0);
all_res_ind = [1:max_depth];

% Only does 4 paired interactions 
num_tested = 0;
for res1 = 1:size(sub_paired,2)
    for res2 = res1+1:size(sub_paired,2)
        if sub_paired(res1,res2) == 1 && not_found == 0 && num_tested < 4
            num_tested = num_tested + 1;
            ind0 = ismember(all_res_ind,[res1,res2]);
            rest_res = all_res_ind(~ind0);
            coord_x = eval(strcat('coord_', num2str(res1)));
            coord_y = eval(strcat('coord_', num2str(res2)));
            
            e_x = coord_x(size(coord_x,1), [1:3:size(coord_x,2)-1]);
            e_y = coord_y(size(coord_y,1), [1:3:size(coord_y,2)-1]);
            
            %Increase the cutoff by a tiny bit, so that you can see what give you same
            %or less than
            mins_x = repmat(min_sofar*1.01,size(e_x,2),1);
            mins_y = repmat(min_sofar*1.01,size(e_y,2),1);
            
            en_loc_x = size(coord_x,1);
            en_loc_y = size(coord_y,1);
            % Loop over all combinations
            tic
            
            x_ind = 1:size(e_x,2);
            
            indx = zeros(size(e_x,2),1);
            for i = 1:size(e_x,2)
                e_sub = e_x(i)+e_y +sum(mins(rest_res));
                indx(i) =  min(e_sub) <= min_sofar;
                %   indx(i) = 1;
                % end
            end
            sub_x = x_ind(logical(indx));
            
           
            counter = 0;
            y_ind = 1:size(e_y,2);
            for i = sub_x%1:size(e_x,2)
                counter = counter + 1;
                this_x = coord_x(1:en_loc_x-1,[(i-1)*3+1:i*3, size(coord_x,2)]);
                nn = e_x(i) + e_y+sum(mins(rest_res));
                if min(nn) <= min_sofar
                    ind0 = nn <= min_sofar;
                    sub_y = y_ind(ind0');
                    for j = sub_y
                        this_y = coord_y(1:en_loc_y-1,[(j-1)*3+1:j*3, size(coord_y,2)]);
                        this_e = e_x(i) + e_y(j) + sum(mins(rest_res));
                        if this_e <= min_sofar
                            this_e = this_e + check_clash(this_x, this_y);
                        end
                        if this_e <= mins_x(i)
                            mins_x(i) = this_e;
                        end
                        if this_e <= mins_y(j)
                            mins_y(j) = this_e;
                        end
                    end
                end
                if mod(counter,1000) == 0
                    counter
                end
            end
            
            %Make changes to the coord
            current_coord = eval(strcat('coord_', num2str(res1)));
            ind1 = find(mins_x' <= min_energy);
            ind2 = sort([(ind1-1)*3+1, (ind1-1)*3+2, (ind1-1)*3+3]);
            load(strcat(folder_data, pdb_names, '_', Residue_names{res1}, num2str(Residue_id(res1)), '_dihedrals_fast2_', num2str(variant), '.mat'));
            dihedrals = dihedrals(ind1,:);
            if size(dihedrals,1) == 0
                not_found = 1;
            end
            save(strcat(folder_data, pdb_names, '_', Residue_names{res1}, num2str(Residue_id(res1)), '_dihedrals_fast2_', num2str(variant), '.mat'), 'dihedrals');

            switch res1
                case 1
                    coord_1 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                case 2
                    coord_2 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                    
                case 3
                    coord_3 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                case 4
                    coord_4 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                case 5
                    coord_5 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                case 6
                    coord_6 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                case 7
                    coord_7 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                case 8
                    coord_8 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                case 9
                    coord_9 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                case 10
                    coord_10 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
            end
            %Make changes to the coord
            current_coord = eval(strcat('coord_', num2str(res2)));
            ind1 = find(mins_y' <= min_energy);
            ind2 = sort([(ind1-1)*3+1, (ind1-1)*3+2, (ind1-1)*3+3]);
            load(strcat(folder_data, pdb_names, '_', Residue_names{res2}, num2str(Residue_id(res2)), '_dihedrals_fast2_', num2str(variant), '.mat'));
            dihedrals = dihedrals(ind1,:);
            if size(dihedrals,1) == 0
                not_found = 1;
            end
            save(strcat(folder_data, pdb_names, '_', Residue_names{res2}, num2str(Residue_id(res2)), '_dihedrals_fast2_', num2str(variant), '.mat'), 'dihedrals');
            switch res2
                case 1
                    coord_1 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                case 2
                    coord_2 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                    
                case 3
                    coord_3 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                case 4
                    coord_4 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                case 5
                    coord_5 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                case 6
                    coord_6 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                case 7
                    coord_7 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                case 8
                    coord_8 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                case 9
                    coord_9 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
                case 10
                    coord_10 = [current_coord(:,ind2),current_coord(:,size(current_coord,2))];
            end
        end
    end
end

%% Reassess Mins
if not_found == 0
    for res1 = 1:size(mins,2)
        coord_x = eval(strcat('coord_', num2str(res1)));
        e_x = coord_x(size(coord_x,1), [1:3:size(coord_x,2)-1]);
        
        mins(res1) = min(e_x);
    end
end
