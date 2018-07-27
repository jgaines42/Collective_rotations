%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [energy, coordinates,dihedrals] = Rotate_residue_in_protein_energycuttoff(Position, res_name, Atom_sizes,  rest_of_pro, next_pro, energy_cutoff)
% 
% Input:
%   Position: x,y,z coordinates of the dipeptide
%   res_name: 3 letter name of the residues
%   Atom_size: Atomic radii of the atoms in the dipeptide
%   rest_of_pro: The rest of the protein (full data), minus other residues
%   being rotate
%   next_pro: 0 or 1 indicating if the next residue is a proline
%   energy_cutoff: Energy cutoff to use in calculations
%   
% Returns:
%   energy: energy of each configuraiton returned
%   coordinates: coordinates of the residue at all configurations with E <
%   energy_cutoff
%   dihedrals: Dihedral angles of each configuration. Includes CH3 and OH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [energy, coordinates,dihedrals] = Rotate_residue_in_protein_energycuttoff(Position, res_name, Atom_sizes,  rest_of_pro, next_pro, energy_cutoff)


rest_pro_position = cell2mat(rest_of_pro(:,8:10));
XtalPos = Position;

%Run setup script to get all of the variables for the particular amino acid
switch_residue_setup;
set_up_clashArrays;

numAtom = size(Position,1);

%set up initial variables to be returned
energy = zeros(10000,1);
coordinates = zeros(size(moveAtomID2,2),3*10000);
dihedrals = zeros(10000,DOF+CH3);
count_lowE = 1;

%Prepare to rotate Chi1
subtract_array_1 = repmat(Position(iChi1Array(2),:),numAtom,1);
delta_term_1 =  pi*sign(InitChi1)*InitChi1/180;

%Rotate chi1
for chi1 = 1:72
    Position=XtalPos;
    setChi1 = chi1*5;
    Position = Rotate_DA(Position, setChi1, subtract_array_1, delta_term_1, iChi1Array, moveAtomID2);
    [chi1_energy] = get_energy_wProtein(c1_Clash, Position,  protein_clash_c1, rest_pro_position);
    total_energy = chi1_energy;

    % If has Chi 2
    if DOF >=2 && chi1_energy < energy_cutoff
        Pos_b4_Chi2 = Position;
        subtract_array_2 = repmat(Position(iChi2Array(2),:),numAtom,1);
        delta_term_2 =  pi*sign(InitChi2)*InitChi2/180;
        %Rotate Chi 2
        for chi2 = 1:max_Chi2
            Position=Pos_b4_Chi2;
            setChi2 = chi2*5;
            Position = Rotate_DA(Position, setChi2, subtract_array_2, delta_term_2, iChi2Array, moveAtomID);
            [chi2_energy] = get_energy_wProtein(c2_Clash, Position,  protein_clash_c2, rest_pro_position);
            
            total_energy = chi1_energy + chi2_energy;

            % If has Chi 3
            if DOF >= 3 && total_energy < energy_cutoff
                Pos_b4_Chi3 = Position;
                subtract_array_3 = repmat(Position(iChi3Array(2),:),numAtom,1);
                delta_term_3 =  pi*sign(InitChi3)*InitChi3/180;
                
                % Rotate Chi 3
                for chi3 = 1:72
                    Position=Pos_b4_Chi3;
                    setChi3 = chi3*5;
                    Position = Rotate_DA(Position, setChi3, subtract_array_3, delta_term_3, iChi3Array, moveAtomID3);
                    [chi3_energy] = get_energy_wProtein(c3_Clash, Position,  protein_clash_c3, rest_pro_position);
                    total_energy = chi1_energy + chi2_energy + chi3_energy;                    
                    
                    % If lower than energy cutoff
                    if DOF == 3 && total_energy < energy_cutoff
                        index = (chi1-1)*72*72 + (chi2-1)*72 + chi3;
                        
                        % If residue has CH3 group
                        if CH3 == 1
                            subtract_array_HG1 = repmat(Position(HG_Array_1(2),:),numAtom,1);
                            delta_term_HG1 =  pi*sign(InitChi_HG1)*InitChi_HG1/180;
                            min_energy_CH3 = 0;
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % rotate_CH3_group_1_LJ_protein_returnAll uses Position,subtract_array_HG1, delta_term_HG1, HG_Array_1,
                            % moveAtomID_HG1, HG1_Clash, numAtom, rest_pro_position,protein_clash_HG1)
                            % and sets coordinates, energy and dihedrals;
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            rotate_CH3_group_1_LJ_protein_returnAll;
                           
                        else % If no CH3 group, save the data
                            energy(count_lowE,1) = total_energy;
                            coordinates(:,(count_lowE-1)*3+1:count_lowE*3) = Position(moveAtomID2,:);
                            dihedrals(count_lowE,1:3) = [setChi1, setChi2, setChi3];
                            count_lowE = count_lowE + 1;
                            if size(dihedrals,1) == (count_lowE -1)
                                energy = [energy; zeros(10000,1)];
                                coordinates = [coordinates,zeros(size(moveAtomID2,2),3*10000)];
                                dihedrals = [dihedrals;zeros(10000,DOF+CH3)];
                            end
                        end                     
                    end
                end%Chi3 loop
            end
            
            % If there wasn't a Chi 3, then do more Chi 2 stuff
            
            if DOF == 2 && total_energy < energy_cutoff 
                index = (chi1-1)*72 + chi2;
                
                if CH3 == 1 % If 1 CH3 group
                    
                    subtract_array_HG1 = repmat(Position(HG_Array_1(2),:),numAtom,1);
                    delta_term_HG1 =  pi*sign(InitChi_HG1)*InitChi_HG1/180;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % rotate_CH3_group_1_LJ_protein_returnAll uses Position,subtract_array_HG1, delta_term_HG1, HG_Array_1,
                    % moveAtomID_HG1, HG1_Clash, numAtom, rest_pro_position,protein_clash_HG1)
                    % and sets coordinates, energy and dihedrals;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    rotate_CH3_group_1_LJ_protein_returnAll;
                    
                   
                elseif CH3 == 2 % If 2 CH3 groups
                    
                    subtract_array_HG1 = repmat(Position(HG_Array_1(2),:),numAtom,1);
                    delta_term_HG1 =  pi*sign(InitChi_HG1)*InitChi_HG1/180;
                    min_energy_CH3 =0;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % rotate_CH3_group_2_LJ_protein_returnAll uses a bunch
                    % of stuff
                    % and sets coordinates, energy and dihedrals;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    rotate_CH3_group_2_LJ_protein_returnAll;
                                      
                elseif DOF == 2 && OH == 1 && CH3 == 0 && total_energy < energy_cutoff
                    % This runs for Tyr
                    Pos_b4_OH = Position;
                    subtract_array_OH = repmat(Position(iOHArray(2),:),numAtom,1);
                    delta_term_OH =  pi*sign(InitOH)*InitOH/180;
                    
                    % Loop OH to find all positions with energy below the cutoff
                    OH_val = 0;
                    while  OH_val < 72
                        OH_val = OH_val +1;
                        setOH = OH_val*5;
                        Position = Pos_b4_OH;
                        Position = Rotate_DA(Position, setOH, subtract_array_OH, delta_term_OH, iOHArray, moveAtomOH);
                        [OH_energy] = get_energy_wProtein(OH_Clash, Position,  protein_clash_OH, rest_pro_position);
                        
                        if total_energy + OH_energy < energy_cutoff
                            energy(count_lowE,1) = total_energy+ OH_energy;
                            coordinates(:,(count_lowE-1)*3+1:count_lowE*3) = Position(moveAtomID2,:);
                            dihedrals(count_lowE,1:3) = [setChi1,setChi2, setOH];
                            count_lowE = count_lowE+1;
                            if size(dihedrals,1) == (count_lowE -1)
                                energy = [energy; zeros(10000,1)];
                                coordinates = [coordinates,zeros(size(moveAtomID2,2),3*10000)];
                                dihedrals = [dihedrals;zeros(10000,DOF+CH3+OH)];
                            end
                        end
                    end
                    %%
                else %If no CH3 groups and no OH groups (Phe, etc)
                    if total_energy < energy_cutoff
                        energy(count_lowE,1) = total_energy;
                        coordinates(:,(count_lowE-1)*3+1:count_lowE*3) = Position(moveAtomID2,:);
                        dihedrals(count_lowE,1:2) = [setChi1, setChi2];
                        count_lowE = count_lowE+1;
                        if size(dihedrals,1) == (count_lowE -1)
                            energy = [energy; zeros(10000,1)];
                            coordinates = [coordinates,zeros(size(moveAtomID2,2),3*10000)];
                            dihedrals = [dihedrals;zeros(10000,DOF+CH3)];
                        end
                    end
                end
            end
            
        end
    end
    
    %% If no Chi 2
    
    %If 1 CH3 group
    if CH3 == 1 &&  OH == 0 && DOF == 1 &&  chi1_energy < energy_cutoff
        
        subtract_array_HG1 = repmat(Position(HG_Array_1(2),:),numAtom,1);
        delta_term_HG1 =  pi*sign(InitChi_HG1)*InitChi_HG1/180;
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % rotate_CH3_group_1_LJ_protein_returnAll uses Position,subtract_array_HG1, delta_term_HG1, HG_Array_1,
        % moveAtomID_HG1, HG1_Clash, numAtom, rest_pro_position,protein_clash_HG1)
        % and sets coordinates, energy and dihedrals;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rotate_CH3_group_1_LJ_protein_returnAll;
    
    % If 2 CH3 groups
    elseif OH == 0 && DOF == 1 && CH3 == 2 && chi1_energy < energy_cutoff 
        subtract_array_HG1 = repmat(Position(HG_Array_1(2),:),numAtom,1);
        delta_term_HG1 =  pi*sign(InitChi_HG1)*InitChi_HG1/180;
       
        total_energy = chi1_energy;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % rotate_CH3_group_2_LJ_protein_returnAll uses a bunch
        % of stuff
        % and sets coordinates, energy and dihedrals;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rotate_CH3_group_2_LJ_protein_returnAll;

     % If OH, but no CH3   
    elseif DOF == 1 && OH == 1 && CH3 == 0 && chi1_energy < energy_cutoff
        % This runs for Ser
        total_energy = chi1_energy;
        Pos_b4_OH = Position;
        subtract_array_OH = repmat(Position(iOHArray(2),:),numAtom,1);
        delta_term_OH =  pi*sign(InitOH)*InitOH/180;
        
        % Loop OH to find all positions with energy below the cutoff
        OH_val = 0;
        while  OH_val < 72
            OH_val = OH_val +1;
            setOH = OH_val*5;
            Position = Pos_b4_OH;
            Position = Rotate_DA(Position, setOH, subtract_array_OH, delta_term_OH, iOHArray, moveAtomOH);
            [OH_energy] = get_energy_wProtein(OH_Clash, Position,  protein_clash_OH, rest_pro_position);
            
            if total_energy + OH_energy < energy_cutoff
                energy(count_lowE,1) = total_energy+ OH_energy;
                coordinates(:,(count_lowE-1)*3+1:count_lowE*3) = Position(moveAtomID2,:);
                dihedrals(count_lowE,1:2) = [setChi1, setOH];
                count_lowE = count_lowE+1;
                if size(dihedrals,1) == (count_lowE -1)
                    energy = [energy; zeros(10000,1)];
                    coordinates = [coordinates,zeros(size(moveAtomID2,2),3*10000)];
                    dihedrals = [dihedrals;zeros(10000,DOF+CH3+OH)];
                end
            end
        end
        
     % If CH3 and OH   
    elseif  chi1_energy < energy_cutoff && DOF == 1 && OH == 1 && CH3 == 1
        total_energy = chi1_energy;
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % rotate_CH3_OH_group_LJ_returnAll uses a bunch
        % of stuff
        % and sets coordinates, energy and dihedrals;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rotate_CH3_OH_group_LJ_returnAll;

    end
    
    
    if DOF == 1 && CH3 == 0 && OH == 0
        % If CH3 has a value then this section has already been done
        
        
        if chi1_energy < energy_cutoff
            energy(count_lowE,1) = total_energy;
            coordinates(:,(count_lowE-1)*3+1:count_lowE*3) = Position(moveAtomID2,:);
            dihedrals(count_lowE,1) = [setChi1];
            count_lowE = count_lowE+1;
            if size(dihedrals,1) == (count_lowE -1)
                energy = [energy; zeros(10000,1)];
                coordinates = [coordinates,zeros(size(moveAtomID2,2),3*10000)];
                dihedrals = [dihedrals;zeros(10000,DOF+CH3)];
            end
        end
        
        
    end
    
end

coordinates = coordinates(:,1:(count_lowE-1)*3);
dihedrals = dihedrals(1:(count_lowE-1),:);
energy = energy(1:(count_lowE-1),:);

coordinates = [coordinates,Atom_sizes(moveAtomID2)];


end