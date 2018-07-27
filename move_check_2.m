%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [are_hitting] = move_check_2(allDipeptide1,allDipeptide2,check_all, rest_backbone)
%
% Determines if two residues can hit each other
%
% Input:
%   allDipeptide1: Dipeptide of first amino acid
%   allDipeptide2: Dipeptide of second amino acid
%   check_all: 0 if just checking for clashes without rotating
%              1 if rotating the amino acids
%   rest_backbone: Matrix containing the rest of the protein backbone
%
% Output:
%   are_hitting: 0 if amino acids can't hit without hitting the backbone
%               1 if amino acids can hit without hitting the backbone
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [are_hitting] = move_check_2(allDipeptide1,allDipeptide2,check_all, rest_backbone)

overlap_amount = 0;
Back_Pos = cell2mat(rest_backbone(:,8:10));
%Extract Residue names
resiName1 = allDipeptide1{4,4};
resiName1(2:3) = lower(resiName1(2:3));
resiName2 = allDipeptide2{4,4};
resiName2(2:3) = lower(resiName2(2:3));


%Step size will be 10 degrees
max_Chi2_1 = 36;
max_Chi2_2 = 36;
if resiName1 == 'Phe'
    max_Chi2_1 = 36/2;
end
if resiName2 == 'Phe'
    max_Chi2_2 = 36/2;
end

resiID2 = cell2mat(allDipeptide2(4,6));
Position_1 = cell2mat(allDipeptide1(:,8:10));
Position_2 = cell2mat(allDipeptide2(:,8:10));
Xtal1 = Position_1;
Xtal2 = Position_2;

%% Set up Residue 1
Amino_Acid_1 = switch_residue_setup(resiName1);
%We don't need the clash array so just pass empty arrays
[Init_dihedrals_1, ~] = set_up_clashArrays_dipeptide( resiName1,0,zeros(size(Position_1,1),1), Amino_Acid_1, Position_1);
numAtom1 =Amino_Acid_1.numAtom;
DOF1 =Amino_Acid_1.DOF;
all_move_1 = Amino_Acid_1.moveAtomID2; %Which atoms will move

%% Set up Residue 2
Amino_Acid_2 = switch_residue_setup(resiName2);
%We don't need the clash array so just pass empty arrays
[Init_dihedrals_2, ~] = set_up_clashArrays_dipeptide( resiName2,0,zeros(size(Position_2,1),1), Amino_Acid_2, Position_2);
numAtom2 =Amino_Acid_2.numAtom;
DOF2 =Amino_Acid_2.DOF;
all_move_2 = Amino_Acid_2.moveAtomID2;

%Check that you have the right number of atoms
if size(Position_1,1) ~= numAtom1 || size(Position_2,1) ~= numAtom2
    are_hitting = -1;
else
    
    %Make clash list
    res_clash= [];
    for clash_res = 1:size(all_move_2,2)
        res_clash = [res_clash;all_move_1', repmat(all_move_2(clash_res),size(all_move_1,2),1)];
    end
    res_clash(:,3) = cell2mat(allDipeptide1(res_clash(:,1),12))  + cell2mat(allDipeptide2(res_clash(:,2),12));
    res_clash(:,4) = res_clash(:,3).^2;
    

    
    %% Run rotations
    subtract_array_1_1 = repmat(Position_1(Amino_Acid_1.iChi1Array(2),:),numAtom1,1);
    delta_term_1_1 =  pi*sign(Init_dihedrals_1.InitChi1)*Init_dihedrals_1.InitChi1/180;
    subtract_array_1_2 = repmat(Position_2(Amino_Acid_2.iChi1Array(2),:),numAtom2,1);
    delta_term_1_2 =  pi*sign(Init_dihedrals_2.InitChi1)*Init_dihedrals_2.InitChi1/180;
    
    if check_all == 1
        %Make clash list for backbone
        back_clash1 = zeros(size(all_move_1,2)*size(rest_backbone,1),4);
        for clash_res = 1:size(all_move_1,2)
            back_clash1((clash_res-1)*size(rest_backbone,1)+1 : clash_res*size(rest_backbone,1),1:2) = [repmat(all_move_1(clash_res),size(rest_backbone,1),1),[1:size(rest_backbone,1)]'];
        end
        back_clash1(:,3) = cell2mat(allDipeptide1(back_clash1(:,1),12))+cell2mat(rest_backbone(back_clash1(:,2),12));
        back_clash1(:,4) = back_clash1(:,3).^2;
        
        %Limit to atoms within 10A of any atom
        dist = Position_1(back_clash1(:,1),1:3)-Back_Pos(back_clash1(:,2),1:3);
        distemp = sum(dist.^2,2);
        clash = distemp - back_clash1(:,4);
        ind0 = clash < 100;
        x = back_clash1(ind0,:);
        y = unique(x(:,2));
        ind0 = ismember(back_clash1(:,2),y);
        back_clash1 = back_clash1(ind0,:);
        
        back_clash2 = zeros(size(all_move_2,2)*size(rest_backbone,1),4);
        for clash_res = 1:size(all_move_2,2)
            back_clash2((clash_res-1)*size(rest_backbone,1)+1 : clash_res*size(rest_backbone,1),1:2) = [repmat(all_move_2(clash_res),size(rest_backbone,1),1),[1:size(rest_backbone,1)]'];
        end
        back_clash2(:,3) = cell2mat(allDipeptide2(back_clash2(:,1),12))+cell2mat(rest_backbone(back_clash2(:,2),12));
        back_clash2(:,4) = back_clash2(:,3).^2;
        dist = Position_2(back_clash2(:,1),1:3)-Back_Pos(back_clash2(:,2),1:3);
        distemp = sum(dist.^2,2);
        clash = distemp - back_clash2(:,4);
        ind0 = clash < 100;
        x = back_clash2(ind0,:);
        y = unique(x(:,2));
        ind0 = ismember(back_clash2(:,2),y);
        back_clash2 = back_clash2(ind0,:);
        
        
        %% Do single residue of each, get sublist of things to check
        D1 = [];
        iChi1Array1 = Amino_Acid_1.iChi1Array;
        moveAtomID2_1 = Amino_Acid_1.moveAtomID2;
        if DOF1 >=2
            iChi2Array1 = Amino_Acid_1.iChi2Array;
            moveAtomID_1 = Amino_Acid_1.moveAtomID;
            delta_term_2_1 =  pi*sign(Init_dihedrals_1.InitChi2)*Init_dihedrals_1.InitChi2/180;
        end
        if DOF1 >= 3
            iChi3Array1 = Amino_Acid_1.iChi3Array;
            moveAtomID3_1 = Amino_Acid_1.moveAtomID3;
            delta_term_3_1 =  pi*sign(Init_dihedrals_1.InitChi3)*Init_dihedrals_1.InitChi3/180;
        end
        
         for chi1_1 = 1:36 % Rotate Chi1 of Residue 1
            Position_1 = Xtal1;
            setChi1 = chi1_1*10;
            Position_1 = Rotate_DA(Position_1, setChi1, subtract_array_1_1, delta_term_1_1, iChi1Array1, moveAtomID2_1);
            if DOF1 == 1
               energy = get_energy_wProtein([], Position_1,  back_clash1(:,[1,2,4]), Back_Pos);
                    if energy == 0
                      D1 = [D1;setChi1];
                  end
            elseif DOF1 >= 2
                Pos_b4_Chi2_1 = Position_1;
                subtract_array_2_1 = repmat(Position_1(iChi2Array1(2),:),numAtom1,1);
                
                for chi2_1 = 1:max_Chi2_1
                    Position_1 = Pos_b4_Chi2_1;
                    setChi2 = chi2_1*10;
                    Position_1 = Rotate_DA(Position_1, setChi2, subtract_array_2_1, delta_term_2_1, iChi2Array1, moveAtomID_1);
                    if DOF1 == 2
                        %Check for backbone clash of residue 1
                        energy = get_energy_wProtein([], Position_1,  back_clash1(:,[1,2,4]), Back_Pos);
                        if energy == 0
                            D1 = [D1;setChi1,setChi2];
                        end
                    else
                        Pos_b4_Chi3 = Position_1;
                         subtract_array_3_1 = repmat(Position_1(iChi3Array1(2),:),numAtom1,1);
                         
                         for chi3_1 = 1:36 %Rotate Chi3 of Residue 2
                             Position_1 = Pos_b4_Chi3;
                             setChi3 = chi3_1*10;
                             Position_1 = Rotate_DA(Position_1, setChi3, subtract_array_3_1, delta_term_3_1, iChi3Array1, moveAtomID3_1);
                             
                             energy = get_energy_wProtein([], Position_1,  back_clash1(:,[1,2,4]), Back_Pos);
                             if energy == 0
                                 D1 = [D1;setChi1, setChi2, setChi3];
                             end
                         end
                    end
                        
                end
            end
         end
         
         D2 = [];
        iChi1Array2 = Amino_Acid_2.iChi1Array;
        moveAtomID2_2 = Amino_Acid_2.moveAtomID2;
        if DOF2 >=2
            iChi2Array2 = Amino_Acid_2.iChi2Array;
            moveAtomID_2 = Amino_Acid_2.moveAtomID;
            delta_term_2_2 =  pi*sign(Init_dihedrals_2.InitChi2)*Init_dihedrals_2.InitChi2/180;
        end
        if DOF2 >= 3
            iChi3Array2 = Amino_Acid_2.iChi3Array;
            moveAtomID3_2 = Amino_Acid_2.moveAtomID3;
            delta_term_3_2 =  pi*sign(Init_dihedrals_2.InitChi3)*Init_dihedrals_2.InitChi3/180;
        end
        
         for chi1_2 = 1:36 %Rotate Chi1 of Residue 2
             Position_2 = Xtal2;
             setChi1 = chi1_2*10;
             Position_2 = Rotate_DA(Position_2, setChi1, subtract_array_1_2, delta_term_1_2, iChi1Array2, moveAtomID2_2);
             
             if DOF2 == 1 % Residue 2 only has 1 DOF
                 % Check for backbone clash of residue 2
                 energy = get_energy_wProtein([], Position_2,  back_clash2(:,[1,2,4]), Back_Pos);
                 if energy == 0
                     D2 = [D2;setChi1];
                 end
             else
                 Pos_b4_Chi2 = Position_2;
                 subtract_array_2_2 = repmat(Position_2(iChi2Array2(2),:),numAtom2,1);
                 
                 for chi2_2 = 1:max_Chi2_2 %Rotate Residue 2 Chi2
                     Position_2 = Pos_b4_Chi2;
                     setChi2 = chi2_2*10;
                     Position_2 = Rotate_DA(Position_2, setChi2, subtract_array_2_2, delta_term_2_2, iChi2Array2, moveAtomID_2);
                     if DOF2 == 2 %If Residue 2 only has 2 DOF
                         %Check for backbone clash of Residue 2
                         energy = get_energy_wProtein([], Position_2,  back_clash2(:,[1,2,4]), Back_Pos);
                         if energy == 0
                             D2 = [D2;setChi1,setChi2];
                         end
                     else %DOF == 3
                         Pos_b4_Chi3 = Position_2;
                         subtract_array_3_2 = repmat(Position_2(iChi3Array2(2),:),numAtom2,1);
                         
                         for chi3_2 = 1:36 %Rotate Chi3 of Residue 2
                             Position_2 = Pos_b4_Chi3;
                             setChi3 = chi3_2*10;
                             Position_2 = Rotate_DA(Position_2, setChi3, subtract_array_3_2, delta_term_3_2, iChi3Array2, moveAtomID3_2);
                             
                             %Check for backbone clash of Residue 2
                             energy = get_energy_wProtein([], Position_2,  back_clash2(:,[1,2,4]), Back_Pos);
                             if energy == 0
                                 D2 = [D2;setChi1,setChi2, setChi3];
                             end
                         end
                     end
                 end
             end
         end

        %% Now loop through combined rotations to find clashes
        are_hitting = 0;
        for loop_1 = 1:size(D1,1)
            
            %Set first residue to a new position
            Position_1 = Xtal1;
            setChi1 = D1(loop_1,1);
            Position_1 = Rotate_DA(Position_1, setChi1, subtract_array_1_1, delta_term_1_1, iChi1Array1, moveAtomID2_1);
            if DOF1 >1
                setChi2 = D1(loop_1,2);
                subtract_array_2_1 = repmat(Position_1(iChi2Array1(2),:),numAtom1,1);
                Position_1 = Rotate_DA(Position_1, setChi2, subtract_array_2_1, delta_term_2_1, iChi2Array1, moveAtomID_1);
                if DOF1 > 2
                    setChi3 = D1(loop_1,3);
                    subtract_array_3_1 = repmat(Position_1(iChi3Array1(2),:),numAtom1,1);
                    Position_1 = Rotate_DA(Position_1, setChi3, subtract_array_3_1, delta_term_3_1, iChi3Array1, moveAtomID3_1);
                end
            end
            
            %Set second residue to a new position
            for loop_2 = 1:size(D2,1)
                Position_2 = Xtal2;
                setChi1 = D2(loop_2,1);
                Position_2 = Rotate_DA(Position_2, setChi1, subtract_array_1_2, delta_term_1_2, iChi1Array2, moveAtomID2_2);
                if DOF2 >1
                    setChi2 = D2(loop_2,2);
                    subtract_array_2_2 = repmat(Position_2(iChi2Array2(2),:),numAtom2,1);
                    Position_2 = Rotate_DA(Position_2, setChi2, subtract_array_2_2, delta_term_2_2, iChi2Array2, moveAtomID_2);
                    if DOF2 > 2
                        setChi3 = D2(loop_2,3);
                        subtract_array_3_2 = repmat(Position_2(iChi3Array2(2),:),numAtom2,1);
                        Position_2 = Rotate_DA(Position_2, setChi3, subtract_array_3_2, delta_term_3_2, iChi3Array2, moveAtomID3_2);
                    end
                end
                
                %Calculate energy
                dist = Position_1(res_clash(:,1),1:3)-Position_2(res_clash(:,2),1:3);
                distemp = sum(dist.^2,2);
                clash = distemp - res_clash(:,4);
                ind_clash = find(clash < 0);
                if size(ind_clash,1) > 0
                    are_hitting = 1;
                    break;
                end
            end
            if are_hitting == 1
                break
            end
        end
    else %just check initial
        dist = Position_1(res_clash(:,1),1:3)-Position_2(res_clash(:,2),1:3);
        distemp = sum(dist.^2,2);
        clash = distemp - res_clash(:,4);
        ind_clash = find(clash < 0);
        if size(ind_clash,1) > 0
            are_hitting = 1;
        else
            are_hitting = 0;
            
        end
    end
end
end
          
           