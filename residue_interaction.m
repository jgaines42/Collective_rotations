%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [are_hitting] = residue_interaction(tempModel2,res1, res2,check_all, all_moving)
%
% Determines if two residues are able to hit each other
%
% Input:
%   tempModel2: Cell array of the entire protein
%   res1: index of first residue in all_moving
%   res2: index of second residue in all moving
%   check_all: if 0, just check initial positions for overlap, if 1, rotate
%       the 2 residues (allows for increased speed)
%   all_moving: All residue Ids that will be rotated in find_interactions
%
% Output:
%   Returns 1 if the residues can hit without hitting a backbone first.
%   Returns 0 if they cannont hit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [are_hitting] = residue_interaction(tempModel2,res1, res2,check_all, all_moving)
are_hitting = 0;

%Extract sizes and residue id for each atom
sizes_all = cell2mat(tempModel2(:,12));
res_ids = cell2mat(tempModel2(:,6));

%Extract first residue to look at
ind0 = find(res_ids == all_moving(res1));
resiName1 = tempModel2{ind0(1),4}; %Get name of interface resiude
resiName1(2:3) = lower(resiName1(2:3));
resiId1 = double(all_moving(res1));

switch (resiName1)
    case 'Ala'
        numAtom = 16;
        DOF = 1;
    case 'Arg'
        numAtom = 30;
        DOF = 3;
    case 'Asn'
        numAtom = 20;
        DOF = 2;
    case 'Asp'
        numAtom = 18;
        DOF = 2;
    case 'Cys'
        DOF = 1;
        numAtom = 17;
    case 'Gln'
        numAtom = 23;
        DOF = 3;
    case 'Glu'
        DOF = 3;
        numAtom = 21;
    case 'Gly'
        DOF = 0;
        numAtom = 13;
    case 'His'
        numAtom = 22;
        DOF = 2;
    case 'Ile'
        numAtom =25;
        DOF = 2;
    case 'Leu'
        numAtom =25;
        DOF = 2;
    case 'Lys'
        numAtom = 28;
        DOF = 4;
    case 'Met'
        numAtom = 23;
        DOF = 3;
    case 'Phe'
        numAtom = 26;
        DOF =2;
    case 'Pro'
        fprintf('Not yet supported\n' );
        DOF = 0;
    case 'Ser'
        numAtom = 17;
        DOF = 1;
    case 'Thr'
        numAtom = 20;
        DOF = 1;
    case 'Trp'
        numAtom = 30;
        DOF = 2;
    case 'Tyr'
        numAtom = 27;
        DOF = 2;
    case 'Val'
        numAtom = 22;
        DOF = 1;
    otherwise
        fprintf('Invalid amino acid\n' );
        DOF = 0;
end


%Isolate dipeptide
orig = [];
[allDipeptide1,next_pro] = isolate_dipeptide(tempModel2, res_ids, resiId1, resiName1);
if DOF >0 && (size(allDipeptide1,1) == numAtom || (size(allDipeptide1,1) == numAtom-1 && next_pro ==1))
    if next_pro == 0
        [new_Dipeptide, correct_now, ind] = check_Dipeptide_order(allDipeptide1(4:size(allDipeptide1,1)-3,:), resiName1);
        allDipeptide1(4:size(allDipeptide1,1)-3,:)= new_Dipeptide;
    else
        [new_Dipeptide, correct_now, ind] = check_Dipeptide_order(allDipeptide1(4:size(allDipeptide1,1)-2,:), resiName1);
        allDipeptide1(4:size(allDipeptide1,1)-2,:)= new_Dipeptide;
        allDipeptide1(size(allDipeptide1,1)+1,:) = allDipeptide1(size(allDipeptide1,1),:);
    end

    if correct_now == 1 && DOF >0
        
        %Extract second residue
        ind0 = find(res_ids == all_moving(res2));
        resiName2 = tempModel2{ind0(1),4}; %Get name of interface resiude
        resiName2(2:3) = lower(resiName2(2:3));
        resiId2 = double(all_moving(res2));
        
        switch (resiName2)
            case 'Ala'
                numAtom = 16;
                DOF = 1;
            case 'Arg'
                numAtom = 30;
                DOF = 3;
            case 'Asn'
                numAtom = 20;
                DOF = 2;
            case 'Asp'
                numAtom = 18;
                DOF = 2;
            case 'Cys'
                DOF = 1;
                numAtom = 17;
            case 'Gln'
                numAtom = 23;
                DOF = 3;
            case 'Glu'
                DOF = 3;
                numAtom = 21;
            case 'Gly'
                DOF = 0;
                numAtom = 13;
            case 'His'
                numAtom = 22;
                DOF = 2;
            case 'Ile'
                numAtom =25;
                DOF = 2;
            case 'Leu'
                numAtom =25;
                DOF = 2;
            case 'Lys'
                numAtom = 28;
                DOF = 4;
            case 'Met'
                numAtom = 23;
                DOF = 3;
            case 'Phe'
                numAtom = 26;
                DOF =2;
            case 'Pro'
                fprintf('Not yet supported\n' );
                DOF = 0;
            case 'Ser'
                numAtom = 17;
                DOF = 1;
            case 'Thr'
                numAtom = 20;
                DOF = 1;
            case 'Trp'
                numAtom = 30;
                DOF = 2;
            case 'Tyr'
                numAtom = 27;
                DOF = 2;
            case 'Val'
                numAtom = 22;
                DOF = 1;
            otherwise
                fprintf('Invalid amino acid\n' );
                DOF = 0;
        end

        
        
        %Isolate dipeptide
        orig = [];
        [allDipeptide2,next_pro] = isolate_dipeptide(tempModel2, res_ids, resiId2, resiName2);
        if DOF >0 && (size(allDipeptide2,1) == numAtom || (size(allDipeptide2,1) == numAtom-1 && next_pro ==1))
            
            if next_pro == 0
                [new_Dipeptide, correct_now, ind] = check_Dipeptide_order(allDipeptide2(4:size(allDipeptide2,1)-3,:), resiName2);
                allDipeptide2(4:size(allDipeptide2,1)-3,:)= new_Dipeptide;
            else
                [new_Dipeptide, correct_now, ind] = check_Dipeptide_order(allDipeptide2(4:size(allDipeptide2,1)-2,:), resiName2);
                allDipeptide2(4:size(allDipeptide2,1)-2,:)= new_Dipeptide;
                allDipeptide2(size(allDipeptide2,1)+1,:) = allDipeptide2(size(allDipeptide2,1),:);
                
            end
            
            %Get rest of protein backbone
            %Changed on 9/12 to include backbone of adjacent residues
            %(that was grouped in dipeptide before).
            ind0 = ismember(cell2mat(tempModel2(:,6)), [resiId1, resiId2]);
            ind1 = ismember(tempModel2(:,2), {'N', 'H', 'CA', 'C', 'O', 'HA'});
            rest_backbone = tempModel2((~ind0)&ind1,:);
            
            %% Save dipeptide positions
            
            %Check that Cbeta atoms (atom 8) are within 10A, if they are
            %further than that, the atoms can't hit
            res_clash = [8 8];
            dist1 = cell2mat(allDipeptide1(res_clash(:,1),8:10))-cell2mat(allDipeptide2(res_clash(:,2),8:10));
            dist2 = sum(dist1.^2,2);
            if min(dist2) < 10^2 % If within 10A, run check
                [are_hitting] = move_check_2(allDipeptide1,allDipeptide2,check_all, rest_backbone);
            else
                are_hitting = 0;
            end
        end
    end
end
