%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [all_coordinates, dihedrals] = set_up_recursion_SingleRes_rot(resiId, resiName, tempModel2, all_moving, energy_cutoff)
%
% Finds all states of the given residue with an energy <= energy_cutoff.
% Looks for clashes within the peptides as well as with all non-moving
% atoms in the protein
%
% Input:
%   resiId: Residue ID
%   resiName: 3 letter code
%   tempModel2: cell array of entire protein
%   all_moving: Residue Ids of all moving residues in the cluster
%   energy_cutoff
%
% Output: 
%   all_coordinates: xyz coordinates for all dihedral placements <=
%       energy_cutoff
%   dihedrals: The dihedral angles corresponding to the coordinates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [all_coordinates, dihedrals] = set_up_recursion_SingleRes_rot(resiId, resiName, tempModel2, all_moving, energy_cutoff)

ind0 = ismember(all_moving, resiId);
rest_of_moving = all_moving(~ind0);
res_ids = cell2mat(tempModel2(:,6)); %Ids of the whole protein

switch (resiName)
    case 'Ala'
        numAtom = 16;
        DOF = 1;
    case 'Ile'
        numAtom =25;
        DOF = 2;
        
    case 'Leu'
        numAtom =25;
        DOF = 2;
        
    case 'Val'
        numAtom = 22;
        DOF = 1;
        
    case 'Phe'
        numAtom = 26;
        DOF =2;
        
    case 'Trp'
        numAtom = 30;
        DOF = 2;
    case 'Tyr'
        numAtom = 27;
        DOF = 2;
    case 'Asn'
        fprintf('Not yet supported\n' );
        DOF = 0;
    case 'Cys'
        DOF = 1;
        numAtom = 17;
        
    case 'Glu'
        DOF = 3;
        numAtom = 21;
        DOF = 0;
    case 'Met'
        numAtom = 23;
        DOF = 3;
        
    case 'Ser'
        numAtom = 17;
        DOF = 1;
    case 'Thr'
        numAtom = 20;
        DOF = 1;
    case 'Asp'
        numAtom = 18;
        DOF = 2;
        DOF = 0;
        
    case 'Gln'
        fprintf('Not yet supported\n' );
        DOF = 0;
    case 'Arg'
        fprintf('Not yet supported\n' );
        DOF = 0;
    case 'His'
        numAtom = 22;
        DOF = 2;
    case 'Lys'
        numAtom = 28;
        DOF = 4;
        DOF = 0;
    case 'Gly'
        DOF = 0;
        fprintf('Not yet supported\n' );
    case 'Pro'
        fprintf('Not yet supported\n' );
        DOF = 0;
    otherwise
        fprintf('Invalid amino acid\n' );
        DOF = 0;
end


%Isolate dipeptide
[allDipeptide,next_pro] = isolate_dipeptide(tempModel2, res_ids, resiId);
if DOF >0 && (size(allDipeptide,1) == numAtom || (size(allDipeptide,1) == numAtom-1 && next_pro ==1))
    
    if next_pro == 0
        [new_Dipeptide, correct_now, ind] = check_Dipeptide_order(allDipeptide(4:size(allDipeptide,1)-3,:), resiName);
        allDipeptide(4:size(allDipeptide,1)-3,:)= new_Dipeptide;
    else
        [new_Dipeptide, correct_now, ind] = check_Dipeptide_order(allDipeptide(4:size(allDipeptide,1)-2,:), resiName);
        allDipeptide(4:size(allDipeptide,1)-2,:)= new_Dipeptide;
    end
    ind0 = ismember(cell2mat(tempModel2(:,1)), cell2mat(allDipeptide(:,1)));
    rest_of_pro = tempModel2(~ind0,:);
    ind0 = ismember(cell2mat(rest_of_pro(:,6)), resiId);
    rest_of_pro = rest_of_pro(~ind0,:);
    
    
    %Only include backbone of the rest of the moving residues
    ind0 = ismember(rest_of_pro(:,3), {'','A'});
    rest_of_pro = rest_of_pro(ind0,:);
    ind2 = ismember(cell2mat(rest_of_pro(:,6)), rest_of_moving);
    ind0a = strcmp(rest_of_pro(:,2), 'CA');
    ind0b = strcmp(rest_of_pro(:,2),'C');
    ind0c  = strcmp(rest_of_pro(:,2),'O');
    ind0d = strcmp(rest_of_pro(:,2),'N');
    ind0e = strcmp(rest_of_pro(:,2),'H');
    ind0f = strcmp(rest_of_pro(:,2),'CB');
    ind0g = strcmp(rest_of_pro(:,2),'HA');
    rest_of_pro = rest_of_pro(~ind2|(ind2 & (ind0a|ind0b|ind0c|ind0d|ind0e|ind0f|ind0g)),:);
    
    
    %% Save dipeptide positions
    XtalPosition = cell2mat(allDipeptide(:,8:10));
    Position = XtalPosition;
    Atom_sizes = cell2mat(allDipeptide(:,12));
    
    %% Limit Rest of Pro to only be atoms within 10A of Cb\
    % basic checks have shown that Phe, Met and Leu always have all
    % atoms stay within 6A of Cb, so add 4 just in case
   %{
    pro_pos = cell2mat(rest_of_pro(:,8:10));
    distp = repmat(Position(8,:),size(pro_pos,1),1)-pro_pos;
    distempP = sqrt(sum(distp.^2,2));
    ind1 = distempP <= 10;
    rest_of_pro = rest_of_pro(ind1,:);
    %}
    if correct_now == 1 && DOF >0
       
        [energy,coordinates, dihedrals] = Rotate_residue_in_protein_energycuttoff(Position, resiName, Atom_sizes,  rest_of_pro, next_pro, energy_cutoff);
        all_coordinates = [coordinates;zeros(1,size(coordinates,2))];
        for i = 1:size(energy)
            all_coordinates(size(all_coordinates,1), (i-1)*3+1)= energy(i);
        end
        
    end
end

end