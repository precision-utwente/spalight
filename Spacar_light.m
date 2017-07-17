function [Results] = Spacar_light(Nodes, Elements, Node_props, Elem_props, Rlse, Optional)
% SPACAR_LIGHT(Nodes, Elements, Node_props, Elem_props, Rlse, Optional)
% runs spacar simulation with reduced set of input arguments. See
% www.spacar.nl for more information.
%
% Created by: M. Nijenhuis and M. Naves
% Contact: m.naves@utwente.nl
%
% CONTRIBUTIONS
% J. Meijaard (complm)
% S.E. Boer (calc_stiffness, calc_inertia, calcTorsStiff and Spavisual functions)
% D. H. Wiersma (CWvalues)
%
% LIMITATIONS
% On each node only a single rotational input can be prescribed (due to concatenation of spatial rotations)
% Output rotations are provided in quaternions and euler rotations in order z, y, x
%
% NOTES
% Torsional stiffness of beams is corrected for constrained warping by
% multiplying torsional stiffness with the constraint warping value


%% INITIALIZE VARIABLES
%#ok<*AGROW>
warning('off','all');               %suppress spacar warnings
Results     = [];
id_inputx   = false;                %identifier to check for prescribed input displacements/rotations
id_inputf   = false;                %identifier to check for external load
x_count     = size(Nodes,1)*2+1;    %counter for node numbering
e_count     = 1;                    %counter for element numbering
X_list      = [];                   %list with node numbers
E_list      = [];                   %list with element numbers


%% CHECK INPUT
%TO BE DONE
%Check input voor transfer functies, en input momenten. Kan worden
%afgemaakt als de omschrijving voor momenten geimplementeerd is.

%CHECK FOR SIMPLE SIMULATION MODE
if (isfield(Optional,'simple_mode') && Optional.simple_mode==1); simple_mode = 1; else simple_mode = 0; end

%CHECK FILENAME
if      (isfield(Optional,'filename') && length(Optional.filename)>19);     fprinf('filename to long, maximum of 20 characters. Filename spacar_file is used instead.');
elseif  (isfield(Optional,'filename') && ~ischar(Optional.filename));       fprinf('filename is not a string. Filename spacar_file is used instead.');
elseif  (isfield(Optional,'filename') && isempty(Optional.filename));       fprinf('filename is empty. Filename spacar_file is used instead.');
else     isfield(Optional,'filename');                                      filename = Optional.filename; end
if      ~exist('filename','var');                                           filename = 'spacar_file'; end

if ~simple_mode %skip checks for simple mode
    %CHECK EXISTENSE OF REQUIRED FUNCTIONS
    if ~exist('spavisual','file')   ==2;   error('spavisual.m not in your path');                               end
    if ~exist('stressbeam','file')  ==2;   error('stressbeam.m not in your path (part of spavisual install)');  end
    if ~exist('spacar','file')      ==3;   error('spacar not in your path');                                    end
    
    %CHECK NODES INPUT VARIABLE
    validateattributes(Nodes,   {'double'},{'ncols',3,'ndims',2},'','Nodes')
    
    %CHECK ELEMENTS INPUT VARIABLE
    validateattributes(Elements,{'double'},{'ncols',2,'ndims',2},'','Elements')
    
    %CHECK NODE_PROPS INPUT VARIABLE
    validateattributes(Node_props,{'struct'},{'nonempty'},'','Node_props')
    for i=1:size(Node_props,2)
        Node_fields = fields(Node_props(i));
        for j=1:length(Node_fields)
            switch Node_fields{j}
                case {'fix_all','fix_xyz','fix_x','fix_y','fix_z','fix_rxyz','fix_rx','fix_ry','fix_rz'}
                    if ~isempty(getfield(Node_props(i),Node_fields{j}));    validateattributes(getfield(Node_props(i),Node_fields{j}),{'logical'},{'scalar'},'',            sprintf('fix property in Node_props(%u)',i));       end %#ok<*GFLD>
                case {'force','force_initial'}
                    if ~isempty(getfield(Node_props(i),Node_fields{j}));     validateattributes(getfield(Node_props(i),Node_fields{j}),{'double'},{'vector','numel',3},'',   sprintf('force property in Node_props(%u)',i));     end
                case {'disp_x','disp_y','disp_z','disp_initial_x','disp_initial_y','disp_initial_z','disp_rx','disp_ry','disp_rz','disp_initial_rx','disp_initial_ry','disp_initial_rz'}
                    if ~isempty(getfield(Node_props(i),Node_fields{j}));     validateattributes(getfield(Node_props(i),Node_fields{j}),{'double'},{'scalar'},'',             sprintf('disp property in Node_props(%u)',i));      end
                case 'mass'
                    if ~isempty(getfield(Node_props(i),Node_fields{j}));     validateattributes(getfield(Node_props(i),Node_fields{j}),{'double'},{'scalar'},'',             sprintf('mass property in Node_props(%u)',i));      end
                case 'inertia'
                    if ~isempty(getfield(Node_props(i),Node_fields{j}));     validateattributes(getfield(Node_props(i),Node_fields{j}),{'double'},{'vector','numel',6},'',   sprintf('inertia property in Node_props(%u)',i));   end
            end
        end
    end
    
    %CHECK ELEM_PROPS INPUT VARIABLE
    validateattributes(Elem_props,{'struct'},{'nonempty'},'','Elem_props')
    for i=1:size(Elem_props,2)
        
        %mandatory fields (%EL_Nrs)
        if ~(isfield(Elem_props(i),'El_Nrs') && ~isempty(Elem_props(i).El_Nrs)); error('Property El_Nrs is not defined in Elem_props(%u)',i);   end
        validateattributes(Elem_props(i).El_Nrs,{'double'},{'vector'},'',sprintf('El_Nrs property in Elem_props(%u)',i));
        
        %mandatory fields when properties are flexible
        if (isfield(Elem_props(i),'flex') && ~isempty(Elem_props(i).flex))
            
            %Check if field exist in structure
            if ~(isfield(Elem_props(i),'flex') && ~isempty(Elem_props(i).flex));    error('Property flex is not defined in Elem_props(%u)',i);  end
            if ~(isfield(Elem_props(i),'E') && ~isempty(Elem_props(i).E));          error('Property E is not defined in Elem_props(%u)',i);     end
            if ~(isfield(Elem_props(i),'G') && ~isempty(Elem_props(i).G));          error('Property G is not defined in Elem_props(%u)',i);     end
            if ~(isfield(Elem_props(i),'rho') && ~isempty(Elem_props(i).rho));      error('Property rho is not defined in Elem_props(%u)',i);   end
            if ~(isfield(Elem_props(i),'type') && ~isempty(Elem_props(i).type));    error('Property type is not defined in Elem_props(%u)',i);  end
            if ~(isfield(Elem_props(i),'dim') && ~isempty(Elem_props(i).dim));      error('Property dim is not defined in Elem_props(%u)',i);   end
            
            %Check if values are valid
            validateattributes(Elem_props(i).flex,{'double'},{'vector'},'',  sprintf('flex property in Elem_props(%u)',i));
            if any(((Elem_props(i).flex==1)+(Elem_props(i).flex==2)+(Elem_props(i).flex==3)+(Elem_props(i).flex==4)+(Elem_props(i).flex==5)+(Elem_props(i).flex==6))==0)
                error('Invalid deformation mode in Elem_props(%u).flex',i)
            end
            validateattributes(Elem_props(i).E,{'double'},{'scalar'},'',                    sprintf('E property in Elem_props(%u)',i));
            validateattributes(Elem_props(i).G,{'double'},{'scalar'},'',                    sprintf('G property in Elem_props(%u)',i));
            validateattributes(Elem_props(i).rho,{'double'},{'scalar'},'',                  sprintf('rho property in Elem_props(%u)',i));
            validateattributes(Elem_props(i).type,{'char'},{'nonempty'},'',                 sprintf('type property in Elem_props(%u)',i));
            validateattributes(Elem_props(i).dim,{'double'},{'vector'},'',   sprintf('dim property in Elem_props(%u)',i));
        end
        
        %Check optional fields
        Elem_fields = fields(Elem_props(i));
        for j=1:length(Elem_fields)
            switch Elem_fields{j}
                case 'orien'
                    if ~isempty(getfield(Elem_props(i),Elem_fields{j}));     validateattributes(getfield(Elem_props(i),Elem_fields{j}),{'double'},{'vector','numel',3},'',             sprintf('orien property in Elem_props(%u)',i));      end
                case 'n_beams'
                    if ~isempty(getfield(Elem_props(i),Elem_fields{j}));     validateattributes(getfield(Elem_props(i),Elem_fields{j}),{'double'},{'scalar'},'',                       sprintf('n_beams property in Elem_props(%u)',i));    end
                case 'color'
                    if ~isempty(getfield(Elem_props(i),Elem_fields{j}));     validateattributes(getfield(Elem_props(i),Elem_fields{j}),{'double'},{'vector','numel',3},'',             sprintf('color property in Elem_props(%u)',i));      end
                case 'hide'
                    if ~isempty(getfield(Elem_props(i),Elem_fields{j}));     validateattributes(getfield(Elem_props(i),Elem_fields{j}),{'logical'},{'scalar'},'',                      sprintf('hide property in Elem_props(%u)',i));      end
            end
        end
    end
    
    %CHECK RLSE INPUT VARIABLE
    if ~isempty(Rlse)
        for i=1:size(Rlse,2)
            if ~isempty(Rlse(i).def)
                validateattributes(Rlse(i).def,{'double'},{'vector'},'',   sprintf('def property in Rlse(%u)',i));
                if any(((Rlse(i).def==1)+(Rlse(i).def==2)+(Rlse(i).def==3)+(Rlse(i).def==4)+(Rlse(i).def==5)+(Rlse(i).def==6))==0)
                    error('Invalid deformation mode in Rlse(%u)',i)
                end
            end
        end
    end
    
    %CHECK OPTIONAL ARGUMENTS
    if (isfield(Optional,'filename') && ~isempty(Optional.filename));
        validateattributes(Optional.filename,{'char'},{'vector'},'',             'filename property in Optional');
    end; if (isfield(Optional,'simple_mode') && ~isempty(Optional.simple_mode));
        validateattributes(Optional.simple_mode,{'logical'},{'scalar'},'',      'simple_mode property in Optional');
    end; if (isfield(Optional,'buck_load') && ~isempty(Optional.buck_load));
        validateattributes(Optional.buck_load,{'logical'},{'scalar'},'',        'buck_load property in Optional');
    end; if (isfield(Optional,'gravity') && ~isempty(Optional.gravity));
        validateattributes(Optional.gravity,{'double'},{'vector','numel',3},'', 'gravity property in Optional');
    end
end



%% START CREATING DATFILE
fileID = fopen([filename '.dat'],'w');

%% USERDEFINED NODES
fprintf(fileID,'#NODES\n');
%print all nodes provides by user
for i=1:size(Nodes,1)
    fprintf(fileID,'X       %3u  %6f  %6f  %6f      #node %u\n',(i-1)*2+1,Nodes(i,1),Nodes(i,2),Nodes(i,3),i);
end


%% ELEMENTS
fprintf(fileID,'\n\n#ELEMENTS\n');

for i=1:size(Elements,1)
    fprintf(fileID,'#ELEMENT %u\n',i);
    
    % get element property set corresponding to element i
    for j=1:size(Elem_props,2)
        if sum(Elem_props(j).El_Nrs==i)>0 %element i is in property set j
            
            %element information
            N           = Elem_props(j).n_beams;    %number of beams
            if isempty(N); N=1; end                 %if number of beams not given, N = 1
            Orien       = Elem_props(j).orien;      %orientation local y-vector
            Flex        = Elem_props(j).flex;       %flexibility of this element
            N_p         = Elements(i,1);            %p-node nodenumber
            N_q         = Elements(i,2);            %q-node nodenumber
            X_list(i,1) = N_p;                      %store p-node in X_list
            
            if N>1 %if more then 1 beam
                
                X_p = Nodes(N_p,1:3);   %Location p-node
                X_q = Nodes(N_q,1:3);   %Location q-node
                V   = X_q - X_p;        %Vector from p to q-node
                
                %create additional intermediate nodes
                for k = 1:N-1
                    X = X_p+V/N*k;              %intermediate node position
                    fprintf(fileID,'X       %3u    %6f  %6f  %6f                        #intermediate node\n',x_count,X(1),X(2),X(3));
                    X_list(i,k+1) = x_count;    %add intermediate node to X_list
                    
                    if k==1 %if the first beam, connect to p-node and first intermediate node
                        fprintf(fileID,'BEAM    %3u  %3u  %3u  %3u  %3u  %6f  %6f  %6f      #beam %u\n',e_count,(N_p-1)*2+1,(N_p-1)*2+2,x_count,x_count+1,Orien(1),Orien(2),Orien(3),k);
                    else    %if not the first beam, connect to two intermediate nodes
                        fprintf(fileID,'BEAM    %3u  %3u  %3u  %3u  %3u  %6f  %6f  %6f      #beam %u\n',e_count,x_count-2,x_count-1,x_count,x_count+1,Orien(1),Orien(2),Orien(3),k);
                    end
                    
                    if ~isempty(Flex)        %if element has flexibility, add dyne (no rlse, rlse is only added to last beam in element i)
                        fprintf(fileID,'dyne    %3u',e_count);
                        for m=1:length(Flex) %loop over all flexible deformation modes
                            fprintf(fileID,'  %3u',Flex(m));
                        end
                        fprintf(fileID,'\n');
                    end
                    
                    E_list(i,k) = e_count;      %add beam number to E_list
                    e_count     = e_count+1;    %increase beam counter by 1
                    x_count     = x_count+2;    %increase node counter by 2 (+1 for rotation node)
                end
                
                %for the last beam in element i, connect to last intermediate node and q-node
                fprintf(fileID,'BEAM    %3u  %3u  %3u  %3u  %3u  %6f  %6f  %6f      #beam %u\n',e_count,x_count-2,x_count-1,(N_q-1)*2+1,(N_q-1)*2+2,Orien(1),Orien(2),Orien(3),k+1);
                
                X_list(i,k+2) = N_q;        %add q-node to X_list
                E_list(i,k+1) = e_count;    %add beam number to E_list
                
            else %if only a single beam is used, directly connect to p and q-node without intermediate noodes
                fprintf(fileID,'BEAM    %3u  %3u  %3u  %3u  %3u  %6f  %6f  %6f      #beam\n',e_count,(N_p-1)*2+1,(N_p-1)*2+2,(N_q-1)*2+1,(N_q-1)*2+2,Orien(1),Orien(2),Orien(3));
                
                X_list(i,2) = N_q;          %add q-node to X_list
                E_list(i,1) = e_count;      %add beam number to E_list
            end
            
            %for the last beam only, add dyne and/or rlse
            if isempty(Rlse)            %if no rlse, add all flexible deformation modes as dyne
                fprintf(fileID,'dyne    %3u',e_count);
                for m=1:length(Flex)    %loop over all flexible deformation modes
                    fprintf(fileID,'  %3u',Flex(m));
                end
                fprintf(fileID,'\n');
                
            else%if some rlse are specified
                %compensate size of Rlse if size is smaller then element list
                if i>size(Rlse,2)
                    Rlse(i).def = [];
                end
                
                % add dyne
                if ~isempty(Flex)                           %if some flexibility is specified
                    dyn_added = false;                      %reset identifier to check if string 'dyne' is added
                    for m=1:length(Flex)                    %loop over all flexible deformation modes
                        if ~(sum(Rlse(i).def==Flex(m))>0)   %if flexible deformation mode is not a rlse, it is dyne
                            if ~dyn_added                   %only add string 'dyne' if it is not yet added
                                fprintf(fileID,'dyne    %3u',e_count);
                                dyn_added = true;           %set 'dyne' identifier
                            end
                            fprintf(fileID,'  %3u',Flex(m));
                        end
                    end
                    fprintf(fileID,'\n');
                end
                
                % add rlse
                rlse_added = false;                     %reset identifier to check if string 'rlse' is added
                for m=1:length(Rlse(i).def)             %loop over all released deformation modes
                    if ~rlse_added                      %only add string 'rlse' if it is not yet added
                        fprintf(fileID,'rlse    %3u',e_count);
                        rlse_added = true;
                    end
                    fprintf(fileID,'  %3u',Rlse(i).def(m));
                end
                fprintf(fileID,'\n');
            end
            e_count = e_count+1; %increase beam counter by 1 for last beam in the element
        end
    end
    fprintf(fileID,'\n');
end


%% NODE FIXES AND INPUTS
fprintf(fileID,'\n\n#NODE FIXES AND INPUTS\n');
for i=1:size(Node_props,2)
    %fixes
    if(isfield(Node_props(i),'fix_all') && ~isempty(Node_props(i).fix_all));    fprintf(fileID,'FIX      %u  \n',(i-1)*2+1);
        fprintf(fileID,'FIX      %u  \n',(i-1)*2+2);   end
    if(isfield(Node_props(i),'fix_xyz') && ~isempty(Node_props(i).fix_xyz));    fprintf(fileID,'FIX      %u  \n',(i-1)*2+1);   end
    if(isfield(Node_props(i),'fix_x') && ~isempty(Node_props(i).fix_x));        fprintf(fileID,'FIX      %u 1 \n',(i-1)*2+1);  end
    if(isfield(Node_props(i),'fix_y') && ~isempty(Node_props(i).fix_y));        fprintf(fileID,'FIX      %u 2 \n',(i-1)*2+1);  end
    if(isfield(Node_props(i),'fix_z') && ~isempty(Node_props(i).fix_z));        fprintf(fileID,'FIX      %u 3 \n',(i-1)*2+1);  end
    if(isfield(Node_props(i),'fix_rxyz') && ~isempty(Node_props(i).fix_rxyz));  fprintf(fileID,'FIX      %u  \n',(i-1)*2+2);   end
    if(isfield(Node_props(i),'fix_rx') && ~isempty(Node_props(i).fix_rx));      fprintf(fileID,'FIX      %u 4  \n',(i-1)*2+2); end
    if(isfield(Node_props(i),'fix_ry') && ~isempty(Node_props(i).fix_ry));      fprintf(fileID,'FIX      %u 3  \n',(i-1)*2+2); end
    if(isfield(Node_props(i),'fix_rz') && ~isempty(Node_props(i).fix_rz));      fprintf(fileID,'FIX      %u 2  \n',(i-1)*2+2); end
    
    %input displacements
    if((isfield(Node_props(i),'disp_x') && ~isempty(Node_props(i).disp_x)) ||...
            (isfield(Node_props(i),'disp_initial_x') && ~isempty(Node_props(i).disp_initial_x)));   fprintf(fileID,'INPUTX      %3u     1\n',(i-1)*2+1);id_inputx = true;    end
    if((isfield(Node_props(i),'disp_y') && ~isempty(Node_props(i).disp_y)) ||...
            (isfield(Node_props(i),'disp_initial_y') && ~isempty(Node_props(i).disp_initial_y)));   fprintf(fileID,'INPUTX      %3u     2\n',(i-1)*2+1);id_inputx = true;    end
    if((isfield(Node_props(i),'disp_z') && ~isempty(Node_props(i).disp_z)) ||...
            (isfield(Node_props(i),'disp_initial_z') && ~isempty(Node_props(i).disp_initial_z)));   fprintf(fileID,'INPUTX      %3u     3\n',(i-1)*2+1);id_inputx = true;    end
    
    %input rotations
    id_inputr = 0; %identifier to count the number of input rotations
    if((isfield(Node_props(i),'disp_rx') && ~isempty(Node_props(i).disp_rx)) ||...
            (isfield(Node_props(i),'disp_initial_rx') && ~isempty(Node_props(i).disp_initial_rx))); fprintf(fileID,'INPUTX      %3u     4\n',(i-1)*2+2);id_inputx = true; id_inputr=id_inputr+1; end
    if((isfield(Node_props(i),'disp_ry') && ~isempty(Node_props(i).disp_ry)) ||...
            (isfield(Node_props(i),'disp_initial_ry') && ~isempty(Node_props(i).disp_initial_ry))); fprintf(fileID,'INPUTX      %3u     3\n',(i-1)*2+2);id_inputx = true; id_inputr=id_inputr+1; end
    if((isfield(Node_props(i),'disp_rz') && ~isempty(Node_props(i).disp_rz)) ||...
            (isfield(Node_props(i),'disp_initial_rz') && ~isempty(Node_props(i).disp_initial_rz))); fprintf(fileID,'INPUTX      %3u     2\n',(i-1)*2+2);id_inputx = true; id_inputr=id_inputr+1; end
    if id_inputr>1 %if multiple rotations are prescribed, problems can arise with quaternion<->euler conversion
        error('Only a single input rotation can be added on a single node. Multiple rotational inputs defined for node %u',i)
    end
end
fprintf(fileID,'\n\nEND\nHALT\n\n');


%% STIFFNESS/INERTIA PROPS
for i=1:size(Elem_props,2) %loop over each element property set
    
    %only write stiffness values when deformations are flexible
    if (isfield(Elem_props(i),'flex') && ~isempty(Elem_props(i).flex))
        stiffness = calc_stiffness(Elem_props(i)); %calculate stiffness values
        
        for j=1:length(Elem_props(i).El_Nrs)                        %loop over all elemenents in element property set
            L   = norm(Nodes(Elements(Elem_props(i).El_Nrs(j),2),:)...
                - Nodes(Elements(Elem_props(i).El_Nrs(j),1),:));    %calculate flexure length for constrained warping values
            cw  = CWvalues(L,Elem_props(i));                        %calculate constrained warping values
            for k=1:size(E_list,2)                                  %loop over all beams in the element
                El = E_list(Elem_props(i).El_Nrs(j),k);
                if El>0
                    fprintf(fileID,'ESTIFF      %u %f %f %f %f %f %f \n',El,stiffness(1),cw*stiffness(2),stiffness(3),stiffness(4),stiffness(5),stiffness(6));
                end
            end
        end
    end
    
    %write mass/inertia values
    inertia = calc_inertia(Elem_props(i));          %calculate mass properties
    for j=1:length(Elem_props(i).El_Nrs)            %loop over all element property sets
        for k=1:size(E_list,2)
            El = E_list(Elem_props(i).El_Nrs(j),k); %loop over all beams in element set
            if El>0
                fprintf(fileID,'EM      %u %f %f %f %f %f \n',El,inertia(1),inertia(2),inertia(3),inertia(4),inertia(5));
            end
        end
    end
end


%% FORCES/MOMENTS/NODAL MASSES
fprintf(fileID,'\n\n#FROCES/MOMENTS\n');
id_ini  = false; %check for initial loading or displacement
id_add  = false; %check for aditional loading or displacement

for i=1:size(Node_props,2) %loop over all user defined nodes
    %forces
    if(isfield(Node_props(i),'force') && ~isempty(Node_props(i).force));                    fprintf(fileID,'DELXF   %3u %6f %6f %6f  \n',(i-1)*2+1,Node_props(i).force(1),Node_props(i).force(2),Node_props(i).force(3));                           id_add = true;  id_inputf=true; end
    if(isfield(Node_props(i),'force_initial') && ~isempty(Node_props(i).force_initial));    fprintf(fileID,'XF      %3u %6f %6f %6f  \n',(i-1)*2+1,Node_props(i).force_initial(1),Node_props(i).force_initial(2),Node_props(i).force_initial(3));   id_ini = true;  id_inputf=true; end
    
    %TO BE DONE
    %if(isfield(Node_props(i),'moment') && ~isempty(Node_props(i).moment));
    %   moments = e2q(Node_props(i).moment);
    %                                                                                        fprintf(fileID,'DELXF   %3u %6f %6f %6f %6f  \n',(i-1)*2+2,moments(1),moments(2),moments(3),moments(4));                                                id_add = true;  id_inputf=true; end
    %if(isfield(Node_props(i),'moment_initial') && ~isempty(Node_props(i).moment_initial));
    %   moments = e2q(Node_props(i).moment_initial) ;
    %                                                                                        fprintf(fileID,'XF      %3u %6f %6f %6f %6f  \n',(i-1)*2+2,moments(1),moments(2),moments(3),moments(4));                                                id_ini = true;  id_inputf=true; end
    
    %displacements
    if(isfield(Node_props(i),'disp_x') && ~isempty(Node_props(i).disp_x));                  fprintf(fileID,'DELINPX  %3u  1  %6f  \n',(i-1)*2+1,Node_props(i).disp_x(1));                       id_add = true; end
    if(isfield(Node_props(i),'disp_y') && ~isempty(Node_props(i).disp_y));                  fprintf(fileID,'DELINPX  %3u  2  %6f  \n',(i-1)*2+1,Node_props(i).disp_y(1));                       id_add = true; end
    if(isfield(Node_props(i),'disp_z') && ~isempty(Node_props(i).disp_z));                  fprintf(fileID,'DELINPX  %3u  3  %6f  \n',(i-1)*2+1,Node_props(i).disp_z(1));                       id_add = true; end
    if(isfield(Node_props(i),'disp_initial_x') && ~isempty(Node_props(i).disp_initial_x));  fprintf(fileID,'INPUTX   %3u  1  %6f  \n',(i-1)*2+1,Nodes(i,1) + Node_props(i).disp_initial_x(1));  id_ini = true; end
    if(isfield(Node_props(i),'disp_initial_y') && ~isempty(Node_props(i).disp_initial_y));  fprintf(fileID,'INPUTX   %3u  2  %6f  \n',(i-1)*2+1,Nodes(i,2) + Node_props(i).disp_initial_y(1));  id_ini = true; end
    if(isfield(Node_props(i),'disp_initial_z') && ~isempty(Node_props(i).disp_initial_z));  fprintf(fileID,'INPUTX   %3u  3  %6f  \n',(i-1)*2+1,Nodes(i,3) +Node_props(i).disp_initial_z(1));   id_ini = true; end
    
    %rotations
    if(isfield(Node_props(i),'disp_rx') && ~isempty(Node_props(i).disp_rx));                rot = e2q([Node_props(i).disp_rx(1) 0 0]);
        fprintf(fileID,'DELINPX     %3u     4   %6f  \n',(i-1)*2+2,rot(4)); id_add = true; end
    if(isfield(Node_props(i),'disp_ry') && ~isempty(Node_props(i).disp_ry));                rot = e2q([0 Node_props(i).disp_ry(1) 0]);
        fprintf(fileID,'DELINPX     %3u     3   %6f  \n',(i-1)*2+2,rot(3)); id_add = true; end
    if(isfield(Node_props(i),'disp_rz') && ~isempty(Node_props(i).disp_rz));                rot = e2q([0 0 Node_props(i).disp_rz(1)]);
        fprintf(fileID,'DELINPX     %3u     2   %6f  \n',(i-1)*2+2,rot(2)); id_add = true; end
    if(isfield(Node_props(i),'disp_initial_rx') && ~isempty(Node_props(i).disp_initial_rx));rot = e2q([Node_props(i).disp_initial_rx(1) 0 0]);
        fprintf(fileID,'INPUTX     %3u     4   %6f  \n',(i-1)*2+2,rot(4));  id_ini = true; end
    if(isfield(Node_props(i),'disp_initial_ry') && ~isempty(Node_props(i).disp_initial_ry));rot = e2q([0 Node_props(i).disp_initial_ry(1) 0]);
        fprintf(fileID,'INPUTX     %3u     3   %6f  \n',(i-1)*2+2,rot(3));  id_ini = true; end
    if(isfield(Node_props(i),'disp_initial_rz') && ~isempty(Node_props(i).disp_initial_rz));rot = e2q([0 0 Node_props(i).disp_initial_rz(1)]);
        fprintf(fileID,'INPUTX     %3u     2   %6f  \n',(i-1)*2+2,rot(2));  id_ini = true; end
    
    %nodal masses/inertia
    if(isfield(Node_props(i),'mass') && ~isempty(Node_props(i).mass));                      fprintf(fileID,'XM      %3u %6f  \n',(i-1)*2+1,Node_props(i).mass); end
    if(isfield(Node_props(i),'inertia') && ~isempty(Node_props(i).inertia));                fprintf(fileID,'XM      %3u %6f %6f %6f %6f %6f %6f\n',(i-1)*2+2,Node_props(i).inertia(1),Node_props(i).inertia(2),Node_props(i).inertia(3),...
            Node_props(i).inertia(4),Node_props(i).inertia(5),Node_props(i).inertia(6)); end
end


%% ADITIONAL OPTIONS
%GRAVITY
if(isfield(Optional,'gravity') && ~isempty(Optional.gravity)); fprintf(fileID,'\nGRAVITY  %6f %6f %6f',Optional.gravity(1),Optional.gravity(2),Optional.gravity(3)); end

%ITERSTEP SETTINGS
if      (id_ini && id_add);      fprintf(fileID,'\nITERSTEP 10 10 0.0000005 1 3 10');    %if initial and aditional loading/displacement
elseif  (id_ini && ~id_add);     fprintf(fileID,'\nITERSTEP 10 1 0.0000005  1 1 10');    %if initial loading/displacement
elseif  (~id_ini && id_add);     fprintf(fileID,'\nITERSTEP 10 10 0.0000005 1 3 0');     %if initial loading/displacement
else                             fprintf(fileID,'\nITERSTEP 10 1  0.0000005 1 1 0');  end%no loading/displacement

%TRANSFER FUNCTION INPUT/OUTPUT
if ((isfield(Optional,'transfer_in') && ~isempty(Optional.transfer_in)) ||  (isfield(Optional,'transfer_out') && ~isempty(Optional.transferout)))
    if id_inputx
        disp('Warning: input displacement is prediscribed, possibly affecting input/ouput transferfunction')
    end
    
    fprintf(fileID,'\n\nEND\nHALT\n\n');
    for i=1:size(Optional.transfer_in,2) %add inputs
        switch Optional.transfer_in(i).type
            case 'force_x';     fprintf(fileID,'\nINPUTF %2u %3u 1',i,(Optional.transfer_in(i).node-1)*2+1);
            case 'force_y';     fprintf(fileID,'\nINPUTF %2u %3u 2',i,(Optional.transfer_in(i).node-1)*2+1);
            case 'force_z';     fprintf(fileID,'\nINPUTF %2u %3u 3',i,(Optional.transfer_in(i).node-1)*2+1);
                %TO BE DONE
                %case 'moment_x';    fprintf(fileID,'\nINPUTF %2u %3u 4',i,(Optional.transfer_in(i).node-1)*2+2);
                %case 'moment_y';    fprintf(fileID,'\nINPUTF %2u %3u 3',i,(Optional.transfer_in(i).node-1)*2+2);
                %case 'moment_z';    fprintf(fileID,'\nINPUTF %2u %3u 2',i,(Optional.transfer_in(i).node-1)*2+2);
                
            case 'disp_x';      fprintf(fileID,'\nINX %2u %3u 1',i,(Optional.transfer_in(i).node-1)*2+1);
            case 'disp_y';      fprintf(fileID,'\nINX %2u %3u 2',i,(Optional.transfer_in(i).node-1)*2+1);
            case 'disp_z';      fprintf(fileID,'\nINX %2u %3u 3',i,(Optional.transfer_in(i).node-1)*2+1);
                %case 'rot_x';       fprintf(fileID,'\nINX %2u %3u 4',i,(Optional.transfer_in(i).node-1)*2+2);
                %case 'rot_y';       fprintf(fileID,'\nINX %2u %3u 3',i,(Optional.transfer_in(i).node-1)*2+2);
                %case 'rot_z';       fprintf(fileID,'\nINX %2u %3u 2',i,(Optional.transfer_in(i).node-1)*2+2);
        end
    end
    for i=1:size(Optional.transfer_out,2) %add outputs
        switch Optional.transfer_out(i).type
            case 'force_x';     fprintf(fileID,'\nOUTF %2u %3u 1',i,(Optional.transfer_out(i).node-1)*2+1);
            case 'force_y';     fprintf(fileID,'\nOUTF %2u %3u 2',i,(Optional.transfer_out(i).node-1)*2+1);
            case 'force_z';     fprintf(fileID,'\nOUTF %2u %3u 3',i,(Optional.transfer_out(i).node-1)*2+1);
                %TO BE DONE
                %case 'moment_x';    fprintf(fileID,'\nOUTF %2u %3u 4',i,(Optional.transfer_out(i).node-1)*2+2);
                %case 'moment_y';    fprintf(fileID,'\nOUTF %2u %3u 3',i,(Optional.transfer_out(i).node-1)*2+2);
                %case 'moment_z';    fprintf(fileID,'\nOUTF %2u %3u 2',i,(Optional.transfer_out(i).node-1)*2+2);
                
            case 'disp_x';      fprintf(fileID,'\nOUTX %2u %3u 1',i,(Optional.transfer_out(i).node-1)*2+1);
            case 'disp_y';      fprintf(fileID,'\nOUTX %2u %3u 2',i,(Optional.transfer_out(i).node-1)*2+1);
            case 'disp_z';      fprintf(fileID,'\nOUTX %2u %3u 3',i,(Optional.transfer_out(i).node-1)*2+1);
                %case 'rot_x';       fprintf(fileID,'\nOUTX %2u %3u 4',i,(Optional.transfer_out(i).node-1)*2+2);
                %case 'rot_y';       fprintf(fileID,'\nOUTX %2u %3u 3',i,(Optional.transfer_out(i).node-1)*2+2);
                %case 'rot_z';       fprintf(fileID,'\nOUTX %2u %3u 2',i,(Optional.transfer_out(i).node-1)*2+2);
        end
    end
end

fprintf(fileID,'\n\nEND\nEND\n\n');


%% VISUALIZATION
fprintf(fileID,'\n\nVISUALIZATION');
for i=1:size(Elem_props,2) %loop over all element property sets
    
    %CROSSECTIONAL PROPERTIES
    fprintf(fileID,'\n\nBEAMPROPS ');
    for j=1:length(Elem_props(i).El_Nrs) %loop over all elemenents in element set i
        for k=1:size(E_list,2)
            El = E_list(Elem_props(i).El_Nrs(j),k);
            if El>0
                fprintf(fileID,' %u ',El);
            end
        end
    end
    switch Elem_props(i).type
        case {'leafspring','rigid'} %if leafspring or rigid, rect crossection
            fprintf(fileID,'\nCROSSTYPE  RECT');
            fprintf(fileID,'\nCROSSDIM  %f  %f',Elem_props(i).dim(1),Elem_props(i).dim(2));
        case 'wire'                 %if wire, circular crossection
            fprintf(fileID,'\nCROSSTYPE  CIRC');
            fprintf(fileID,'\nCROSSDIM  %f ',Elem_props(i).dim(1));
    end
    
    %COLOR
    if (isfield(Elem_props(i),'color') && ~isempty(Elem_props(i).color))
        fprintf(fileID,'\nGRAPHICS ');
        for j=1:length(Elem_props(i).El_Nrs)
            for k=1:size(E_list,2)
                El = E_list(Elem_props(i).El_Nrs(j),k);
                if El>0
                    fprintf(fileID,' %u ',El);
                end
            end
        end
        fprintf(fileID,'\nFACECOLOR  %f %f %f ',Elem_props(i).color(1),Elem_props(i).color(2),Elem_props(i).color(3));
    end
    
    %VISIBILITY
    if (isfield(Elem_props(i),'hide') && ~isempty(Elem_props(i).hide) && Elem_props(i).hide==1)
        fprintf(fileID,'\nDONOTDRAW ');
        for j=1:length(Elem_props(i).El_Nrs)
            for k=1:size(E_list,2)
                El = E_list(Elem_props(i).El_Nrs(j),k);
                if El>0
                    fprintf(fileID,' %u ',El);
                end
            end
        end
    end
end
fclose(fileID); %datfile finished!


%% SIMULATE CONSTRAINTS
if simple_mode==0
    spacar(0,filename)
    
    %CHECK CONSTRAINTS
    sbd     = [filename '.sbd'];
    nep     = getfrsbf(sbd,'nep');
    nxp     = getfrsbf(sbd,'nxp');
    le      = getfrsbf(sbd,'le');
    BigD    = getfrsbf(sbd,'bigd',1);
    Dcc     = BigD( 1:(nep(1)+nep(3)+nep(4)) , nxp(1)+(1:nxp(2)) );
    [ U, s, V ] = svd(Dcc);
    s       = diag(s);
    
    nsing = length(find(s<sqrt(eps)*s(1))); %number of near zero singular values
    if length(U)>length(V)
        nover  = length(U)-length(V)+nsing;
        nunder = nsing;
    elseif length(U)<length(V)
        nover  = nsing;
        nunder = length(V)-length(U)+nsing;
    else
        nover  = nsing;
        nunder = nsing;
    end
    
    if nunder>0 %underconstrained
        fprintf('\nSystem is underconstrained. Check element conectivity and fixes\n')
        warning('on','all')
        return
    elseif nover>0 %overconstrained
        overconstraint = U(:,end-nover+1:end);
        oc = overconstraint(:,1);
        [oc_sort,order] = sort(oc.^2,1,'descend');
        idx = find(cumsum(oc_sort)>sqrt(0.95),1,'first');% select only part that explains 95% (or more) of singular vector's length
        idx = order(1:idx);
        sel = (1:numel(oc))';
        sel = sel(idx);
        
        listData = zeros(numel(sel),3);
        listData(:,3) = oc(idx);
        for i=1:size(listData,1)
            [elnr,defpar] = find(le==sel(i));
            listData(i,1:2) = [elnr, defpar]; %put overconstrained element numbers and deformations in listData
        end
        
        %Reshape rlse suggestions according to user defined elements
        OC_el= [];
        OC_defs = [];
        for i=1:size(E_list,1)
            list = [];
            for j=1:size(E_list,2)
                list = [list; sort(listData(find(listData(:,1)==E_list(i,j)),2))]; %#ok<FNDSB>
            end
            red_list=[];
            for j=1:6
                if sum(list==j)>0
                    red_list(end+1) = j;
                end
            end
            if ~isempty(red_list)
                OC_defs(end+1,1:length(red_list)) = red_list;
                OC_el(end+1,1) = i;
            end
        end
        
        Results.overconstraints = [OC_el OC_defs];
        fprintf('\nSystem is overconstrained, releases are required in order to run static simulation.\nA suggestion for possible releases is given in Results.overconstraints in the workspace and the table below.\n')
        fprintf('\nNumber of overconstraints: %u\n\n',nover);
        disp(table(OC_el,sum((OC_defs==1),2),sum((OC_defs==2),2),sum((OC_defs==3),2),sum((OC_defs==4),2),sum((OC_defs==5),2),sum((OC_defs==6),2),...
            'VariableNames',{'Element' 'def_1' 'def_2 ' 'def_3' 'def_4' 'def_5' 'def_6'}));
        warning('on','all')
        return
    end
end


%% SIMULATE STATICS
try
    %TO BE DONE
    %of toch mode 10 draaien, mode 9 werkt zeer slecht voor input
    %displacements/rotaties. Mode 10 geeft niet de transfer functies
    spacar(9,filename)
    if simple_mode==0
        spavisual(filename)
    end
    %get results
    Results = calc_Results(filename, E_list, id_inputf, id_inputx, Nodes, Elements, Node_props, Elem_props, Rlse, Optional);
catch
    disp('Spacar simulation failed. Possibly failed to converge to solution. Check magnitude of input displacements, loads and other input data.')
end

%% END OF SPACAR_LIGHT
end






%% AUXILIARY FUNCTIONS
function stiffness = calc_stiffness(Elem_props)
% Compute the stiffness properties for leafspring or wireflexure
type    = Elem_props.type;
dim     = Elem_props.dim;
E       = Elem_props.E;
G       = Elem_props.G;
v       = E/(2*G) - 1;
switch lower(type)
    case 'leafspring'
        t   = dim(1);
        w   = dim(2);
        A   = t*w;
        It 	= calcTorsStiff(t,w);
        Iy  = (1/12)*t*w^3;
        Iz  = (1/12)*w*t^3;
        k   = 10*(1+v)/(12+11*v);
    case 'wire'
        d   = dim(1);
        A   = (pi/4)*d^2;
        It  = (pi/32)*d^4;
        Iy  = (pi/64)*d^4;
        Iz  = Iy;
        k   = 6*(1+v)/(7+6*v);
end
stiffness(1,1) = E*A;
stiffness(1,2) = G*It;
stiffness(1,3) = E*Iy;
stiffness(1,4) = E*Iz;
stiffness(1,5) = stiffness(1,3)/(G*A*k);
stiffness(1,6) = stiffness(1,4)/(G*A*k);
end


function Ip = calcTorsStiff(t, w)
% Compute the polar moment of inertia for rectangular cross-sections
if w > t
    a = t/2;
    b = w/2;
else
    a = w/2;
    b = t/2;
end
sumN_new    = 0;
sumN        = 0;
n           = 1;
while (n < 3 || abs((sumN_new/sumN)) > 0.0001)
    sumN_new    = (1/n^5) * tanh(n*pi*b/(2*a));
    n           = n + 2;
    sumN        = sumN + sumN_new;
end
Ip = 1/3 * (2*a)^3*(2*b) * (1 - (192/pi^5)*(a/b)*sumN);
end



function inertia = calc_inertia(Elem_props)
% Compute the inertia properties for leafspring or wireflexure
type    = Elem_props.type;
dim     = Elem_props.dim;
rho     = Elem_props.rho;
switch lower(type)
    case {'leafspring','rigid'}
        t   = dim(1);
        w   = dim(2);
        A   = t*w;
        Iy  = 1/12 * t*w^3;
        Iz  = 1/12 * w*t^3;
    case 'wire'
        d   = dim(1);
        A   = (pi/4)*d^2;
        Iy  = (pi/64)*d^4;
        Iz  = Iy;
end
inertia(1,1) = rho*A;
inertia(1,2) = rho*(Iy+Iz);
inertia(1,3) = rho*Iy;
inertia(1,4) = rho*Iz;
inertia(1,5) = 0;
end


function Results = calc_Results(filename, E_list, id_inputf, id_inputx, Nodes, ~, ~, Elem_props, ~, Optional)
nddof   = getfrsbf([filename '.sbd'],'nddof'); %number of dynamic DOFs
t_list  =  1:getfrsbf([filename,'.sbd'],'tdef'); %timesteps
lnp     = getfrsbf([filename,'.sbd'],'lnp'); %lnp data
%CHECK BUCKLING SETTINGS
calcbuck = false;
if (isfield(Optional,'buck_load') && Optional.buck_load == 1)
    calcbuck = true;
    if id_inputx
        disp('Warning: input displacement prescribed, buckling load multipliers are also with respect to input reaction forces');
    end
    if ~id_inputf
        disp('Warning: No external forces are prescribed. Buckling values are not calculated');
        calcbuck = false;
    end
end

%EIGENFREQUENCIES and BUCKLING
for i=1:length(t_list)
    M = getfrsbf([filename '.sbm'] ,'m0', t_list(i));
    G = getfrsbf([filename '.sbm'] ,'g0', i);
    K = getfrsbf([filename '.sbm'] ,'k0', t_list(i)) + getfrsbf([filename '.sbm'] ,'n0', t_list(i)) + G;
    %C = getfrsbf([filename '.sbm'] ,'c0', t_list(i)) + getfrsbf([filename '.sbm'] ,'d0', t_list(i));
    
    [~,D]   = eig(K(1:nddof,1:nddof),M(1:nddof,1:nddof));
    D       = diag(D);
    [~,o]   = sort(abs(D(:)));
    d       = D(o);
    Results.step(i).Freq = sqrt(d)*1/(2*pi);
    
    if calcbuck
        [~,Buck] = eig(-K,G);
        Results.step(i).Buck = sort(abs(diag(Buck)));
    end
end

[propcrossect, Sig_nums]  = calc_propcrossect(E_list,Elem_props);
for i=t_list
    x       = getfrsbf([filename '.sbd'] ,'x', i);
    fxtot   = getfrsbf([filename '.sbd'] ,'fxt',i);
    for j=1:size(Nodes,1)
        Results.step(i).node(j).x           = x(lnp((j-1)*2+1,1:3));
        Results.step(i).node(j).rx_eulzyx   = q2e(x(lnp((j-1)*2+2,1:4))');
        Results.step(i).node(j).rx_quat     = (x(lnp((j-1)*2+2,1:4))');
        Results.step(i).node(j).Freac       = fxtot(lnp((j-1)*2+1,1:3)) ;
        %TO BE DONE
        %Results.step(i).node(j).Mreac       = q2e(fxtot(lnp((j-1)*2+2,1:4))');
        [Results.step(i).node(j).CMglob, Results.step(i).node(j).CMloc]  =  complm(filename,(j-1)*2+1,(j-1)*2+2,i); %#ok<*AGROW>
    end
    [~,~,~,stressextrema] = stressbeam([filename,'.sbd'],Sig_nums,i,[],propcrossect);
    Results.step(i).stressmax = stressextrema.max*1e6;
    %  Results.step(i).bode_data =  getss('spacarfile',i);
end
Results.ndof = getfrsbf([filename '.sbd'] ,'ndof');

end


function cw= CWvalues(L,Elem_props)
w = Elem_props.dim(2);
E = Elem_props.E;
G = Elem_props.G;
v = E/(2*G) - 1;

aspect  = L/w;        %- aspect ratio of flexure
c       = sqrt(24/(1+v));
cw      = (aspect*c/(aspect*c-2*tanh(aspect*c/2)));
end


function [CMglob, CMloc] = complm(filename,ntr,nrot,tstp)
% Calculate the compliance matrix in global directions and in body-fixed
% local directions.

tstp    = tstp-1;
CMglob  =zeros(6,6);
CMloc   =zeros(6,6);
if (nargin < 3) || (nargin>4)
    disp('complm needs 3 or 4 input arguments');
    return;
end;
if nargin < 4, tstp=0; end;
% get data from files
DX      =getfrsbf([filename '.sbd'],'dx',tstp+1);
X       =getfrsbf([filename '.sbd'],'x',tstp+1);
lnp     =getfrsbf([filename '.sbd'],'lnp');
nxp     =getfrsbf([filename '.sbd'],'nxp');
nep     =getfrsbf([filename '.sbd'],'nep');
% nddof   =getfrsbf([filename '.sbd'],'nddof');
% locate the place of the coordinates of the points where the compliance
% has to be determined
locv    =[lnp(ntr,1:3), lnp(nrot,1:4)];
% test whether the selected coordinates are feasible
for i=1:7
    if locv(i) <= 0
        disp('ERROR: invalid node number');
        return;
    end;
    if locv(i) <= nxp(1) || ...
            (locv(i)>(nxp(1)+nxp(2)) && locv(i) <= (nxp(1)+nxp(2)+nxp(3)))
        %  disp('WARNING: constrained node');
    end;
end;
% search for the right degrees of freedom
locdof=[nxp(3)+(1:nxp(4)), nxp(3)+nxp(4)+nep(3)+(1:nep(4))];

% reduce DX to the rows needed
DX=DX(locv,locdof);
K0      =getfrsbf([filename '.sbm'],'k0',tstp+1);
G0      =getfrsbf([filename '.sbm'],'g0',tstp+1);
CMlambda=DX*((K0+G0)\(DX'));
% Reduce CMlambda to the correct matrices by the lambda matrices
lambda0=X(locv(4));
lambda1=X(locv(5));
lambda2=X(locv(6));
lambda3=X(locv(7));
lambdabt=[ -lambda1  lambda0 -lambda3  lambda2
    -lambda2  lambda3  lambda0 -lambda1
    -lambda3 -lambda2  lambda1  lambda0 ];
lambdat= [ -lambda1  lambda0  lambda3 -lambda2
    -lambda2 -lambda3  lambda0  lambda1
    -lambda3  lambda2 -lambda1  lambda0 ];
Tglob= [ eye(3)     zeros(3,4)
    zeros(3,3) 2*lambdabt];
Tloc = [ lambdat*(lambdabt') zeros(3,4)
    zeros(3,3) 2*lambdat];
CMglob=Tglob*CMlambda*(Tglob');
CMloc=Tloc*CMlambda*(Tloc');
end


function [propcrossect, Sig_nums]  = calc_propcrossect(E_list,Elem_props)
%restructure crossectional properties to evaluate stresses throuqh
%stressbeam.m

propcrossect = [];Sig_nums = [];

for i=1:size(E_list,1)
    id=[];
    for j=1:size(Elem_props,2)
        for k=1:length(Elem_props(j).El_Nrs)
            if Elem_props(j).El_Nrs(k)==i;
                id=j;
            end
        end
    end
    if ~isempty(id)
        if (isfield(Elem_props(id),'flex') && ~isempty(Elem_props(id).flex))
            Elements = E_list(i,:);
            Elements(Elements==0) = [];
            Sig_nums = [Sig_nums Elements];
            
            switch Elem_props(id).type
                case 'leafspring'
                    for j=1:length(Elements)
                        propcrossect(end+1).CrossSection = 'rect';
                        propcrossect(end).Dimensions = [Elem_props(id).dim(1),Elem_props(id).dim(2)];
                    end
                case 'wire'
                    for j=1:length(Elements)
                        propcrossect(end+1).CrossSection = 'circ';
                        propcrossect(end).Dimensions = Elem_props(id).dim(1);
                    end
            end
        end
    end
end
end


function q = e2q( eul )
%e2q Convert Euler angles to quaternion

% Pre-allocate output
q = zeros(size(eul,1), 4, 'like', eul);

% Compute sines and cosines of half angles
c = cos(eul/2);
s = sin(eul/2);


q = [c(:,1).*c(:,2).*c(:,3)+s(:,1).*s(:,2).*s(:,3), ...
    c(:,1).*c(:,2).*s(:,3)-s(:,1).*s(:,2).*c(:,3), ...
    c(:,1).*s(:,2).*c(:,3)+s(:,1).*c(:,2).*s(:,3), ...
    s(:,1).*c(:,2).*c(:,3)-c(:,1).*s(:,2).*s(:,3)];


end


function eul = q2e( q )
%q2e Convert quaternion to Euler angles

% Normalize the quaternions
q = robotics.internal.normalizeRows(q);

qw = q(:,1);
qx = q(:,2);
qy = q(:,3);
qz = q(:,4);

% Pre-allocate output
eul = zeros(size(q,1), 3, 'like', q);

% The parsed sequence will be in all upper-case letters and validated
aSinInput = -2*(qx.*qz-qw.*qy);
aSinInput(aSinInput > 1) = 1;

eul = [ atan2( 2*(qx.*qy+qw.*qz), qw.^2 + qx.^2 - qy.^2 - qz.^2 ), ...
    asin( aSinInput ), ...
    atan2( 2*(qy.*qz+qw.*qx), qw.^2 - qx.^2 - qy.^2 + qz.^2 )];


% Check for complex numbers
if ~isreal(eul)
    eul = real(eul);
end
end
