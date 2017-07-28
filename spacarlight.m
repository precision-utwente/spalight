function [results] = spacarlight(varargin)
% SPACARLIGHT(nodes, elements, nprops, eprops, rls, opt)
% runs spacar simulation with reduced set of input arguments. See
% www.spacar.nl for more information.
%
% Created by: M. Nijenhuis and M. Naves
% Contact: m.naves@utwente.nl
%
% CONTRIBUTIONS
% J.P. Meijaard (complm)
% S.E. Boer (calc_stiffness, calc_inertia, calcTorsStiff and Spavisual functions)
% D.H. Wiersma (CWvalues)
%
% LIMITATIONS (note that the full Spacar version does allow these things)
% - Boundary conditions: the orientation of a node can either be fixed or
% free. It is not possible to create a pinned boundary condition about a certain axis.

% On each node only a single rotational input can be prescribed
% Output rotations are provided in quaternions and euler rotations in order z, y, x

% For certain desired simulations, the current feature set of 
% spacarlight() is too limited. In that case, the full version of Spacar 
% should be used. It offers *many* more features.

% NOTES
% Constrained warping is included by means of an effective torsional
% stiffness increase.

%% WARNINGS
warning off backtrace

%% Do not allow running this function as a script
ensure(size(dbstack,1)>1,'Call spacarlight() from a script instead.')

%% CHECK FOR INCOMPLETE INPUT
switch nargin
    case 0
        err('No input was provided.');
    case 1
        warning('Incomplete input; no simulation is run.');
        % validate only nodes
        [nodes] = validateInput(varargin{:});
        showInput(nodes);
        results = [];
        return
    case 2
        warning('Incomplete input; no simulation is run.');
        % validate only nodes and elements
        [nodes,elements] = validateInput(varargin{:});
        showInput(nodes,elements);
        results = [];
        return
    case 3
        warning('Incomplete input; no simulation is run.');
        % validate only nodes, elements and nprops
        [nodes,elements,nprops] = validateInput(varargin{:});
        showInput(nodes,elements,nprops);
        results = [];
        return
    case 4
        [nodes,elements,nprops,eprops] = validateInput(varargin{:});
        % attempt simulation
    case 5
        [nodes,elements,nprops,eprops,rls] = validateInput(varargin{:});
        % attempt simulation
    case 6
        [nodes,elements,nprops,eprops,rls,opt] = validateInput(varargin{:});
        if isfield(opt,'showinputonly') && opt.showinputonly == true
           showInput(nodes,elements,nprops,eprops);
           results = [];
           return
        end
        % attempt simulation
    otherwise
        err('Expecting a maximum of 6 input arguments.');
end

%% INITIALIZE VARIABLES, SET DEFAULTS (DO NOT SET DEFAULTS IN VALIDATEINPUT())
results     = [];
id_inputx   = false;                %identifier to check for prescribed input displacements/rotations
id_inputf   = false;                %identifier to check for external load
x_count     = size(nodes,1)*2+1;    %counter for node numbering
e_count     = 1;                    %counter for element numbering
X_list      = [];                   %list with node numbers
E_list      = [];                   %list with element numbers

%% CHECK EXISTENCE OF REQUIRED FUNCTIONS
if exist('spavisual','file')   ~=2;   err('spavisual() is not in your path.');                                         end
if exist('stressbeam','file')  ~=2;   err('stressbeam() is not in your path (typically part of spavisual package).');  end
if exist('spacar','file')      ~=3;   err('spacar() is not in your path.');                                            end
     
%% START CREATING DATFILE
fileID = fopen([opt.filename '.dat'],'w');

%% USERDEFINED NODES
fprintf(fileID,'#NODES\n');
%print all nodes provides by user
for i=1:size(nodes,1)
    fprintf(fileID,'X       %3u  %6f  %6f  %6f      #node %u\n',(i-1)*2+1,nodes(i,1),nodes(i,2),nodes(i,3),i);
end

%% ELEMENTS
fprintf(fileID,'\n\n#ELEMENTS\n');

for i=1:size(elements,1)
    fprintf(fileID,'#ELEMENT %u\n',i);
    
    % get element property set corresponding to element i
    for j=1:size(eprops,2)
        if sum(eprops(j).elems==i)>0 %element i is in property set j
            
            %element information
            N           = eprops(j).nbeams;    %number of beams per userdefined element
            Orien       = eprops(j).orien;      %orientation local y-vector
            Flex        = eprops(j).flex;       %flexibility of this element
            N_p         = elements(i,1);            %p-node nodenumber
            N_q         = elements(i,2);            %q-node nodenumber
            X_list(i,1) = N_p;                      %store p-node in X_list
            
            if N>1 %if more then 1 beam
                
                X_p = nodes(N_p,1:3);   %Location p-node
                X_q = nodes(N_q,1:3);   %Location q-node
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
            if ~exist('rls','var') || isempty(rls)            %if no rlse, add all flexible deformation modes as dyne
                fprintf(fileID,'dyne    %3u',e_count);
                for m=1:length(Flex)    %loop over all flexible deformation modes
                    fprintf(fileID,'  %3u',Flex(m));
                end
                fprintf(fileID,'\n');
                
            else%if some rls are specified
                %compensate size of rls if size is smaller then element list
                if i>size(rls,2)
                    rls(i).def = [];
                end
                
                % add dyne
                if ~isempty(Flex)                           %if some flexibility is specified
                    dyn_added = false;                      %reset identifier to check if string 'dyne' is added
                    for m=1:length(Flex)                    %loop over all flexible deformation modes
                        if ~(sum(rls(i).def==Flex(m))>0)   %if flexible deformation mode is not a rlse, it is dyne
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
                for m=1:length(rls(i).def)             %loop over all released deformation modes
                    if ~rlse_added                      %only add string 'rlse' if it is not yet added
                        fprintf(fileID,'rlse    %3u',e_count);
                        rlse_added = true;
                    end
                    fprintf(fileID,'  %3u',rls(i).def(m));
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
for i=1:size(nprops,2)
    %fixes
    if(isfield(nprops(i),'fix') && ~isempty(nprops(i).fix));    fprintf(fileID,'FIX      %u  \n',(i-1)*2+1);
        fprintf(fileID,'FIX      %u  \n',(i-1)*2+2);   end
    if(isfield(nprops(i),'fix_pos') && ~isempty(nprops(i).fix_pos));    fprintf(fileID,'FIX      %u  \n',(i-1)*2+1);   end
    if(isfield(nprops(i),'fix_x') && ~isempty(nprops(i).fix_x));        fprintf(fileID,'FIX      %u 1 \n',(i-1)*2+1);  end
    if(isfield(nprops(i),'fix_y') && ~isempty(nprops(i).fix_y));        fprintf(fileID,'FIX      %u 2 \n',(i-1)*2+1);  end
    if(isfield(nprops(i),'fix_z') && ~isempty(nprops(i).fix_z));        fprintf(fileID,'FIX      %u 3 \n',(i-1)*2+1);  end
    if(isfield(nprops(i),'fix_orien') && ~isempty(nprops(i).fix_orien));    fprintf(fileID,'FIX      %u  \n',(i-1)*2+2);   end
    
    %input displacements
    if((isfield(nprops(i),'displ_x') && ~isempty(nprops(i).displ_x)) ||...
            (isfield(nprops(i),'displ_initial_x') && ~isempty(nprops(i).displ_initial_x)));   fprintf(fileID,'INPUTX      %3u     1\n',(i-1)*2+1);id_inputx = true;    end
    if((isfield(nprops(i),'displ_y') && ~isempty(nprops(i).displ_y)) ||...
            (isfield(nprops(i),'displ_initial_y') && ~isempty(nprops(i).displ_initial_y)));   fprintf(fileID,'INPUTX      %3u     2\n',(i-1)*2+1);id_inputx = true;    end
    if((isfield(nprops(i),'displ_z') && ~isempty(nprops(i).displ_z)) ||...
            (isfield(nprops(i),'displ_initial_z') && ~isempty(nprops(i).displ_initial_z)));   fprintf(fileID,'INPUTX      %3u     3\n',(i-1)*2+1);id_inputx = true;    end
    
    %input rotations
    id_inputr = 0; %identifier to count the number of input rotations
    if((isfield(nprops(i),'rot_x') && ~isempty(nprops(i).rot_x)) ||...
            (isfield(nprops(i),'rot_initial_x') && ~isempty(nprops(i).rot_initial_x))); fprintf(fileID,'INPUTX      %3u     4\n',(i-1)*2+2);id_inputx = true; id_inputr=id_inputr+1; end
    if((isfield(nprops(i),'rot_y') && ~isempty(nprops(i).rot_y)) ||...
            (isfield(nprops(i),'rot_initial_y') && ~isempty(nprops(i).rot_initial_y))); fprintf(fileID,'INPUTX      %3u     3\n',(i-1)*2+2);id_inputx = true; id_inputr=id_inputr+1; end
    if((isfield(nprops(i),'rot_z') && ~isempty(nprops(i).rot_z)) ||...
            (isfield(nprops(i),'rot_initial_z') && ~isempty(nprops(i).rot_initial_z))); fprintf(fileID,'INPUTX      %3u     2\n',(i-1)*2+2);id_inputx = true; id_inputr=id_inputr+1; end
    if id_inputr>1 %if multiple rotations are prescribed, problems can arise with quaternion<->euler conversion
        err('Multiple rotational inputs defined for node %u. Only a single input rotation can be added to a node.',i)
    end
end
fprintf(fileID,'\n\nEND\nHALT\n\n');


%% STIFFNESS/INERTIA PROPS
for i=1:size(eprops,2) %loop over each element property set
    
    %only write stiffness and mass values when deformations are flexible
    if (isfield(eprops(i),'flex') && ~isempty(eprops(i).flex))
        stiffness = calc_stiffness(eprops(i)); %calculate stiffness values
        inertia = calc_inertia(eprops(i));     %calculate mass properties
        for j=1:length(eprops(i).elems)                        %loop over all elemenents in element property set
            L   = norm(nodes(elements(eprops(i).elems(j),2),:)...
                - nodes(elements(eprops(i).elems(j),1),:));    %calculate flexure length for constrained warping values
            cw  = CWvalues(L,eprops(i));                        %calculate constrained warping values
            for k=1:size(E_list,2)                                  %loop over all beams in the element
                El = E_list(eprops(i).elems(j),k);
                if El>0
                    fprintf(fileID,'ESTIFF      %u %f %f %f %f %f %f \n',El,stiffness(1),cw*stiffness(2),stiffness(3),stiffness(4),stiffness(5),stiffness(6));
                end
            end
            for k=1:size(E_list,2) %write mass/inertia values
                El = E_list(eprops(i).elems(j),k); %loop over all beams in element set
                if El>0
                    fprintf(fileID,'EM      %u %f %f %f %f %f \n',El,inertia(1),inertia(2),inertia(3),inertia(4),inertia(5));
                end
            end
        end
    end
   
end


%% FORCES/MOMENTS/NODAL MASSES
fprintf(fileID,'\n\n#FROCES/MOMENTS\n');
id_ini  = false; %check for initial loading or displacement
id_add  = false; %check for aditional loading or displacement

for i=1:size(nprops,2) %loop over all user defined nodes
    %forces
    if(isfield(nprops(i),'force') && ~isempty(nprops(i).force));                    fprintf(fileID,'DELXF   %3u %6f %6f %6f  \n',(i-1)*2+1,nprops(i).force(1),nprops(i).force(2),nprops(i).force(3));                           id_add = true;  id_inputf=true; end
    if(isfield(nprops(i),'force_initial') && ~isempty(nprops(i).force_initial));    fprintf(fileID,'XF      %3u %6f %6f %6f  \n',(i-1)*2+1,nprops(i).force_initial(1),nprops(i).force_initial(2),nprops(i).force_initial(3));   id_ini = true;  id_inputf=true; end
    
    %moments
    if(isfield(nprops(i),'moment') && ~isempty(nprops(i).moment)) %#ok<*ALIGN>
      moments = nprops(i).moment;
                                                                                    fprintf(fileID,'DELXF   %3u %6f %6f %6f %6f  \n',(i-1)*2+2,0,2*moments(1),2*moments(2),2*moments(3));                                       id_add = true;  id_inputf=true; end
    if(isfield(nprops(i),'moment_initial') && ~isempty(nprops(i).moment_initial))
      moments_i = nprops(i).moment_initial;
                                                                                    fprintf(fileID,'XF      %3u %6f %6f %6f %6f  \n',(i-1)*2+2,0,2*moments_i(1),2*moments_i(2),2*moments_i(3));                                 id_ini = true;  id_inputf=true; end
    
    %displacements
    if(isfield(nprops(i),'displ_x') && ~isempty(nprops(i).displ_x));                  fprintf(fileID,'DELINPX  %3u  1  %6f  \n',(i-1)*2+1,nprops(i).displ_x(1));                       id_add = true; end
    if(isfield(nprops(i),'displ_y') && ~isempty(nprops(i).displ_y));                  fprintf(fileID,'DELINPX  %3u  2  %6f  \n',(i-1)*2+1,nprops(i).displ_y(1));                       id_add = true; end
    if(isfield(nprops(i),'displ_z') && ~isempty(nprops(i).displ_z));                  fprintf(fileID,'DELINPX  %3u  3  %6f  \n',(i-1)*2+1,nprops(i).displ_z(1));                       id_add = true; end
    if(isfield(nprops(i),'displ_initial_x') && ~isempty(nprops(i).displ_initial_x));  fprintf(fileID,'INPUTX   %3u  1  %6f  \n',(i-1)*2+1,nodes(i,1) + nprops(i).displ_initial_x(1));  id_ini = true; end
    if(isfield(nprops(i),'displ_initial_y') && ~isempty(nprops(i).displ_initial_y));  fprintf(fileID,'INPUTX   %3u  2  %6f  \n',(i-1)*2+1,nodes(i,2) + nprops(i).displ_initial_y(1));  id_ini = true; end
    if(isfield(nprops(i),'displ_initial_z') && ~isempty(nprops(i).displ_initial_z));  fprintf(fileID,'INPUTX   %3u  3  %6f  \n',(i-1)*2+1,nodes(i,3) +nprops(i).displ_initial_z(1));   id_ini = true; end
    
    %rotations
    if(isfield(nprops(i),'rot_x') && ~isempty(nprops(i).rot_x));                rot = eul2quat([nprops(i).rot_x(1) 0 0]);
        fprintf(fileID,'DELINPX     %3u     4   %6f  \n',(i-1)*2+2,rot(4)); id_add = true; end
    if(isfield(nprops(i),'rot_y') && ~isempty(nprops(i).rot_y));                rot = eul2quat([0 nprops(i).rot_y(1) 0]);
        fprintf(fileID,'DELINPX     %3u     3   %6f  \n',(i-1)*2+2,rot(3)); id_add = true; end
    if(isfield(nprops(i),'rot_z') && ~isempty(nprops(i).rot_z));                rot = eul2quat([0 0 nprops(i).rot_z(1)]);
        fprintf(fileID,'DELINPX     %3u     2   %6f  \n',(i-1)*2+2,rot(2)); id_add = true; end
    if(isfield(nprops(i),'rot_initial_x') && ~isempty(nprops(i).rot_initial_x));rot = eul2quat([nprops(i).rot_initial_x(1) 0 0]);
        fprintf(fileID,'INPUTX     %3u     4   %6f  \n',(i-1)*2+2,rot(4));  id_ini = true; end
    if(isfield(nprops(i),'rot_initial_y') && ~isempty(nprops(i).rot_initial_y));rot = eul2quat([0 nprops(i).rot_initial_y(1) 0]);
        fprintf(fileID,'INPUTX     %3u     3   %6f  \n',(i-1)*2+2,rot(3));  id_ini = true; end
    if(isfield(nprops(i),'rot_initial_z') && ~isempty(nprops(i).rot_initial_z));rot = eul2quat([0 0 nprops(i).rot_initial_z(1)]);
        fprintf(fileID,'INPUTX     %3u     2   %6f  \n',(i-1)*2+2,rot(2));  id_ini = true; end
    
    %nodal masses/inertia
    if(isfield(nprops(i),'mass') && ~isempty(nprops(i).mass));                      fprintf(fileID,'XM      %3u %6f  \n',(i-1)*2+1,nprops(i).mass); end
    if(isfield(nprops(i),'mominertia') && ~isempty(nprops(i).mominertia));                fprintf(fileID,'XM      %3u %6f %6f %6f %6f %6f %6f\n',(i-1)*2+2,nprops(i).mominertia(1),nprops(i).mominertia(2),nprops(i).mominertia(3),...
            nprops(i).mominertia(4),nprops(i).mominertia(5),nprops(i).mominertia(6)); end
end


%% ADITIONAL OPTIONS
%GRAVITY
if(fieldexist('opt','gravity') && ~isempty(opt.gravity)); fprintf(fileID,'\nGRAVITY  %6f %6f %6f',opt.gravity(1),opt.gravity(2),opt.gravity(3)); end

%ITERSTEP SETTINGS
if      (id_ini && id_add);      fprintf(fileID,'\nITERSTEP 10 10 0.0000005 1 3 10');    %if initial and aditional loading/displacement
elseif  (id_ini && ~id_add);     fprintf(fileID,'\nITERSTEP 10 1 0.0000005  1 1 10');    %if initial loading/displacement
elseif  (~id_ini && id_add);     fprintf(fileID,'\nITERSTEP 10 10 0.0000005 1 3 0');     %if initial loading/displacement
else                             fprintf(fileID,'\nITERSTEP 10 1  0.0000005 1 1 0');  end%#ok<SEPEX> %no loading/displacement

% %TRANSFER FUNCTION INPUT/OUTPUT
% if ((isfield(opt,'transfer_in') && ~isempty(opt.transfer_in)) ||  (isfield(opt,'transfer_out') && ~isempty(opt.transferout)))
%     if id_inputx
%         disp('Warning: input displacement is prediscribed, possibly affecting input/ouput transfer function.')
%     end
%     
%     fprintf(fileID,'\n\nEND\nHALT\n\n');
%     for i=1:size(opt.transfer_in,2) %add inputs
%         switch opt.transfer_in(i).type
%             case 'force_x';     fprintf(fileID,'\nINPUTF %2u %3u 1',i,(opt.transfer_in(i).node-1)*2+1);
%             case 'force_y';     fprintf(fileID,'\nINPUTF %2u %3u 2',i,(opt.transfer_in(i).node-1)*2+1);
%             case 'force_z';     fprintf(fileID,'\nINPUTF %2u %3u 3',i,(opt.transfer_in(i).node-1)*2+1);
%                 %TO BE DONE
%                 %case 'moment_x';    fprintf(fileID,'\nINPUTF %2u %3u 4',i,(opt.transfer_in(i).node-1)*2+2);
%                 %case 'moment_y';    fprintf(fileID,'\nINPUTF %2u %3u 3',i,(opt.transfer_in(i).node-1)*2+2);
%                 %case 'moment_z';    fprintf(fileID,'\nINPUTF %2u %3u 2',i,(opt.transfer_in(i).node-1)*2+2);
%                 
%             case 'displ_x';      fprintf(fileID,'\nINX %2u %3u 1',i,(opt.transfer_in(i).node-1)*2+1);
%             case 'displ_y';      fprintf(fileID,'\nINX %2u %3u 2',i,(opt.transfer_in(i).node-1)*2+1);
%             case 'displ_z';      fprintf(fileID,'\nINX %2u %3u 3',i,(opt.transfer_in(i).node-1)*2+1);
%                 %case 'rot_x';       fprintf(fileID,'\nINX %2u %3u 4',i,(opt.transfer_in(i).node-1)*2+2);
%                 %case 'rot_y';       fprintf(fileID,'\nINX %2u %3u 3',i,(opt.transfer_in(i).node-1)*2+2);
%                 %case 'rot_z';       fprintf(fileID,'\nINX %2u %3u 2',i,(opt.transfer_in(i).node-1)*2+2);
%         end
%     end
%     for i=1:size(opt.transfer_out,2) %add outputs
%         switch opt.transfer_out(i).type
%             case 'force_x';     fprintf(fileID,'\nOUTF %2u %3u 1',i,(opt.transfer_out(i).node-1)*2+1);
%             case 'force_y';     fprintf(fileID,'\nOUTF %2u %3u 2',i,(opt.transfer_out(i).node-1)*2+1);
%             case 'force_z';     fprintf(fileID,'\nOUTF %2u %3u 3',i,(opt.transfer_out(i).node-1)*2+1);
%                 %TO BE DONE
%                 %case 'moment_x';    fprintf(fileID,'\nOUTF %2u %3u 4',i,(opt.transfer_out(i).node-1)*2+2);
%                 %case 'moment_y';    fprintf(fileID,'\nOUTF %2u %3u 3',i,(opt.transfer_out(i).node-1)*2+2);
%                 %case 'moment_z';    fprintf(fileID,'\nOUTF %2u %3u 2',i,(opt.transfer_out(i).node-1)*2+2);
%                 
%             case 'displ_x';      fprintf(fileID,'\nOUTX %2u %3u 1',i,(opt.transfer_out(i).node-1)*2+1);
%             case 'displ_y';      fprintf(fileID,'\nOUTX %2u %3u 2',i,(opt.transfer_out(i).node-1)*2+1);
%             case 'displ_z';      fprintf(fileID,'\nOUTX %2u %3u 3',i,(opt.transfer_out(i).node-1)*2+1);
%                 %case 'rot_x';       fprintf(fileID,'\nOUTX %2u %3u 4',i,(opt.transfer_out(i).node-1)*2+2);
%                 %case 'rot_y';       fprintf(fileID,'\nOUTX %2u %3u 3',i,(opt.transfer_out(i).node-1)*2+2);
%                 %case 'rot_z';       fprintf(fileID,'\nOUTX %2u %3u 2',i,(opt.transfer_out(i).node-1)*2+2);
%         end
%     end
% end

fprintf(fileID,'\n\nEND\nEND\n\n');


%% VISUALIZATION
fprintf(fileID,'\n\nVISUALIZATION');
for i=1:size(eprops,2) %loop over all element property sets
    
    %CROSSECTIONAL PROPERTIES
    fprintf(fileID,'\n\nBEAMPROPS ');
    for j=1:length(eprops(i).elems) %loop over all elemenents in element set i
        for k=1:size(E_list,2)
            El = E_list(eprops(i).elems(j),k);
            if El>0
                fprintf(fileID,' %u ',El);
            end
        end
    end
    
    if isempty(eprops(i).type)
        fprintf(fileID,'\nCROSSTYPE  RECT');
        fprintf(fileID,'\nCROSSDIM 0.05 0.05');
    else
        switch eprops(i).type
            case {'leafspring','rigid'} %if leafspring or rigid, do rect crossection
                fprintf(fileID,'\nCROSSTYPE  RECT');
                fprintf(fileID,'\nCROSSDIM  %f  %f',eprops(i).dim(1),eprops(i).dim(2));
            case 'wire'                 %if wire, do circular crossection
                fprintf(fileID,'\nCROSSTYPE  CIRC');
                fprintf(fileID,'\nCROSSDIM  %f ',eprops(i).dim(1));
        end
    end
    
    %COLOR
    if (isfield(eprops(i),'color') && ~isempty(eprops(i).color))
        fprintf(fileID,'\nGRAPHICS ');
        for j=1:length(eprops(i).elems)
            for k=1:size(E_list,2)
                El = E_list(eprops(i).elems(j),k);
                if El>0
                    fprintf(fileID,' %u ',El);
                end
            end
        end
        fprintf(fileID,'\nFACECOLOR  %f %f %f ',eprops(i).color(1),eprops(i).color(2),eprops(i).color(3));
    end
    
    %VISIBILITY
    if (isfield(eprops(i),'hide') && ~isempty(eprops(i).hide) && eprops(i).hide==1)
        fprintf(fileID,'\nDONOTDRAW ');
        for j=1:length(eprops(i).elems)
            for k=1:size(E_list,2)
                El = E_list(eprops(i).elems(j),k);
                if El>0
                    fprintf(fileID,' %u ',El);
                end
            end
        end
    end
end
fclose(fileID); %datfile finished!


%% SIMULATE CONSTRAINTS
if ~(fieldexist('opt','silent') && opt.silent==1)
    try 
        warning('off','all')
        spacar(0,opt.filename)
        warning('on','all')
    catch
        warning('on','all') %needed here, since a spacar error in the try block would leave warnings off
        err('Connectivity incorrect. Check element lengths, eprops.orien, etc.');
    end
    
    %CHECK CONSTRAINTS
    sbd     = [opt.filename '.sbd'];
    nep     = getfrsbf(sbd,'nep');
    nxp     = getfrsbf(sbd,'nxp');
    nddof   = getfrsbf(sbd,'nddof');
    le      = getfrsbf(sbd,'le');
    BigD    = getfrsbf(sbd,'bigd',1);
    Dcc     = BigD( 1:(nep(1)+nep(3)+nep(4)) , nxp(1)+(1:nxp(2)) );
    [ U, s, V ] = svd(Dcc);
    s       = diag(s);
    
    if isempty(s) %%% TO BE DONE: s can be empty. What does this mean?
        return
    end
    
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
        warning('System is underconstrained. Check element connectivity, boundary conditions and releases.')
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

        %Reshape rls suggestions according to user defined elements
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

        results.overconstraints = [OC_el OC_defs];
        warning('System is overconstrained; releases are required in order to run static simulation.\nA suggestion for possible releases is given in results.overconstraints in the workspace and the table below.\n')
        fprintf('Number of overconstraints: %u\n\n',nover);
        disp(table(OC_el,sum((OC_defs==1),2),sum((OC_defs==2),2),sum((OC_defs==3),2),sum((OC_defs==4),2),sum((OC_defs==5),2),sum((OC_defs==6),2),...
            'VariableNames',{'Element' 'def_1' 'def_2 ' 'def_3' 'def_4' 'def_5' 'def_6'}));
        warning('on','all')
        return
    end
    
    if nddof == 0
        warning('The system has no degrees of freedom (so no Spacar simulation will be performed). Check eprops.flex and rls.')
        return;
    end
        
end

%% SIMULATE STATICS
try
    warning('off','all')
    spacar(-10,opt.filename)
    warning('on','all')
    if ~(fieldexist('opt','silent') && opt.silent==1)
        spavisual(opt.filename)
    end
    disp('Spacar simulation succeeded.')
catch
    warning('on','all') %needed here, since a spacar error in the try block would leave warnings off
    warning('Spacar simulation failed. Possibly failed to converge to solution. Check magnitude of input displacements, loads and other input data.')
end
try
    %get results
    results = calc_Results(opt.filename, E_list, id_inputf, id_inputx, nodes, elements, nprops, eprops, rls, opt);
catch
    err('A problem occurred processing simulation results.')
end

%% WARNINGS
warning backtrace on

%% END OF SPACAR_LIGHT
end

function varargout = validateInput(varargin)

    %TO BE DONE
    %Check input voor input momenten
    switch nargin
        case 1
            nodes = varargin{1};
        case 2
            nodes = varargin{1};
            elements = varargin{2};
        case 3
            nodes = varargin{1};
            elements = varargin{2};
            nprops = varargin{3};
        case 4
            nodes = varargin{1};
            elements = varargin{2};
            nprops = varargin{3};
            eprops = varargin{4};
        case 5
            nodes = varargin{1};
            elements = varargin{2};
            nprops = varargin{3};
            eprops = varargin{4};
            rls = varargin{5};
        case 6
            nodes = varargin{1};
            elements = varargin{2};
            nprops = varargin{3};
            eprops = varargin{4};
            rls = varargin{5};
            opt = varargin{6};
    end
    
    %DO NOT PERFORM CHECKS IN SILENT MODE
    if ~(fieldexist('opt','silent') && opt.silent==1)

        %CHECK NODES INPUT VARIABLE
        validateattributes(nodes,   {'double'},{'ncols',3,'ndims',2},'','nodes')

        %CHECK ELEMENTS INPUT VARIABLE
        if exist('elements','var')
            nno = size(nodes,1);
            validateattributes(elements,{'double'},{'ncols',2,'ndims',2},'','elements')
            ensure(all(elements(:) == floor(elements(:))),'Entries in elements seem to be non-integers.')
            
            ensure(all(elements(:)>0),'Element seems connected to node number <=0.')
            ensure(max(elements(:))<=nno,'Element seems connected to node that does not exist.')
            if max(elements(:))<nno; warning('Node seems not attached to element.'); end
            ensure(~any(abs(elements(:,1)-elements(:,2))==0),('Both sides of element seem connected to the same node.'))
            
            %check if unique pairs node numbers (independent of p/q order)
            ensure(size(unique(sort(elements,2),'rows'),1)==size(elements,1),'Multiple elements seem connected between the same node pair.')
            
            ensure(all(sqrt(sum((nodes(elements(:,1),:) - nodes(elements(:,2),:)).^2,2))>1e-5),'Element length seems smaller than 0.00001.')
            
        end

        %CHECK NPROPS INPUT VARIABLE
        if exist('nprops','var')
            allowed_nprops = {'fix','fix_x','fix_y','fix_z','fix_pos','fix_orien','displ_x','displ_y','displ_z','rot_x','rot_y','rot_z','force','moment','mass','mominertia','force_initial','moment_initial', ...
                'displ_initial_x','displ_initial_y','displ_initial_z','rot_initial_rx','rot_initial_ry','rot_initial_rz'};
            supplied_nprops = fieldnames(nprops);
            unknown_nprops_i = ~ismember(supplied_nprops,allowed_nprops);
            if any(unknown_nprops_i)
                err('Unknown nprops field %s',supplied_nprops{unknown_nprops_i});
            end
            
            validateattributes(nprops,{'struct'},{'nonempty'},'','nprops')
            count_bcs = 0; %counter for total number of constraints (there should be at least 6)
            Node_fields = fieldnames(nprops);
            for i=1:size(nprops,2)
                for j=1:length(Node_fields)
                    switch Node_fields{j}
                        case 'fix'
                            if ~isempty(nprops(i).(Node_fields{j}));    validateattributes(nprops(i).(Node_fields{j}),{'logical'},{'scalar'},'',            sprintf('fix property in nprops(%u)',i));       end
                        case {'fix_pos','fix_orien'}
                            if ~isempty(nprops(i).(Node_fields{j}));    validateattributes(nprops(i).(Node_fields{j}),{'logical'},{'scalar'},'',            sprintf('fix_pos/orien property in nprops(%u)',i)); end
                        case {'fix_x','fix_y','fix_z'}
                            if ~isempty(nprops(i).(Node_fields{j}));    validateattributes(nprops(i).(Node_fields{j}),{'logical'},{'scalar'},'',            sprintf('fix_x/y/z property in nprops(%u)',i)); end
                        case {'force','force_initial'}
                            if ~isempty(nprops(i).(Node_fields{j}));     validateattributes(nprops(i).(Node_fields{j}),{'double'},{'vector','numel',3},'',   sprintf('force property in nprops(%u)',i));     end
                        case {'moment','moment_initial'}
                            if ~isempty(nprops(i).(Node_fields{j}));     validateattributes(nprops(i).(Node_fields{j}),{'double'},{'vector','numel',3},'',   sprintf('moment property in nprops(%u)',i));     end    
                        case {'displ_x','displ_y','displ_z','displ_initial_x','displ_initial_y','displ_initial_z'}
                            if ~isempty(nprops(i).(Node_fields{j}));     validateattributes(nprops(i).(Node_fields{j}),{'double'},{'scalar'},'',             sprintf('displ property in nprops(%u)',i));      end
                        case {'rot_x','rot_y','rot_z','rot_initial_x','rot_initial_y','rot_initial_z'}
                            if ~isempty(nprops(i).(Node_fields{j}));     validateattributes(nprops(i).(Node_fields{j}),{'double'},{'scalar'},'',             sprintf('rot property in nprops(%u)',i));      end
                        case 'mass'
                            if ~isempty(nprops(i).(Node_fields{j}));     validateattributes(nprops(i).(Node_fields{j}),{'double'},{'scalar'},'',             sprintf('mass property in nprops(%u)',i));      end
                        case 'mominertia'
                            if ~isempty(nprops(i).(Node_fields{j}));     validateattributes(nprops(i).(Node_fields{j}),{'double'},{'vector','numel',6},'',   sprintf('mominertia property in nprops(%u)',i));   end
                    end
                end
                
                %check for total number of constraints supplied (there should be at least 6)
                    %checks for fixes
                    if(isfield(nprops(i),'fix') && ~isempty(nprops(i).fix) && nprops(i).fix == true);                   count_bcs = count_bcs + 6;   end
                    if(isfield(nprops(i),'fix_pos') && ~isempty(nprops(i).fix_pos) && nprops(i).fix_pos == true);       count_bcs = count_bcs + 3;   end
                    if(isfield(nprops(i),'fix_orien') && ~isempty(nprops(i).fix_orien) && nprops(i).fix_orien == true); count_bcs = count_bcs + 3;   end
                    if(isfield(nprops(i),'fix_x') && ~isempty(nprops(i).fix_x) && nprops(i).fix_x == true);             count_bcs = count_bcs + 1;   end
                    if(isfield(nprops(i),'fix_y') && ~isempty(nprops(i).fix_y) && nprops(i).fix_y == true);             count_bcs = count_bcs + 1;   end
                    if(isfield(nprops(i),'fix_z') && ~isempty(nprops(i).fix_z) && nprops(i).fix_z == true);             count_bcs = count_bcs + 1;   end

                    %checks for input displacements
                    if((isfield(nprops(i),'displ_x') && ~isempty(nprops(i).displ_x)) ||...
                            (isfield(nprops(i),'displ_initial_x') && ~isempty(nprops(i).displ_initial_x)));     count_bcs = count_bcs + 1;   end 
                    if((isfield(nprops(i),'displ_y') && ~isempty(nprops(i).displ_y)) ||...
                            (isfield(nprops(i),'displ_initial_y') && ~isempty(nprops(i).displ_initial_y)));     count_bcs = count_bcs + 1;   end
                    if((isfield(nprops(i),'displ_z') && ~isempty(nprops(i).displ_z)) ||...
                            (isfield(nprops(i),'displ_initial_z') && ~isempty(nprops(i).displ_initial_z)));     count_bcs = count_bcs + 1;   end

                    %checks for input rotations
                    if((isfield(nprops(i),'rot_x') && ~isempty(nprops(i).rot_x)) ||...
                            (isfield(nprops(i),'rot_initial_x') && ~isempty(nprops(i).rot_initial_x)));         count_bcs = count_bcs + 1;   end
                    if((isfield(nprops(i),'rot_y') && ~isempty(nprops(i).rot_y)) ||...
                            (isfield(nprops(i),'rot_initial_y') && ~isempty(nprops(i).rot_initial_y)));         count_bcs = count_bcs + 1;   end
                    if((isfield(nprops(i),'rot_z') && ~isempty(nprops(i).rot_z)) ||...
                            (isfield(nprops(i),'rot_initial_z') && ~isempty(nprops(i).rot_initial_z)));         count_bcs = count_bcs + 1;   end
                    
                %no combination of fix, force, displ on one node in the same direction (translational along x)
                ensure(sum([ ...
                (isfield(nprops(i),'fix_x') && ~isempty(nprops(i).fix_x) && nprops(i).fix_x == true) ...
                ((isfield(nprops(i),'force') && ~isempty(nprops(i).force) && nprops(i).force(1)~=0) || (isfield(nprops(i),'force_initial') && ~isempty(nprops(i).force_initial) && nprops(i).force_initial(1)~=0 )) ...
                ((isfield(nprops(i),'displ_x') && ~isempty(nprops(i).displ_x)) || (isfield(nprops(i),'displ_initial_x') && ~isempty(nprops(i).displ_initial_x))) ...
                    ])<=1,'There is a combination of fix_x, force(1) and displ_(initial_)x on node %i.',i);
                
                %no combination of fix, force, displ on one node in the same direction (translational along y)
                ensure(sum([ ...
                (isfield(nprops(i),'fix_y') && ~isempty(nprops(i).fix_y) && nprops(i).fix_y == true) ...
                ((isfield(nprops(i),'force') && ~isempty(nprops(i).force) && nprops(i).force(2)~=0) || (isfield(nprops(i),'force_initial') && ~isempty(nprops(i).force_initial) && nprops(i).force_initial(2)~=0 )) ...
                ((isfield(nprops(i),'displ_y') && ~isempty(nprops(i).displ_y)) || (isfield(nprops(i),'displ_initial_y') && ~isempty(nprops(i).displ_initial_y))) ...
                    ])<=1,'There is a combination of fix_y, force(2) and displ_(initial_)y on node %i.',i);
                
                %no combination of fix, force, displ on one node in the same direction (translational along z)
                ensure(sum([ ...
                (isfield(nprops(i),'fix_z') && ~isempty(nprops(i).fix_z) && nprops(i).fix_z == true) ...
                ((isfield(nprops(i),'force') && ~isempty(nprops(i).force) && nprops(i).force(3)~=0) || (isfield(nprops(i),'force_initial') && ~isempty(nprops(i).force_initial) && nprops(i).force_initial(3)~=0 )) ...
                ((isfield(nprops(i),'displ_z') && ~isempty(nprops(i).displ_z)) || (isfield(nprops(i),'displ_initial_z') && ~isempty(nprops(i).displ_initial_z))) ...
                    ])<=1,'There is a combination of fix_z, force(3) and displ_(initial_)z on node %i.',i);
                
                %no combination of fix_orien, moment, rot on one node
                ensure(sum([ ...
                (isfield(nprops(i),'fix_orien') && ~isempty(nprops(i).fix_orien) && nprops(i).fix_orien == true) ...
                ((isfield(nprops(i),'moment') && ~isempty(nprops(i).moment) && any(nprops(i).moment~=0)) || (isfield(nprops(i),'moment_initial') && ~isempty(nprops(i).moment_initial) && any(nprops(i).moment_initial~=0) )) ...
                ( (isfield(nprops(i),'rot_x') && ~isempty(nprops(i).rot_x)) || (isfield(nprops(i),'rot_y') && ~isempty(nprops(i).rot_y)) || (isfield(nprops(i),'rot_z') && ~isempty(nprops(i).rot_z)) || ...
                  (isfield(nprops(i),'rot_initial_x') && ~isempty(nprops(i).rot_initial_x)) || (isfield(nprops(i),'rot_initial_y') && ~isempty(nprops(i).rot_initial_y)) || (isfield(nprops(i),'rot_initial_z') && ~isempty(nprops(i).rot_initial_z)) ) ...
                    ])<=1,'There is a combination of fix_orien, moment and rot_x/y/z on node %i.',i);
                
                %no combination of fix (6 constraints) and (mass or mominertia)
                if (isfield(nprops(i),'fix') && ~isempty(nprops(i).fix) && nprops(i).fix == true) && ...
                        (   (isfield(nprops(i),'mass') && ~isempty(nprops(i).mass) && nprops(i).mass~=0) || ...
                            (isfield(nprops(i),'mominertia') && ~isempty(nprops(i).mominertia) && any(nprops(i).mominertia~=0)) ...
                        )
                        
                    warning('Inertia associated with fixed node %i.',i);
                end
                
                %no combination of fix_pos and mass
                if (isfield(nprops(i),'fix_pos') && ~isempty(nprops(i).fix_pos) && nprops(i).fix_pos == true) && ...
                   (isfield(nprops(i),'mass') && ~isempty(nprops(i).mass) && nprops(i).mass~=0)
                        
                    warning('Mass associated with position-fixed node %i.',i);
                end                
                
                %no combination of fix_orien and mominertia
                if (isfield(nprops(i),'fix_orien') && ~isempty(nprops(i).fix_orien) && nprops(i).fix_orien == true) && ...
                   (isfield(nprops(i),'mominertia') && ~isempty(nprops(i).mominertia) && any(nprops(i).mominertia~=0))
                        
                    warning('Moment of inertia associated with orientation-fixed node %i.',i);
                end   
            end
            ensure(count_bcs >= 6,'The nodes seem to have insufficient (%i<6) constraints (fix, displ, or rot).',count_bcs);
        end
        
        %CHECK EPROPS INPUT VARIABLE
        if exist('eprops','var')
            allowed_eprops = {'elems','emod','smod','dens','type','dim','orien','nbeams','flex','color','hide'};
            supplied_eprops = fieldnames(eprops);
            unknown_eprops_i = ~ismember(supplied_eprops,allowed_eprops);
            if any(unknown_eprops_i)
               err('Unknown eprops field %s.',supplied_eprops{unknown_eprops_i});
            end
            
            el_nr_doubles_check = []; %filling this with user defined element numbers (with .elems) to check for doubles and missing elements
            for i=1:size(eprops,2)

                %the only mandatory field is elems
                if ~(isfield(eprops(i),'elems') && ~isempty(eprops(i).elems))
                    warning('Property elems is not defined in eprops(%u); ignoring eprops(%i).',i,i);
                    eprops(i) = []; %note, this deletes the current property set from the list
                else
                    %%%%%%%%%%%%
                    %elems exists, so validate elems
                    validateattributes(eprops(i).elems,{'double'},{'vector'},'',sprintf('elems property in eprops(%u)',i));
                    validateattributes(eprops(i).elems,{'double'},{'positive'},'',sprintf('elems property in eprops(%u)',i));

                    %check if elems for set i are unique
                    if length(unique(eprops(i).elems)) < length(eprops(i).elems)
                        err('Property elems of eprops(%i) contains non-unique element numbers.',i)
                    end

                    %check if elems only contains defined elements
                    undef_el_index = ~ismember(eprops(i).elems,1:size(elements,1));
                    if any(undef_el_index)
                        err('eprops(%i).elems contains undefined element number(s).',i)
                    end

                    %check if elements in elems are not already defined previously
                    double_el_index = ismember(eprops(i).elems,el_nr_doubles_check);
                    if any(double_el_index)
                        err('eprops(%i).elems contains element(s) whose properties have already been defined.',i)
                    end
                    el_nr_doubles_check = [el_nr_doubles_check; eprops(i).elems(:)];
                    %end validation of elems
                    %%%%%%%%%%%
                    
                    %%%%%%%%%%%%%%%
                    %validate other properties that are specified
                    if (isfield(eprops(i),'emod') && ~isempty(eprops(i).emod));     validateattributes(eprops(i).emod,{'double'},{'scalar'},'',   sprintf('emod property in eprops(%u)',i)); end
                    if (isfield(eprops(i),'smod') && ~isempty(eprops(i).smod));     validateattributes(eprops(i).smod,{'double'},{'scalar'},'',   sprintf('smod property in eprops(%u)',i)); end
                    if (isfield(eprops(i),'dens') && ~isempty(eprops(i).dens));     validateattributes(eprops(i).dens,{'double'},{'scalar'},'',   sprintf('dens property in eprops(%u)',i)); end
                    if (isfield(eprops(i),'dim') && ~isempty(eprops(i).dim));       validateattributes(eprops(i).dim,{'double'},{'vector'},'',    sprintf('dim property in eprops(%u)',i));  end
                    if (isfield(eprops(i),'color') && ~isempty(eprops(i).color));   validateattributes(eprops(i).color,{'double'},{'vector','numel',3},'',sprintf('color property in eprops(%u)',i)); end
                    if (isfield(eprops(i),'hide') && ~isempty(eprops(i).hide));     validateattributes(eprops(i).hide,{'logical'},{'scalar'},'',sprintf('hide property in eprops(%u)',i)); end
                    
                    if (isfield(eprops(i),'type') && ~isempty(eprops(i).type))
                        validateattributes(eprops(i).type,{'char'},{'nonempty'},'',sprintf('type property in eprops(%u)',i));
                        if ~any(strcmp(eprops(i).type,{'wire','leafspring','rigid'})), err('Element type should be either leafspring, wire or rigid.'); end
                    end
                    
                    if (isfield(eprops(i),'nbeams') && ~isempty(eprops(i).nbeams))
                        validateattributes(eprops(i).nbeams,{'double'},{'scalar','>=',1},'',sprintf('nbeams property in eprops(%u)',i));
                    else
                        eprops(i).nbeams = 1; 
                    end

                    if (isfield(eprops(i),'orien') && ~isempty(eprops(i).orien))   
                        validateattributes(eprops(i).orien,{'double'},{'vector','numel',3},'',sprintf('orien property in eprops(%u)',i));
                        
                        for j=1:length(eprops(i).elems)
                            xp = nodes(elements(eprops(i).elems(j),1),:);
                            xq = nodes(elements(eprops(i).elems(j),2),:);
                            
                            %check if supplied local y vector works
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %this directly from SPACAR 2015 source:
                            ex = xq(:) - xp(:);
                            ex = ex/norm(ex);
                            ey_input = eprops(i).orien;
                            ey = ey_input(:)/norm(ey_input);
                            ex_proj = dot(ey,ex);
                            noemer = sqrt(1-ex_proj^2);
                            if noemer < 1e-5
                                err('Orien property of element %i does not work, because (almost) parallel to element axis.',eprops(i).elems(j))
                            end
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            %check if normal to local x axis (else warn)
                            if dot(ex,ey) > 1e-3
                                warning('Note that local y-axis of element %i might be different than expected, because orien property is not normal to element axis.',eprops(i).elems(j));
                            end
                        end               
                    else
                        %check if default works
                        orien_def = [0 1 0];
                        for j=1:length(eprops(i).elems)
                            xp = nodes(elements(eprops(i).elems(j),1),:);
                            xq = nodes(elements(eprops(i).elems(j),2),:);
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %this directly from SPACAR 2015 source:
                            ex = xq - xp;
                            ex = ex/norm(ex);
                            ey_input = orien_def;
                            ey = ey_input(:)/norm(ey_input);
                            ex_proj = dot(ey,ex);
                            noemer = sqrt(1-ex_proj^2);
                            if noemer < 1e-5
                                err('No orien property specified for element %i. Default value [0 1 0] does not work, because (almost) parallel to element axis.',eprops(i).elems(j))
                            end
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        end  
                        eprops(i).orien = orien_def; %default setting
                    end

                    %check for mandatory fields when flex 
                    if (isfield(eprops(i),'flex') && ~isempty(eprops(i).flex))
                        validateattributes(eprops(i).flex,{'double'},{'vector','positive'},'',sprintf('flex property in eprops(%u)',i));  
                        
                        %Check if values for flex are valid
                        validateattributes(eprops(i).flex,{'double'},{'vector'},'',  sprintf('flex property in eprops(%u)',i));
                        if any(((eprops(i).flex==1)+(eprops(i).flex==2)+(eprops(i).flex==3)+(eprops(i).flex==4)+(eprops(i).flex==5)+(eprops(i).flex==6))==0)
                            err('Invalid deformation mode in eprops(%u).flex.',i)
                        end
                        
                        %Check if field exist in structure
                        if ~(isfield(eprops(i),'emod') && ~isempty(eprops(i).emod)); err('Property emod is not defined in eprops(%u)',i);     end
                        if ~(isfield(eprops(i),'smod') && ~isempty(eprops(i).smod)); err('Property smod is not defined in eprops(%u)',i);     end
                        if ~(isfield(eprops(i),'dens') && ~isempty(eprops(i).dens)); err('Property dens is not defined in eprops(%u)',i);     end
                        if ~(isfield(eprops(i),'type') && ~isempty(eprops(i).type)); err('Property type is not defined in eprops(%u)',i);     end
                        if ~(isfield(eprops(i),'dim') && ~isempty(eprops(i).dim));   err('Property dim is not defined in eprops(%u)',i);      end
                        if ~(isfield(eprops(i),'orien') && ~isempty(eprops(i).orien));   err('Property orien is not defined in eprops(%u)',i);      end

                    else
                        eprops(i).flex = [];
                        if (isfield(eprops(i),'emod') && ~isempty(eprops(i).emod)); warning('Property eprops(%u).emod is redundant without the flex property.',i);     end
                        if (isfield(eprops(i),'smod') && ~isempty(eprops(i).smod)); warning('Property eprops(%u).smod is redundant without the flex property.',i);     end
                    end

                end  
            end
            %warn user if elements are defined without adding properties to them
            el_without_prop_index = ~ismember(1:size(elements,1),el_nr_doubles_check);
            if any(el_without_prop_index)
                el_without_prop = find(el_without_prop_index);
                if length(el_without_prop) == 1
                    el_without_prop_str = num2str(el_without_prop);
                    warning('Element %s has no user-defined properties. Defaults (rigid massless elements) are used.',el_without_prop_str)
                else
                    el_without_prop_str = [num2str(el_without_prop(1)) sprintf(', %i',el_without_prop(2:end))];
                    warning('Elements %s have no user-defined properties. Defaults (rigid massless elements) are used.',el_without_prop_str)
                end
            end
            %warn user if no element has flexibility
            if ~isfield(eprops,'flex') || cellfun(@isempty,{eprops(:).flex})
                warning('No element seems to have the flex property. Simulation does not seem useful.')
            end
        end

        %CHECK RLS INPUT VARIABLE
        if exist('rls','var') && ~isempty(rls)
            ensure(all(ismember(fieldnames(rls),{'def'})),'Unknown field in rls; only def field is allowed.');
            
            for i=1:size(rls,2)
                if ~isempty(rls(i).def)
                    validateattributes(rls(i).def,{'double'},{'vector'},'',   sprintf('def property in rls(%u)',i));
                    if any(((rls(i).def==1)+(rls(i).def==2)+(rls(i).def==3)+(rls(i).def==4)+(rls(i).def==5)+(rls(i).def==6))==0)
                        err('Invalid deformation mode in rls(%u).',i)
                    end
                end
            end
        end

        %CHECK OPTIONAL ARGUMENTS
        if exist('opt','var')
            allowed_opts = {'filename','gravity','silent','buckload','showinputonly'};
            supplied_opts = fieldnames(opt);
            unknown_opts_i = ~ismember(supplied_opts,allowed_opts);
            if any(unknown_opts_i)
               err('Unknown opt field %s.',supplied_opts{unknown_opts_i});
            end
            
            if isfield(opt,'filename')
                if isempty(opt.filename)
                    warning('Filename cannot be empty. Filename spacar_file is used instead.');
                        opt.filename = 'spacar_file';
                else
                    validateattributes(opt.filename,{'char'},{'vector'},'',            'filename property in opt');  
                    if length(opt.filename) > 19
                        warning('Filename too long: maximum of 20 characters. Filename spacar_file is used instead.');
                        opt.filename = 'spacar_file';
                    end
                end
            end
            if (isfield(opt,'silent') && ~isempty(opt.silent))
                validateattributes(opt.silent,{'logical'},{'scalar'},'',            'silent property in opt');   end
            if (isfield(opt,'buckload') && ~isempty(opt.buckload))
                validateattributes(opt.buckload,{'logical'},{'scalar'},'',          'buckload property in opt'); end
            if (isfield(opt,'gravity') && ~isempty(opt.gravity))
                validateattributes(opt.gravity,{'double'},{'vector','numel',3},'',  'gravity property in opt');  end
            
        end
    end %END NOT-SILENT MODE BLOCK
    
    if ~fieldexist('opt','filename') || isempty(opt.filename)
        opt.filename = 'spacar_file';
    end
    
    % assign output 
    switch nargin
        case 1
            varargout{1} = nodes;
        case 2
            varargout{1} = nodes;
            varargout{2} = elements;
        case 3
            varargout{1} = nodes;
            varargout{2} = elements;
            varargout{3} = nprops;
        case 4
            varargout{1} = nodes;
            varargout{2} = elements;
            varargout{3} = nprops;
            varargout{4} = eprops;
        case 5
            varargout{1} = nodes;
            varargout{2} = elements;
            varargout{3} = nprops;
            varargout{4} = eprops;
            varargout{5} = rls;
        case 6
            varargout{1} = nodes;
            varargout{2} = elements;
            varargout{3} = nprops;
            varargout{4} = eprops;
            varargout{5} = rls;
            varargout{6} = opt;
    end

end

function showInput(varargin)
    disp('Showing input geometry');
end

function err(msg,varargin)
    %custom error function to hide the backtrace stuff in command window
    errorstruct.stack.file = '';
    errorstruct.stack.name = 'spacarlight';
    errorstruct.stack.line = 1;
    errorstruct.identifier = '';
 
    if nargin > 1
        errorstruct.message = sprintf(msg,varargin{:});
    else
        errorstruct.message = msg;
    end
    
    error(errorstruct)
    
end

function ensure(cond,msg,varargin)
    %custom assert function to hide the backtrace stuff in command window
    try
       %execute assert, but hide output
       assert(cond); 
    catch
        %if assert failed, throw error
        err(msg,varargin{:});
    end
    %if assert succeeded, do nothing
end

%% AUXILIARY FUNCTIONS
function stiffness = calc_stiffness(eprops)
% Compute the stiffness properties for leafspring or wireflexure
type    = eprops.type;
dim     = eprops.dim;
E       = eprops.emod;
G       = eprops.smod;
v       = E/(2*G) - 1;
switch lower(type)
    case {'leafspring','rigid'}
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



function inertia = calc_inertia(eprops)
% Compute the inertia properties for leafspring or wireflexure
type    = eprops.type;
dim     = eprops.dim;
rho     = eprops.dens;
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


function Results = calc_Results(filename, E_list, id_inputf, id_inputx, nodes, ~, ~, eprops, ~, opt)
nddof   = getfrsbf([filename '.sbd'],'nddof'); %number of dynamic DOFs
t_list  =  1:getfrsbf([filename,'.sbd'],'tdef'); %timesteps
lnp     = getfrsbf([filename,'.sbd'],'lnp'); %lnp data

if nddof == 0
    err('No dynamic degrees of freedom.')
end

%CHECK BUCKLING SETTINGS
calcbuck = false;
if (isfield(opt,'buckload') && opt.buckload == 1)
    calcbuck = true;
    if id_inputx
        warning('Input displacement prescribed; buckling load multipliers are also with respect to reaction forces due to this input.');
    end
    if ~id_inputf
        warning('No external forces are prescribed. Buckling values are not calculated.');
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

[propcrossect, Sig_nums]  = calc_propcrossect(E_list,eprops);
for i=t_list
    x       = getfrsbf([filename '.sbd'] ,'x', i);
    fxtot   = getfrsbf([filename '.sbd'] ,'fxt',i);
    for j=1:size(nodes,1)
        Results.step(i).node(j).x           = x(lnp((j-1)*2+1,1:3));
%         Results.step(i).node(j).rx_eulzyx   = quat2eul(x(lnp((j-1)*2+2,1:4))');
        Results.step(i).node(j).rx_quat     = (x(lnp((j-1)*2+2,1:4))');
        Results.step(i).node(j).Freac       = fxtot(lnp((j-1)*2+1,1:3)) ;
        
        Results.step(i).node(j).Mreac       = quat2eul(fxtot(lnp((j-1)*2+2,1:4))');
        [Results.step(i).node(j).CMglob, Results.step(i).node(j).CMloc]  =  complm(filename,(j-1)*2+1,(j-1)*2+2,i); %#ok<*AGROW>
    end
    [~,~,~,stressextrema] = stressbeam([filename,'.sbd'],Sig_nums,i,[],propcrossect);
    Results.step(i).stressmax = stressextrema.max*1e6;
    %  Results.step(i).bode_data =  getss('spacarfile',i);
end
Results.ndof = getfrsbf([filename '.sbd'] ,'ndof');

end


function cw= CWvalues(L,eprops)
w = eprops.dim(2);
E = eprops.emod;
G = eprops.smod;
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
    warning('complm() needs 3 or 4 input arguments');
    return;
end
if nargin < 4, tstp=0; end
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
        warning('Invalid node number for complm().');
        return;
    end
    if locv(i) <= nxp(1) || ...
            (locv(i)>(nxp(1)+nxp(2)) && locv(i) <= (nxp(1)+nxp(2)+nxp(3)))
        %  disp('WARNING: constrained node');
    end
end
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


function [propcrossect, Sig_nums]  = calc_propcrossect(E_list,eprops)
%restructure crossectional properties to evaluate stresses throuqh
%stressbeam.m

propcrossect = [];Sig_nums = [];

for i=1:size(E_list,1)
    id=[];
    for j=1:size(eprops,2)
        for k=1:length(eprops(j).elems)
            if eprops(j).elems(k)==i
                id=j;
            end
        end
    end
    if ~isempty(id)
        if (isfield(eprops(id),'flex') && ~isempty(eprops(id).flex))
            Elements = E_list(i,:);
            Elements(Elements==0) = [];
            Sig_nums = [Sig_nums Elements];
            
            switch eprops(id).type
                case 'leafspring'
                    for j=1:length(Elements)
                        propcrossect(end+1).CrossSection = 'rect';
                        propcrossect(end).Dimensions = [eprops(id).dim(1),eprops(id).dim(2)];
                    end
                case 'wire'
                    for j=1:length(Elements)
                        propcrossect(end+1).CrossSection = 'circ';
                        propcrossect(end).Dimensions = eprops(id).dim(1);
                    end
            end
        end
    end
end
end

function out = fieldexist(structname,fieldname)
    
    out = (exist(structname,'var') && isstruct(eval(structname)) && isfield(eval(structname),fieldname));
    
end