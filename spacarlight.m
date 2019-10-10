function results = spacarlight(varargin)
% SPACARLIGHT(nodes, elements, nprops, eprops, opt)
% runs SPACAR simulations with simplified syntax. See
% www.spacar.nl for more information.
%
% Created by: M. Nijenhuis and M. Naves
% Contact: m.naves@utwente.nl
%
% CONTRIBUTIONS
% J.P. Meijaard (complt)
% S.E. Boer (calc_stiffness, calc_inertia, calcTorsStiff and Spavisual functions)
% D.H. Wiersma (CWvalues)
%
% LIMITATIONS (note that the full Spacar version does allow these things)
% - Type of analysis: only static analyses are supported;
% - Boundary conditions: the orientation of a node can either be fixed or
% free. It is not possible to create a pinned boundary condition about a certain axis;
% - Prescribed motion: only displacements can be prescribed, no rotations.
%
% Output rotations are provided in quaternions, axis-angle representation and Euler angles (ZYX).
%
% For certain desired simulations, the current feature set of
% spacarlight() is too limited. In that case, the full version of SPACAR
% should be used. It offers *many* more features.
%
% Version 1.28
% 15-05-2019
sl_version = '1.28';

%% WARNINGS
warning off backtrace

%% CHECK FOR INCOMPLETE INPUT
%Temporary disable support for less then 4 input arguments (geometry visualization is not yet completed)
if nargin==0
    err('No input was provided.');
elseif nargin<4
    err('At least four input arguments are required.');
end

%initialize here since the following switch can already abort further execution
%and the results output needs to exist
results = struct();

switch nargin
    case 0
        err('No input was provided.');
    case 1
        warn('Incomplete input; no simulation is run.');
        % validate only nodes
        [nodes] = validateInput(varargin{:});
        return
    case 2
        warn('Incomplete input; no simulation is run.');
        % validate only nodes and elements
        [nodes,elements] = validateInput(varargin{:});
        return
    case 3
        warn('Incomplete input; no simulation is run.');
        % validate only nodes, elements and nprops
        [nodes,elements,nprops] = validateInput(varargin{:});
        return
    case 4
        % validate nodes, elements, nprops and eprops
        [nodes,elements,nprops,eprops] = validateInput(varargin{:});
        % attempt simulation
    case 5
        % validate all
        [nodes,elements,nprops,eprops,opt] = validateInput(varargin{:});
        if isfield(opt,'showinputonly') && opt.showinputonly == true
            showGeom(nodes,elements,[]);
            return
        end
        % attempt simulation
    otherwise
        err('Expecting a maximum of 5 input arguments.');
end


%% INITIALIZE VARIABLES, SET DEFAULTS (DO NOT SET DEFAULTS IN VALIDATEINPUT())
% note: if validateInput() gives a warning for a property in opt
% it will clear that value and return an empty field. So, **set the default here**
% and **do not set defaults in validateInput()**

%if no opt struct provided, create empty one
if ~(exist('opt','var') && isstruct(opt)); opt=struct(); end

%version number in opt (for further use in spacarlight) and in results (for output to user)
results.version = sl_version;

%set filename, even if not specified
if ~(isfield(opt,'filename') && ~isempty(opt.filename))
    opt.filename = 'spacar_file';
end

%determine whether silent mode
if ~(isfield(opt,'silent') && opt.silent == 1)
    opt.silent = false;
end

%determine whether to calculate compliance matrices
if ~isfield(opt,'calccompl')
    opt.calccompl = true;
end

%determine whether silent mode
if (isfield(opt,'transfer') && opt.transfer{1})
    opt.mode = 9;
elseif ~(isfield(opt,'mode'))
    opt.mode = 10;
elseif opt.mode==3
    opt.filename = [opt.filename '_3'];
else
    err('Unsupported simulation mode, use opt.mode = 10 or 3 instead.')
end

%determine whether to attempt autosolve
if ~isfield(opt,'rls')
    opt.autosolve = true;
    opt.rls = [];
else
    opt.autosolve = false;
end


%% CHECK EXISTENCE OF REQUIRED FUNCTIONS
ensure(exist('spacar','file') == 3,'spacar() is not in your path.');
ensure((exist('spavisual','file') == 2 || exist('spavisual','file') == 6),'spavisual() is not in your path.');
ensure((exist('stressbeam','file') == 2 || exist('stressbeam','file') == 6),'stressbeam() is not in your path (typically part of spavisual package).');
ensure((exist('getss','file') == 2 || exist('getss','file') == 6),'getss() is not in your path.');

%% BUILD DATFILE
[~, ~, E_list,~,~] = build_datfile(nodes,elements,nprops,eprops,opt,0);

%% SIMULATE FOR CHECKING CONSTRAINTS
try %try to run spacar in its silent mode
    warning('off','all')
    [~] = spacar(0,opt.filename);
    warning('on','all')
catch
    try %retry to run spacar in non-silent mode for old spacar versions
        warning('off','all')
        spacar(0,opt.filename);
        warning('on','all')
        warning('Old version of Spacar detected.')
        old_version = true; %#ok<NASGU> %Track if old version to prevent duplicate warning at mode 10 simulation
    catch msg
        switch msg.message
            case 'ERROR in subroutine PRPARE: Too many DOFs.'
                err('Too many degrees of freedom. Decrease the number of elements or the number of flexible deformations.');
            otherwise
                err(['Connectivity incorrect. Check element properties, node properties, element connectivity etc.\nCheck the last line of ' opt.filename '.log for more information.']);
        end
    end
end


%% CHECK CONSTRAINTS
[exactconstr, opt, overconstraints] = check_constraints(opt,E_list,eprops);
if ~exactconstr
    results.overconstraints = overconstraints;
    return;
end

%% RE-BUILD DATFILE
%appropriate releases should now be in opt.rls
[id_inputx, id_inputf, E_list, label_transfer_in, label_transfer_out] = build_datfile(nodes,elements,nprops,eprops,opt,opt.mode);

%% SIMULATE STATICS
try %run spacar in its silent mode
    warning('off','all')
    [~] = spacar(-opt.mode,opt.filename);
    warning('on','all')
    if ~(opt.silent)
        results.fighandle = spavisual(opt.filename);
        results.fighandle.Children.XLabel.String = 'x';
        results.fighandle.Children.YLabel.String = 'y';
        results.fighandle.Children.ZLabel.String = 'z';
        disp('Spacar simulation succeeded.')
    end
catch
    try %retry to run spacar in non-silent mode for old spacar versions
        warning('off','all')
        spacar(-opt.mode,opt.filename);
        warning('on','all')
        if ~(opt.silent)
            results.fighandle = spavisual(opt.filename);
            results.fighandle.Children.XLabel.String = 'x';
            results.fighandle.Children.YLabel.String = 'y';
            results.fighandle.Children.ZLabel.String = 'z';
        end
        disp('Spacar simulation succeeded.')
        if ~exist('old_version','var')
            warning('Old version of Spacar detected.')
        end
    catch msg
        %apparently, spacar mode 10 did not succeed.
        %try to figure out what went wrong:
        
        %1) see if a bigD matrix is available and whether its singular:
        try
            warning('off','all')
            [~] = spacar(0,opt.filename); %to get bigD
            warning('on','all')
            sbd     = [opt.filename '.sbd'];
            nep     = getfrsbf(sbd,'nep');
            nxp     = getfrsbf(sbd,'nxp');
            BigD    = getfrsbf(sbd,'bigd',1);
            Dcc     = BigD( 1:(nep(1)+nep(3)+nep(4)) , nxp(1)+(1:nxp(2)) );
            if (size(Dcc,1) ~= size(Dcc,2) || rank(Dcc) < size(Dcc,1) || det(Dcc) == 0)
                warn('Overconstraints could not be solved automatically; Try setting releases (opt.rls) manually.')
                %get the overconstraints in the system (without any autosolve attempt)
                %so build dat file again without any release attempts, do mode 0, get overconstraints
                if opt.autosolve
                    opt.autosolve = false;
                    opt.rls = [];
                    [~, ~, E_list] = build_datfile(nodes,elements,nprops,eprops,opt);
                    [~] = spacar(0,opt.filename);
                    [~, ~, overconstraints] = check_constraints(opt,E_list,eprops);
                    results.overconstraints = overconstraints;
                end
                return
            end
        catch msg
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        err('Spacar simulation failed. Possibly failed to converge to solution. Check magnitude of input displacements, loads, the number of loadsteps and other input data.')
    end
end
try
    %get results
    %note calc_results needs a results struct as input since it can already contain some fields
    results = calc_results(E_list, id_inputf, id_inputx, nodes, eprops, opt, label_transfer_in, label_transfer_out, results);
catch msg
    err(['A problem occurred processing simulation results. See ' opt.filename '.log for more information.'])
end

%% WARNINGS
warning backtrace on

% END OF SPACAR_LIGHT


    function [id_inputx, id_inputf, E_list, label_transfer_in, label_transfer_out] = build_datfile(nodes,elements,nprops,eprops,opt,mode)
        %returns E_list (amongst others):
        % spalight element i is represented by spacar beams E_list(i,:)
        
        %initialize values
        id_inputx   = false;                %identifier to check for prescribed input displacements/rotations
        id_inputf   = false;                %identifier to check for external load
        x_count     = size(nodes,1)*2+1;    %counter for node numbering
        e_count     = 1;                    %counter for element numbering
        X_list      = [];                   %list with node numbers
        E_list      = [];                   %list with element numbers
        
        %get username
        if ispc
            username = getenv('username');
        elseif ismac
            username = getenv('USER');
        end
        
        %% START CREATING DATFILE
        pr_I = sprintf('#Dat-file generated with SPACAR Light version %s\n#Date: %s\n#User: %s',sl_version,datestr(datetime),username);
        
        if mode==3
            pr_I = sprintf('%s \n%s',pr_I,'OUTLEVEL 0 1');
        end
        
        %% USERDEFINED NODES
        pr_N = sprintf('#NODES\t Nn\t\t\tX\t\t\tY\t\t\tZ');
        %print all nodes provides by user
        for i=1:size(nodes,1)
            pr_N = sprintf('%s\nX\t\t%3u\t\t\t%6f\t%6f\t%6f\t\t#node %u',pr_N,(i-1)*2+1,nodes(i,1),nodes(i,2),nodes(i,3),i);
        end
        
        
        %% ELEMENTS
        pr_E = sprintf('#ELEMENTS\t Ne\t\t Xp\t Rp\t Xq\t Rq\t\tOx\t\t\tOy\t\t\tOz');
        pr_D = sprintf('#DEF#\t\t Ne\t\t d1\t d2\t d3\t d4\t d5\t d6');
        
        for i=1:size(elements,1)
            pr_E = sprintf('%s\n#element %u',pr_E,i);
            
            N_p         = elements(i,1);            %p-node nodenumber
            N_q         = elements(i,2);            %q-node nodenumber
            X_list(i,1) = N_p;                      %#ok<*AGROW> %store p-node in X_list
            
            % get element property set corresponding to element i
            i_set = 0;
            for j=1:size(eprops,2)
                if any(eprops(j).elems==i)
                    %element i is in property set i_set%
                    i_set = j;
                    
                    %get element information
                    Flex        = eprops(j).flex;       %flexibility of this element
                    if isfield(eprops(j),'nbeams')
                        N = eprops(j).nbeams;           %number of beams per userdefined element
                    else
                        N = 1;
                    end
                    if isfield(eprops(j),'orien') && ~isempty(eprops(j).orien)
                        Orien = eprops(j).orien;        %orientation local y-vector
                    else
                        Orien = [0 1 0];
                    end
                    
                end
            end
            if i_set == 0 %defaults for if element does not exist in any element set
                N = 1;
                Orien = [0 1 0];
                Flex = [];
            end
            
            if N>1 %if more than 1 beam
                X_p = nodes(N_p,1:3);   %Location p-node
                X_q = nodes(N_q,1:3);   %Location q-node
                V   = X_q - X_p;        %Vector from p to q-node
                
                %create additional intermediate nodes
                for k = 1:N-1
                    
                    X = X_p+V/N*k;              %intermediate node position
                    pr_N = sprintf('%s\nX\t\t%3u\t\t\t%6f\t%6f\t%6f\t\t#intermediate node',pr_N,x_count,X(1),X(2),X(3));
                    X_list(i,k+1) = x_count;    %add intermediate node to X_list
                    
                    if k==1 %if the first beam, connect to p-node and first intermediate node
                        pr_E = sprintf('%s\nBEAM\t\t%3u\t\t%3u\t%3u\t%3u\t%3u\t\t%6f\t%6f\t%6f\t\t#beam %u',pr_E,e_count,(N_p-1)*2+1,(N_p-1)*2+2,x_count,x_count+1,Orien(1),Orien(2),Orien(3),k);
                    else    %if not the first beam, connect to two intermediate nodes
                        pr_E = sprintf('%s\nBEAM\t\t%3u\t\t%3u\t%3u\t%3u\t%3u\t\t%6f\t%6f\t%6f\t\t#beam %u',pr_E,e_count,x_count-2,x_count-1,x_count,x_count+1,Orien(1),Orien(2),Orien(3),k);
                    end
                    
                    if ~isempty(Flex)        %if element has flexibility, add dyne (no rlse, rlse is only added to last beam in element i)
                        pr_D = sprintf('%s\nDYNE\t\t%3u\t',pr_D,e_count);
                        for m=1:length(Flex) %loop over all flexible deformation modes
                            pr_D = sprintf('%s\t%3u',pr_D,Flex(m));
                        end
                    end
                    
                    E_list(i,k) = e_count;      %add beam number to E_list
                    e_count     = e_count+1;    %increase beam counter by 1
                    x_count     = x_count+2;    %increase node counter by 2 (+1 for rotation node)
                end
                
                %for the last beam in element i, connect to last intermediate node and q-node
                pr_E = sprintf('%s\nBEAM\t\t%3u\t\t%3u\t%3u\t%3u\t%3u\t\t%6f\t%6f\t%6f\t\t#beam %u',pr_E,e_count,x_count-2,x_count-1,(N_q-1)*2+1,(N_q-1)*2+2,Orien(1),Orien(2),Orien(3),k+1);
                
                X_list(i,k+2) = N_q;        %add q-node to X_list
                E_list(i,k+1) = e_count;    %add beam number to E_list
                
            else %if only a single beam is used, directly connect to p and q-node without intermediate noodes
                pr_E = sprintf('%s\nBEAM\t\t%3u\t\t%3u\t%3u\t%3u\t%3u\t\t%6f\t%6f\t%6f\t\t#beam %u',pr_E,e_count,(N_p-1)*2+1,(N_p-1)*2+2,(N_q-1)*2+1,(N_q-1)*2+2,Orien(1),Orien(2),Orien(3));
                
                X_list(i,2) = N_q;          %add q-node to X_list
                E_list(i,1) = e_count;      %add beam number to E_list
            end
            
            %for the last beam only, add dyne and/or rlse
            if ((~isfield(opt,'rls') || isempty(opt.rls)) && ~isempty(Flex)) %if no rlse, add all flexible deformation modes as dyne
                pr_D = sprintf('%s\nDYNE\t\t%3u\t',pr_D,e_count);
                for m=1:length(Flex)    %loop over all flexible deformation modes
                    pr_D = sprintf('%s\t%3u',pr_D,Flex(m));
                end
            else%if some rls are specified
                %compensate size of rls if size is smaller than element list
                if i>size(opt.rls,2)
                    opt.rls(i).def = [];
                end
                
                % add dyne
                if ~isempty(Flex)                           %if some flexibility is specified
                    dyn_added = false;                      %reset identifier to check if string 'dyne' is added
                    for m=1:length(Flex)                    %loop over all flexible deformation modes
                        if ~(sum(opt.rls(i).def==Flex(m))>0)   %if flexible deformation mode is not a rlse, it is dyne
                            if ~dyn_added                   %only add string 'dyne' if it is not yet added
                                pr_D = sprintf('%s\nDYNE\t\t%3u\t',pr_D,e_count);
                                dyn_added = true;           %set 'dyne' identifier
                            end
                            pr_D = sprintf('%s\t%3u',pr_D,Flex(m));
                        end
                    end
                end
                
                
                % add rlse
                rlse_added = false;                    %reset identifier to check if string 'rlse' is added
                for m=1:length(opt.rls(i).def)             %loop over all released deformation modes
                    if ~rlse_added                     %only add string 'rlse' if it is not yet added
                        pr_D = sprintf('%s\nRLSE\t\t%3u\t',pr_D,e_count);
                        rlse_added = true;
                    end
                    pr_D = sprintf('%s\t%3u',pr_D,opt.rls(i).def(m));
                end
            end
            
            e_count = e_count+1; %increase beam counter by 1 for last beam in the element
            x_count = x_count+2; %increase node counter by 2 (+1 for rotation node)
        end
        
        
        %% NODE FIXES AND INPUTS
        pr_fix = sprintf('#FIXES\t Nn');
        pr_input = sprintf('#INPUT\t Nn\t\tdir');
        
        if size(nodes,1)<size(nprops,2)
            err('Node properties applied to non-existing nodes.')
        end
        
        for i=1:size(nprops,2)
            %fixes
            if(isfield(nprops(i),'fix') && ~isempty(nprops(i).fix) &&  nprops(i).fix);            pr_fix = sprintf('%s\nFIX\t\t%3u  \nFIX\t\t%3u',pr_fix,(i-1)*2+1,(i-1)*2+2); end
            if(isfield(nprops(i),'fix_pos') && ~isempty(nprops(i).fix_pos) &&  nprops(i).fix_pos);     pr_fix = sprintf('%s\nFIX\t\t%3u',pr_fix,(i-1)*2+1);   end
            if(isfield(nprops(i),'fix_x') && ~isempty(nprops(i).fix_x) && nprops(i).fix_x);         pr_fix = sprintf('%s\nFIX\t\t%3u\t\t1',pr_fix,(i-1)*2+1);  end
            if(isfield(nprops(i),'fix_y') && ~isempty(nprops(i).fix_y) &&  nprops(i).fix_y);         pr_fix = sprintf('%s\nFIX\t\t%3u\t\t2',pr_fix,(i-1)*2+1);  end
            if(isfield(nprops(i),'fix_z') && ~isempty(nprops(i).fix_z) &&  nprops(i).fix_z);         pr_fix = sprintf('%s\nFIX\t\t%3u\t\t3',pr_fix,(i-1)*2+1);  end
            if(isfield(nprops(i),'fix_orien') && ~isempty(nprops(i).fix_orien) &&  nprops(i).fix_orien); pr_fix= sprintf('%s\nFIX\t\t%3u',pr_fix,(i-1)*2+2);   end
            
            %input displacements
            if (mode~=3 && mode~=9)
                if((isfield(nprops(i),'displ_x') && ~isempty(nprops(i).displ_x)) ||...
                        (isfield(nprops(i),'displ_initial_x') && ~isempty(nprops(i).displ_initial_x)));   pr_input = sprintf('%s\nINPUTX\t%3u\t\t1',pr_input,(i-1)*2+1);id_inputx = true;    end
                if((isfield(nprops(i),'displ_y') && ~isempty(nprops(i).displ_y)) ||...
                        (isfield(nprops(i),'displ_initial_y') && ~isempty(nprops(i).displ_initial_y)));   pr_input = sprintf('%s\nINPUTX\t%3u\t\t2',pr_input,(i-1)*2+1);id_inputx = true;    end
                if((isfield(nprops(i),'displ_z') && ~isempty(nprops(i).displ_z)) ||...
                        (isfield(nprops(i),'displ_initial_z') && ~isempty(nprops(i).displ_initial_z)));   pr_input = sprintf('%s\nINPUTX\t%3u\t\t3',pr_input,(i-1)*2+1);id_inputx = true;    end
                
                %input rotations
                id_inputr = 0; %identifier to count the number of input rotations
                if((isfield(nprops(i),'displ_rot') && ~isempty(nprops(i).displ_rot)) ||...
                        (isfield(nprops(i),'displ_initial_rot') && ~isempty(nprops(i).displ_initial_rot))); pr_input = sprintf('%s\nINPUTX\t%3u\t\t2 3 4',pr_input,(i-1)*2+2);id_inputx = true; id_inputr=id_inputr+1; end
                if((isfield(nprops(i),'rot_x') && ~isempty(nprops(i).rot_x)) ||...
                        (isfield(nprops(i),'rot_initial_x') && ~isempty(nprops(i).rot_initial_x))); pr_input = sprintf('%s\nINPUTX\t%3u\t\t2',pr_input,(i-1)*2+2);id_inputx = true; id_inputr=id_inputr+1; end
                if((isfield(nprops(i),'rot_y') && ~isempty(nprops(i).rot_y)) ||...
                        (isfield(nprops(i),'rot_initial_y') && ~isempty(nprops(i).rot_initial_y))); pr_input = sprintf('%s\nINPUTX\t%3u\t\t3',pr_input,(i-1)*2+2);id_inputx = true; id_inputr=id_inputr+1; end
                if((isfield(nprops(i),'rot_z') && ~isempty(nprops(i).rot_z)) ||...
                        (isfield(nprops(i),'rot_initial_z') && ~isempty(nprops(i).rot_initial_z))); pr_input = sprintf('%s\nINPUTX\t%3u\t\t4',pr_input,(i-1)*2+2);id_inputx = true; id_inputr=id_inputr+1; end
                if id_inputr>1 %if multiple rotations are prescribed, problems can arise with quaternion<->euler conversion
                    err('Multiple rotational inputs defined for node %u. Only a single input rotation can be added to a node.',i)
                end
                
            end
        end
        
        
        %% STIFFNESS/INERTIA PROPS
        pr_stiff = sprintf('#STIFFNESS\t Ne\tEA\t\t\t\t\t\t\tGJ\t\t\t\t\t\tEIy\t\t\t\t\t\tEIz\t\t\t\t\t\tShear Y\t\t\t\t\tShear Z');
        pr_mass = sprintf('#MASS\t\t Ne\t\t\tM/L\t\t\t\t\t\tJxx/L\t\t\t\t\tJyy/L\t\t\t\t\tJzz/L\t\t\t\t\tJyz/L');
        for i=1:size(eprops,2) %loop over each element property set
            for j=1:length(eprops(i).elems) %loop over all elements in element property set
                
                if (isfield(eprops(i),'dens') &&  ~isempty(eprops(i).dens))
                    inertia = calc_inertia(eprops(i));     %calculate mass properties
                    for k=1:size(E_list,2) %write mass/inertia values
                        El = E_list(eprops(i).elems(j),k); %loop over all beams in element set
                        if El>0
                            pr_mass = sprintf('%s\nEM\t\t\t%3u\t\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f',pr_mass,El,inertia(1),inertia(2),inertia(3),inertia(4),inertia(5));
                        end
                    end
                end
                %only write stiffness  when deformations are flexible
                if (isfield(eprops(i),'flex') && ~isempty(eprops(i).flex))
                    stiffness = calc_stiffness(eprops(i)); %calculate stiffness values
                    switch eprops(i).cshape
                        case 'rect'
                            L   = norm(nodes(elements(eprops(i).elems(j),2),:)...
                                - nodes(elements(eprops(i).elems(j),1),:));    %calculate flexure length for constrained warping values
                            
                            [cw_value, aspect] = cw_values(L,eprops(i));%calculate constrained warping values
                            
                            if (isfield(eprops(i),'cw') && ~isempty(eprops(i).cw))
                                if eprops(i).cw == 1
                                    cw = cw_value;
                                else
                                    cw = 1;
                                end
                            else
                                cw = 1;
                                if aspect < 3 && mode ~= 0
                                    warn('Aspect ratio of element %i is < 3; consider (constrained) warping.',eprops(i).elems(j));
                                end
                            end
                        case 'circ'
                            cw = 1;
                    end
                    
                    for k=1:size(E_list,2)                                  %loop over all beams in the element
                        El = E_list(eprops(i).elems(j),k);
                        if El>0
                            pr_stiff = sprintf('%s\nESTIFF\t\t%3u\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f',pr_stiff,El,stiffness(1),cw*stiffness(2),stiffness(3),stiffness(4),stiffness(5),stiffness(6));
                        end
                    end
                end
            end
        end
        
        
        %% FORCES/MOMENTS/NODAL MASSES
        pr_force =  sprintf('#FORCES\t\t Nn\t\tFx\t\t\tFy\t\t\tFz');
        pr_moment =  sprintf('#MOMENTS\t Nn\t\tMx\t\t\tMy\t\t\tMz');
        pr_dispr =  sprintf('#INPUTR\t\t Nn\t\tdir\t\tdr');
        pr_dispx =  sprintf('#INPUTX\t\t Nn\t\tdir\t\tdx');
        pr_transfer_in = sprintf('#TRANSFER IN\t  Index \tnode \tdir');
        pr_transfer_out = sprintf('#TRANSFER OUT\t  Index \tnode \tdir');
        label_transfer_in = []; label_transfer_out = [];
        tf_input_count = 1;
        tf_output_count = 1;
        pr_nm =  sprintf('#MASS\t\t  Nn\tM,Ixx\t\tIxy\t\t\tIxz\t\t\tIyy\t\t\tIyz\t\t\tIzz');
        id_ini  = false; %check for initial loading or displacement
        id_add  = false; %check for aditional loading or displacement
        
        for i=1:size(nprops,2) %loop over all user defined nodes
            if mode~=9;
                %forces
                if(isfield(nprops(i),'force') && ~isempty(nprops(i).force));                    pr_force = sprintf('%s\nDELXF\t\t%3u\t\t%6f\t%6f\t%6f',pr_force,(i-1)*2+1,nprops(i).force(1),nprops(i).force(2),nprops(i).force(3));                           id_add = true;  id_inputf=true; end
                if(isfield(nprops(i),'force_initial') && ~isempty(nprops(i).force_initial));    pr_force = sprintf('%s\nXF\t\t\t%3u\t\t%6f\t%6f\t%6f',pr_force,(i-1)*2+1,nprops(i).force_initial(1),nprops(i).force_initial(2),nprops(i).force_initial(3));   id_ini = true;  id_inputf=true; end
                
                %moments
                if(isfield(nprops(i),'moment') && ~isempty(nprops(i).moment)) %#ok<*ALIGN>
                    moments = nprops(i).moment;
                    pr_moment = sprintf('%s\nDELXF\t\t%3u\t\t%6f\t%6f\t%6f\t%6f',pr_moment,(i-1)*2+2,moments(1),moments(2),moments(3));                                       id_add = true;  id_inputf=true; end
                if(isfield(nprops(i),'moment_initial') && ~isempty(nprops(i).moment_initial))
                    moments_i = nprops(i).moment_initial;
                    pr_moment = sprintf('%s\nXF\t\t\t%3u\t\t%6f\t%6f\t%6f\t%6f',pr_moment,(i-1)*2+2,moments_i(1),moments_i(2),moments_i(3));                                 id_ini = true;  id_inputf=true; end
                
                %displacements
                if(isfield(nprops(i),'displ_x') && ~isempty(nprops(i).displ_x));                  pr_dispx = sprintf('%s\nDELINPX\t\t%3u\t\t1\t\t%6f',pr_dispx,(i-1)*2+1,nprops(i).displ_x(1));                       id_add = true; end
                if(isfield(nprops(i),'displ_y') && ~isempty(nprops(i).displ_y));                  pr_dispx = sprintf('%s\nDELINPX\t\t%3u\t\t2\t\t%6f',pr_dispx,(i-1)*2+1,nprops(i).displ_y(1));                       id_add = true; end
                if(isfield(nprops(i),'displ_z') && ~isempty(nprops(i).displ_z));                  pr_dispx = sprintf('%s\nDELINPX\t\t%3u\t\t3\t\t%6f',pr_dispx,(i-1)*2+1,nprops(i).displ_z(1));                       id_add = true; end
                if(isfield(nprops(i),'displ_initial_x') && ~isempty(nprops(i).displ_initial_x));  pr_dispx = sprintf('%s\nINPUTX\t\t%3u\t\t1\t\t%6f',pr_dispx,(i-1)*2+1,nodes(i,1) + nprops(i).displ_initial_x(1));  id_ini = true; end
                if(isfield(nprops(i),'displ_initial_y') && ~isempty(nprops(i).displ_initial_y));  pr_dispx = sprintf('%s\nINPUTX\t\t%3u\t\t2\t\t%6f',pr_dispx,(i-1)*2+1,nodes(i,2) + nprops(i).displ_initial_y(1));  id_ini = true; end
                if(isfield(nprops(i),'displ_initial_z') && ~isempty(nprops(i).displ_initial_z));  pr_dispx = sprintf('%s\nINPUTX\t\t%3u\t\t3\t\t%6f',pr_dispx,(i-1)*2+1,nodes(i,3) +nprops(i).displ_initial_z(1));   id_ini = true; end
                
                %axang rotations
                if(isfield(nprops(i),'displ_rot') && ~isempty(nprops(i).displ_rot)); rot = axang2quat(nprops(i).displ_rot);
                    pr_dispr = sprintf('%s\nDELINPX\t\t%3u\t\t2\t\t%6f',pr_dispr,(i-1)*2+2,rot(2)); id_add = true; 
                    pr_dispr = sprintf('%s\nDELINPX\t\t%3u\t\t3\t\t%6f',pr_dispr,(i-1)*2+2,rot(3));
                    pr_dispr = sprintf('%s\nDELINPX\t\t%3u\t\t4\t\t%6f',pr_dispr,(i-1)*2+2,rot(4)); end
                if(isfield(nprops(i),'displ_initial_rot') && ~isempty(nprops(i).displ_initial_rot));rot = axang2quat(nprops(i).displ_initial_rot);
                    pr_dispr = sprintf('%s\nINPUTX\t\t%3u\t\t2\t\t%6f',pr_dispr,(i-1)*2+2,rot(2));  id_ini = true;
                    pr_dispr = sprintf('%s\nINPUTX\t\t%3u\t\t3\t\t%6f',pr_dispr,(i-1)*2+2,rot(3));
                    pr_dispr = sprintf('%s\nINPUTX\t\t%3u\t\t4\t\t%6f',pr_dispr,(i-1)*2+2,rot(4)); end
                    
                %rotations
                if(isfield(nprops(i),'rot_x') && ~isempty(nprops(i).rot_x));                rot = eul2quat([0 0 nprops(i).rot_x(1)]);
                    pr_dispr = sprintf('%s\nDELINPX\t\t%3u\t\t2\t\t%6f',pr_dispr,(i-1)*2+2,rot(2)); id_add = true; end
                if(isfield(nprops(i),'rot_y') && ~isempty(nprops(i).rot_y));                rot = eul2quat([0 nprops(i).rot_y(1) 0]);
                    pr_dispr = sprintf('%s\nDELINPX\t\t%3u\t\t3\t\t%6f',pr_dispr,(i-1)*2+2,rot(3)); id_add = true; end
                if(isfield(nprops(i),'rot_z') && ~isempty(nprops(i).rot_z));                rot = eul2quat([nprops(i).rot_z(1) 0 0]);
                    pr_dispr = sprintf('%s\nDELINPX\t\t%3u\t\t4\t\t%6f',pr_dispr,(i-1)*2+2,rot(4)); id_add = true; end
                if(isfield(nprops(i),'rot_initial_x') && ~isempty(nprops(i).rot_initial_x));rot = eul2quat([0 0 nprops(i).rot_initial_x(1)]);
                    pr_dispr = sprintf('%s\nINPUTX\t\t%3u\t\t2\t\t%6f',pr_dispr,(i-1)*2+2,rot(2));  id_ini = true; end
                if(isfield(nprops(i),'rot_initial_y') && ~isempty(nprops(i).rot_initial_y));rot = eul2quat([0 nprops(i).rot_initial_y(1) 0]);
                    pr_dispr = sprintf('%s\nINPUTX\t\t%3u\t\t3\t\t%6f',pr_dispr,(i-1)*2+2,rot(3));  id_ini = true; end
                if(isfield(nprops(i),'rot_initial_z') && ~isempty(nprops(i).rot_initial_z));rot = eul2quat([nprops(i).rot_initial_z(1) 0 0]);
                    pr_dispr = sprintf('%s\nINPUTX\t\t%3u\t\t4\t\t%6f',pr_dispr,(i-1)*2+2,rot(4));  id_ini = true; end
            else
                if (isfield(nprops(i),'transfer_in') && ~isempty(nprops(i).transfer_in))
                    for j=1:length(nprops(i).transfer_in)
                        switch nprops(i).transfer_in{j}
                            case 'force_x'
                                pr_transfer_in = sprintf('%s\nINPUTF\t\t\t%3u\t\t  %3u \t  %3u',pr_transfer_in,tf_input_count,(i-1)*2+1,1);
                                label_transfer_in{tf_input_count} =  sprintf('force-x n%u',i);
                            case 'force_y'
                                pr_transfer_in = sprintf('%s\nINPUTF\t\t\t%3u\t\t  %3u \t  %3u',pr_transfer_in,tf_input_count,(i-1)*2+1,2);
                                label_transfer_in{tf_input_count} =  sprintf('force-y n%u',i);
                            case 'force_z'
                                pr_transfer_in = sprintf('%s\nINPUTF\t\t\t%3u\t\t  %3u \t  %3u',pr_transfer_in,tf_input_count,(i-1)*2+1,3);
                                label_transfer_in{tf_input_count} =  sprintf('force-z n%u',i);
                        end
                        tf_input_count = tf_input_count+1;
                    end
                end
                
                if (isfield(nprops(i),'transfer_out') && ~isempty(nprops(i).transfer_out))
                    for j=1:length(nprops(i).transfer_out)
                        switch nprops(i).transfer_out{j}
                            case 'displ_x'
                                pr_transfer_out = sprintf('%s\nOUTX\t\t\t%3u\t\t  %3u \t  %3u',pr_transfer_out,tf_output_count,(i-1)*2+1,1);
                                label_transfer_out{tf_output_count} =  sprintf('displ-x n%u',i);
                            case 'displ_y'
                                pr_transfer_out = sprintf('%s\nOUTX\t\t\t%3u\t\t  %3u \t  %3u',pr_transfer_out,tf_output_count,(i-1)*2+1,2);
                                label_transfer_out{tf_output_count} =  sprintf('displ-y n%u',i);
                            case 'displ_z'
                                pr_transfer_out = sprintf('%s\nOUTX\t\t\t%3u\t\t  %3u \t  %3u',pr_transfer_out,tf_output_count,(i-1)*2+1,3);
                                label_transfer_out{tf_output_count} =  sprintf('displ-z n%u',i);
                            case 'veloc_x'
                                pr_transfer_out = sprintf('%s\nOUTXP\t\t\t%3u\t\t  %3u \t  %3u',pr_transfer_out,tf_output_count,(i-1)*2+1,1);
                                label_transfer_out{tf_output_count} =  sprintf('veloc-x n%u',i);
                            case 'veloc_y'
                                pr_transfer_out = sprintf('%s\nOUTXP\t\t\t%3u\t\t  %3u \t  %3u',pr_transfer_out,tf_output_count,(i-1)*2+1,2);
                                label_transfer_out{tf_output_count} =  sprintf('veloc-y n%u',i);
                            case 'veloc_z'
                                pr_transfer_out = sprintf('%s\nOUTXP\t\t\t%3u\t\t  %3u \t  %3u',pr_transfer_out,tf_output_count,(i-1)*2+1,3);
                                label_transfer_out{tf_output_count} =  sprintf('veloc-z n%u',i);
                        end
                        tf_output_count = tf_output_count+1;
                    end
                end
                
                
            end
            
            %nodal masses/inertia
            if(isfield(nprops(i),'mass') && ~isempty(nprops(i).mass));                      pr_nm = sprintf('%s\nXM\t\t\t%3u\t\t%6f',pr_nm,(i-1)*2+1,nprops(i).mass); end
            if(isfield(nprops(i),'mominertia') && ~isempty(nprops(i).mominertia));          pr_nm = sprintf('%s\nXM\t\t\t%3u\t\t%6f\t%6f\t%6f\t%6f\t%6f\t%6f',pr_nm,(i-1)*2+2,nprops(i).mominertia(1),nprops(i).mominertia(2),nprops(i).mominertia(3),...
                    nprops(i).mominertia(4),nprops(i).mominertia(5),nprops(i).mominertia(6)); end
        end
        
        
        %% ADITIONAL OPTIONS
        %GRAVITY
        pr_add = sprintf('#ADDITIONAL PROPERTIES');
        if (exist('opt','var') && isstruct(opt) && isfield(opt,'gravity') && ~isempty(opt.gravity))
            pr_add = sprintf('%s\n#\t\t\tGx\t\t\tGy\t\t\tGz\nGRAVITY\t\t%6f\t%6f\t%6f',pr_add,opt.gravity(1),opt.gravity(2),opt.gravity(3));
        end
        
        %ITERSTEP SETTINGS
        if (exist('opt','var') && isstruct(opt) && isfield(opt,'loadsteps'))
            steps = opt.loadsteps-1;
        else
            steps = 9;
        end
        if      (id_ini && id_add);      pr_add = sprintf('%s\n\nITERSTEP\t10\t%3u\t0.0000005\t1\t3\t%3u',pr_add,steps,steps);    %if initial and aditional loading/displacement
        elseif  (id_ini && ~id_add);     pr_add = sprintf('%s\n\nITERSTEP\t10\t1\t0.0000005\t1\t1\t%3u',pr_add,steps);    %if initial loading/displacement
        elseif  (~id_ini && id_add);     pr_add = sprintf('%s\n\nITERSTEP\t10\t%3u\t0.0000005\t1\t3\t0',pr_add,steps);     %if initial loading/displacement
        else                             pr_add = sprintf('%s\n\nITERSTEP\t10\t1\t0.0000005\t1\t1\t0',pr_add);  end  %no loading/displacement
        
        
        %% VISUALIZATION
        pr_vis = sprintf('VISUALIZATION');
        
        %initial configuration color
        ini_color = [0 0 0];
        pr_vis = sprintf('%s\n\nINITIAL\nCOLOR\t\t%.2f\t%.2f\t%.2f',pr_vis,ini_color(1),ini_color(2),ini_color(3));
        
        for i=1:size(eprops,2) %loop over all element property sets
            %CROSSECTIONAL DIMENSIONS
            if (isfield(eprops(i),'cshape') && ~isempty(eprops(i).cshape))
                pr_vis = sprintf('%s\n\nBEAMPROPS',pr_vis);
                for j=1:length(eprops(i).elems) %loop over all elements in element set i
                    for k=1:size(E_list,2)
                        El = E_list(eprops(i).elems(j),k);
                        if El>0
                            pr_vis = sprintf('%s\t%3u',pr_vis,El);
                        end
                    end
                end
                
                switch eprops(i).cshape
                    case 'rect'
                        pr_vis = sprintf('%s\nCROSSTYPE\t  RECT',pr_vis);
                        pr_vis = sprintf('%s\nCROSSDIM\t  %6f\t%6f',pr_vis,eprops(i).dim(1),eprops(i).dim(2));
                    case 'circ'
                        pr_vis = sprintf('%s\nCROSSTYPE\t  CIRC',pr_vis);
                        pr_vis = sprintf('%s\nCROSSDIM\t  %f',pr_vis,eprops(i).dim(1));
                end
            end
            
            %COLOR
            %validateInput should make sure that the eprops.color field always exists for each set
            pr_vis = sprintf('%s\n\nGRAPHICS',pr_vis);
            for j=1:length(eprops(i).elems)
                for k=1:size(E_list,2)
                    El = E_list(eprops(i).elems(j),k);
                    if El>0
                        pr_vis = sprintf('%s\t%3u',pr_vis,El);
                    end
                end
            end
            pr_vis = sprintf('%s\nFACECOLOR\t  %3f %3f %3f',pr_vis,eprops(i).color(1),eprops(i).color(2),eprops(i).color(3));
            
            if isfield(eprops(i),'opacity') && ~isempty(eprops(i).opacity)
                pr_vis = sprintf('%s\nOPACITY\t%.2f',pr_vis,eprops(i).opacity);
            end
            
            %HIDE
            if (isfield(eprops(i),'hide') && ~isempty(eprops(i).hide) && eprops(i).hide==1)
                pr_vis = sprintf('%s\n\nDONOTDRAW',pr_vis);
                for j=1:length(eprops(i).elems)
                    for k=1:size(E_list,2)
                        El = E_list(eprops(i).elems(j),k);
                        if El>0
                            pr_vis = sprintf('%s\t%u',pr_vis,El);
                        end
                    end
                end
            end
        end
        
        %in case an element is not in a set, (still) give it a color
        element_without_set = [];
        no_set_color = [206 109 116]/255;
        for i=1:size(elements,1)
            if isempty(find(arrayfun(@(x) ismember(i,x.elems),eprops),1))
                element_without_set(end+1) = i;
            end
        end
        if any(element_without_set)
            pr_vis = sprintf('%s\n\nGRAPHICS',pr_vis);
            for i=1:length(element_without_set)
                pr_vis = sprintf('%s\t%i',pr_vis,E_list(element_without_set(i),1));
            end
            pr_vis = sprintf('%s\nFACECOLOR\t  %3f %3f %3f',pr_vis,no_set_color(1),no_set_color(2),no_set_color(3));
        end
        try
            if(isfield(opt,'customvis') && ~isempty(opt.customvis));
                for i=1:size(opt.customvis,2)
                    pr_vis = sprintf('%s\n%s',pr_vis,opt.customvis{i});
                end
            end
        catch msg
            err('Invalid opt.customvis input')
        end
        
        fileID = fopen([opt.filename '.dat'],'w');
        print_dat(fileID,'%s\n\n\n\n',pr_I);
        print_dat(fileID,'%s\n\n\n\n',pr_N);
        print_dat(fileID,'%s\n\n\n\n',pr_E);
        print_dat(fileID,'%s\n\n\n\n',pr_D);
        print_dat(fileID,'%s\n\n\n\n',pr_fix);
        print_dat(fileID,'%s\n\n\n',pr_input);
        fprintf(fileID,'\nEND\nHALT\n\n\n');
        print_dat(fileID,'%s\n\n\n\n',pr_stiff);
        print_dat(fileID,'%s\n\n\n\n',pr_mass);
        print_dat(fileID,'%s\n\n\n\n',pr_nm);
        print_dat(fileID,'%s\n\n\n\n',pr_add);
        if mode~=9
            print_dat(fileID,'%s\n\n\n\n',pr_force);
            print_dat(fileID,'%s\n\n\n\n',pr_moment);
            print_dat(fileID,'%s\n\n\n\n',pr_dispr);
            print_dat(fileID,'%s\n\n\n\n',pr_dispx);
        else
            fprintf(fileID,'END\nHALT\n\n\n\n');
            print_dat(fileID,'%s\n\n\n\n',pr_transfer_in);
            print_dat(fileID,'%s\n\n\n\n',pr_transfer_out);
        end
        fprintf(fileID,'END\nEND\n\n\n\n');
        print_dat(fileID,'%s',pr_vis);
        fclose(fileID); %datfile finished!
    end

    function [exactconstr, opt, overconstraints] = check_constraints(opt,E_list,eprops)
        overconstraints = []; %initialize in order to not fail the function varout check
        %in some cases results structure will be filled by
        %a list of overconstraints or a list of partial release solutions
        
        %load data
        sbd     = [opt.filename '.sbd'];
        nep     = getfrsbf(sbd,'nep');
        nxp     = getfrsbf(sbd,'nxp');
        %nddof   = getfrsbf(sbd,'nddof');
        le      = getfrsbf(sbd,'le');
        BigD    = getfrsbf(sbd,'bigd',1);
        Dcc     = BigD( 1:(nep(1)+nep(3)+nep(4)) , nxp(1)+(1:nxp(2)) );
        IDlist = 1:size(Dcc,1);
        [ U, s, V ] = svd(Dcc);
        s       = diag(s);
        
        %if no degrees of freedom
        %if nddof == 0
        %    warn('The system has no degrees of freedom (so no Spacar simulation will be performed). Check eprops.flex and rls.')
        %    exactconstr = false; %(false here so that main function will not proceed but abort instead)
        %    return
        %end
        
        %if empty s
        if isempty(s) %%% TO BE DONE: s can be empty. What does this mean? => no calculable nodes? no freedom?
            exactconstr = false;
            return
        end
        
        % calcualate number of over/under constraints
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
            warn('System is underconstrained. Check element connectivity, boundary conditions and releases.')
            exactconstr = false;
            return
        elseif nover>0 %overconstrained
            
            if ~opt.autosolve %do not autosolve, but calculate overconstraints
                
                %part of U matrix coresponding to all overconstraints
                overconstraint = U(:,end-nover+1:end);
                %select part of U matrix corresponding to first overconstraint
                oc = overconstraint(:,1);
                dofs = (1:length(oc))';
                [oc_sorted,order] = sort(oc.^2,1,'descend');
                dofs_sorted = dofs(order);
                
                %select only part that explains 99% (or more) of singular vector's length
                idx = find(cumsum(oc_sorted)>sqrt(1-1e-6),1,'first');
                %overconstrained dofs (all deformations)
                overconstr_dofs = dofs_sorted(1:idx); %=sel
                overconstr_singvec_part = oc_sorted(1:idx);
                
                %create a list with [1] element number,[2] deformation mode and [3] oc value
                listData = zeros(numel(overconstr_dofs),3);
                listData(:,3) = overconstr_singvec_part;
                for i=1:size(listData,1)
                    [elnr,defpar] = find(le==overconstr_dofs(i));
                    listData(i,1:2) = [elnr, defpar];
                end
                
                %now listData is in terms of spacar beams, not spalight elements.
                %reshape rls suggestions according to user defined elements,
                %since spacar might return a lot more release options for each actual beam element
                %(and spacar light hides the subdivision elements from the user)
                OC_el= [];
                OC_defs = [];
                
                %(spalight element i is represented by spacar beams E_list(i,:))
                for k=1:size(E_list,1) %kth spacarlight element
                    list = [];
                    for j=1:size(E_list,2) %loop over all spacar beams for kth spalight element
                        %find all occurences of spacar beams in listData,
                        %get corresponding deformation modes
                        list = [list; sort(listData(find(listData(:,1)==E_list(k,j)),2))]; %#ok<FNDSB>
                    end
                    %list contains overconstrained deformations of all beams in kth spalight element
                    %(could very well contain multiple 1's, 2's ... 6's .
                    
                    red_list=[]; %reduce deformations to release each deformation once at most
                    for j=1:6 %loop over deformation modes
                        if ismember(j,list) %if mode in overconstrained list
                            %and not already in opt.rls (which could be empty)
                            if size(opt.rls,2)==0 || (size(opt.rls,2) < k || ~ismember(j,opt.rls(k).def))
                                % then add to reduction list
                                red_list(end+1) = j;
                            end
                        end
                    end
                    %store deformations of kth element (red_list) in collection of
                    %defs en els
                    if ~isempty(red_list)
                        OC_defs(end+1,1:length(red_list)) = red_list;
                        OC_el(end+1,1) = k;
                    end
                end
                
                overconstraints = [OC_el OC_defs];
                warn('System is overconstrained; releases are required in order to run a static simulation.\nA suggestion for possible releases is given in the table below.\n')
                fprintf('Number of overconstraints: %u\n\n',nover);
                disp(table(OC_el,sum((OC_defs==1),2),sum((OC_defs==2),2),sum((OC_defs==3),2),sum((OC_defs==4),2),sum((OC_defs==5),2),sum((OC_defs==6),2),...
                    'VariableNames',{'Element' 'def_1' 'def_2 ' 'def_3' 'def_4' 'def_5' 'def_6'}));
                exactconstr = false;
                return
                
            else %autosolve overconstrains
                %create empty rlse matrix used internally for autosolver
                rlse = zeros(size(E_list,1),6);
                %number of initial overconstraints
                n = nover;
                
                for dummy = 1:n %solve for n overconstraints, loop counter is not required (dummy)
                    %redo SVD analyses (has to be redone after each resolved overconstrain)
                    [ U, ~, ~ ] = svd(Dcc);
                    
                    %part of U matrix coresponding to all overconstraints
                    overconstraint = U(:,end-nover+1:end);
                    %select part of U matrix corresponding to first overconstraint
                    oc = overconstraint(:,1);
                    [oc_sort,order] = sort(oc.^2,1,'descend');
                    % select only part that explains 99% (or more) of singular vector's length
                    idx = find(cumsum(oc_sort)>sqrt(1-1e-6),1,'first');
                    idx = order(1:idx);
                    sel = (1:numel(oc))';
                    %overconstrained deformation modes
                    sel = sel(idx);
                    
                    loop = true; %loop until proper solution for overconstraints is obtained
                    j = 1;       %start with deformation mode with highest singular value
                    while loop
                        if j>length(sel) %unsolvable, all deformations modes have been attempted, probably rigid body release required
                            warn('Overconstraints could not be solved automatically; partial release information is provided in the output.')
                            fprintf('Number of overconstraints left: %u\n\n',nover);
                            disp( table((1:size(rlse,1))',rlse(:,1),rlse(:,2),rlse(:,3),rlse(:,4),rlse(:,5),rlse(:,6),'VariableNames',{'Element' 'def_1' 'def_2 ' 'def_3' 'def_4' 'def_5' 'def_6'}))
                            %                     results.rls = restruct_rlse(rlse);
                            %workspace output of //remaining// overconstraints in same style as
                            %//full// overconstraint workspace output
                            for v = 1:size(rlse,1)
                                rlsout(v,1:nnz(rlse(v,:))+1) = [v find(rlse(v,:)==1)];
                            end
                            overconstraints = rlsout;
                            exactconstr = false;
                            return
                        end
                        
                        %Select deformation to be test for release
                        def_id = sel(j);
                        def = IDlist(def_id);
                        %After each release, Dcc matrix is
                        %reduced in size, causing shifting
                        %of deformation modes.
                        %i.e., if deformation mode 2 is
                        %released, deformation 3 shifts to
                        %position 2, deformation 4 to 3,
                        %etc.
                        %The IDlist variable tracks and
                        %compensates for this shift
                        
                        %element number and deformation
                        [elnr,defmode] = find(le==def);
                        
                        for k=1:size(E_list,1) %kth spacar light element
                            if elnr==E_list(k,find(E_list(k,:)~=0,1,'last')) %if element number is in last beam of element k
                                
                                %Check if kth spacar light element is flexible
                                for i=1:size(eprops,2) %loop over all element property sets to check flexibility
                                    if (isfield(eprops(i),'flex') && ~isempty(eprops(i).flex) && any(eprops(i).elems==k)) %if kth element is flexible
                                        if rlse(k,defmode) == 0 %if this deformation is not yet released
                                            IDlist(def_id) = []; %update IDlist for deformation tracking
                                            Dcc(def_id,:) = []; %reduce Dcc matrix to remove this overconstrained
                                            rlse(k,defmode) = 1; %store release informatino in rlse matrix ()
                                            nover = nover-1; %reduce number of overconstraints
                                            loop = false; %terminate while loop and redo from line 235
                                        end
                                    end
                                end
                            end
                        end
                        j=j+1;
                    end
                end
                %release information is in matrix form (rlse), it has to be restructured
                %for ouput in structure
                opt.rls = restruct_rlse(rlse);
                exactconstr = true;
            end
        else
            %no underconstraints and no overconstraints
            exactconstr = true;
        end
    end

    function varargout = validateInput(varargin)
        %this function receives the user-supplied input and only returns
        %that input when it turns out to be valid
        %**** do not set defaults in this function, ****
        %**** but in the main file at the designated place ****%
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
                opt = varargin{5};
        end
        
        %BEGIN NOT-SILENT MODE BLOCK
        if ~(exist('opt','var') && isstruct(opt) && isfield(opt,'silent') && opt.silent==1) %checks are skipped in silent mode
            
            %CHECK NODES INPUT VARIABLE
            validateattributes(nodes,   {'double'},{'ncols',3,'ndims',2},'','nodes')
            
            %CHECK ELEMENTS INPUT VARIABLE
            if exist('elements','var')
                nno = size(nodes,1);
                validateattributes(elements,{'double'},{'ncols',2,'ndims',2},'','elements')
                ensure(all(elements(:) == floor(elements(:))),'Entries in elements seem to be non-integers.')
                
                ensure(all(elements(:)>0),'Element seems connected to node number <=0.')
                ensure(max(elements(:))<=nno,'Element seems connected to node that does not exist.')
                if max(elements(:))<nno; warn('Node seems not attached to element.'); end
                ensure(~any(abs(elements(:,1)-elements(:,2))==0),('Both sides of element seem connected to the same node.'))
                
                %check if unique pairs node numbers (independent of p/q order)
                ensure(size(unique(sort(elements,2),'rows'),1)==size(elements,1),'Multiple elements seem connected between the same node pair.')
                
                ensure_idelret(sqrt(sum((nodes(elements(:,1),:) - nodes(elements(:,2),:)).^2,2))>1e-5,'length seems smaller than 0.00001.')
                
                maxlength = max(sqrt(sum((nodes(elements(:,1),:) - nodes(elements(:,2),:)).^2,2)));
                minlength = min(sqrt(sum((nodes(elements(:,1),:) - nodes(elements(:,2),:)).^2,2)));
                ensure(maxlength/minlength<=1000,'Ratio between element lengths seems larger than 1000.');
                
            end
            
            %CHECK NPROPS INPUT VARIABLE
            if exist('nprops','var')
                
                allowed_nprops = {'fix','fix_x','fix_y','fix_z','fix_pos','fix_orien','displ_x','displ_y','displ_z','displ_rot','rot_x','rot_y','rot_z','force','moment','mass','mominertia','force_initial','moment_initial', ...
                    'displ_initial_x','displ_initial_y','displ_initial_z','displ_initial_rot','rot_initial_x','rot_initial_y','rot_initial_z','transfer_in','transfer_out'};
                supplied_nprops = fieldnames(nprops);
                ensure(size(supplied_nprops,1)>0,'Node properties seem empty.')
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
                            case {'displ_rot','displ_initial_rot'}
                                if ~isempty(nprops(i).(Node_fields{j}));     validateattributes(nprops(i).(Node_fields{j}),{'double'},{'vector','numel',4},'',  sprintf('displ_rot property in nprops(%u)',i));      end
                            case 'mass'
                                if ~isempty(nprops(i).(Node_fields{j}));     validateattributes(nprops(i).(Node_fields{j}),{'double'},{'scalar'},'',             sprintf('mass property in nprops(%u)',i));      end
                            case 'mominertia'
                                if ~isempty(nprops(i).(Node_fields{j}));     validateattributes(nprops(i).(Node_fields{j}),{'double'},{'vector','numel',6},'',   sprintf('mominertia property in nprops(%u)',i));   end
                                
                            case 'transfer_in'
                                if ~isempty(nprops(i).(Node_fields{j}))
                                    validateattributes(nprops(i).(Node_fields{j}),{'char','cell'},{'nonempty'},'',   sprintf('transfer_in property in nprops(%u)',i));
                                    %if single string (so char), convert to a cell
                                    if ischar(nprops(i).transfer_in)
                                        nprops(i).transfer_in = cellstr(nprops(i).transfer_in);
                                    end
                                    %since transfer_in field can be a cell of multiple char entries, loop over them and check if all of them are allowed
                                    for k=1:length(nprops(i).transfer_in)
                                        ensure(any(strcmp(nprops(i).transfer_in{k},{'force_x','force_y','force_z'})),'Unknown transfer_in{%u} value in nprops(%u)',k,i);
                                    end
                                    %check if no double entries
                                    [~,ia,ib] = unique(nprops(i).transfer_in);
                                    ensure(length(ia)==length(ib),'nprops(%u).transfer_in should contain unique values.',i);
                                end
                            case 'transfer_out'
                                if ~isempty(nprops(i).(Node_fields{j}))
                                    validateattributes(nprops(i).(Node_fields{j}),{'char','cell'},{'nonempty'},'',   sprintf('transfer_out property in nprops(%u)',i));
                                    %if single string (so char), convert to a cell
                                    if ischar(nprops(i).transfer_out)
                                        nprops(i).transfer_out = cellstr(nprops(i).transfer_out);
                                    end
                                    %since transfer_out field can be a cell of multiple char entries, loop over them and check if they are allowed
                                    for k=1:length(nprops(i).transfer_out)
                                        ensure(any(strcmp(nprops(i).transfer_out{k},{'displ_x','displ_y','displ_z','veloc_x','veloc_y','veloc_z'})),'Unknown transfer_out{%u} value in nprops(%u)',k,i);
                                    end
                                    %check if no double entries
                                    [~,ia,ib] = unique(nprops(i).transfer_out);
                                    ensure(length(ia)==length(ib),'nprops(%u).transfer_out should contain unique values.',i);
                                end
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
                    if((isfield(nprops(i),'displ_rot') && ~isempty(nprops(i).displ_rot)) ||...
                            (isfield(nprops(i),'displ_initial_rot') && ~isempty(nprops(i).displ_initial_rot))); count_bcs = count_bcs + 3;   end
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
                        ((isfield(nprops(i),'rot_x') && ~isempty(nprops(i).rot_x)) || (isfield(nprops(i),'rot_y') && ~isempty(nprops(i).rot_y)) || (isfield(nprops(i),'rot_z') && ~isempty(nprops(i).rot_z)) || (isfield(nprops(i),'displ_rot') && ~isempty(nprops(i).displ_rot)) || (isfield(nprops(i),'displ_initial_rot') && ~isempty(nprops(i).displ_initial_rot)) || (isfield(nprops(i),'rot_initial_x') && ~isempty(nprops(i).rot_initial_x)) || (isfield(nprops(i),'rot_initial_y') && ~isempty(nprops(i).rot_initial_y)) || (isfield(nprops(i),'rot_initial_z') && ~isempty(nprops(i).rot_initial_z)) ) ...
                        ])<=1,'There is a combination of fix_orien, moment and rot_x/y/z on node %i.',i);
                    
                    %no combination of fix (6 constraints) and (mass or mominertia)
                    if (isfield(nprops(i),'fix') && ~isempty(nprops(i).fix) && nprops(i).fix == true) && ...
                            (   (isfield(nprops(i),'mass') && ~isempty(nprops(i).mass) && nprops(i).mass~=0) || ...
                            (isfield(nprops(i),'mominertia') && ~isempty(nprops(i).mominertia) && any(nprops(i).mominertia~=0)) ...
                            )
                        
                        warn('Inertia associated with fixed node %i.',i);
                    end
                    
                    %no combination of fix_pos and mass
                    if (isfield(nprops(i),'fix_pos') && ~isempty(nprops(i).fix_pos) && nprops(i).fix_pos == true) && ...
                            (isfield(nprops(i),'mass') && ~isempty(nprops(i).mass) && nprops(i).mass~=0)
                        
                        warn('Mass associated with position-fixed node %i.',i);
                    end
                    
                    %no combination of fix_orien and mominertia
                    if (isfield(nprops(i),'fix_orien') && ~isempty(nprops(i).fix_orien) && nprops(i).fix_orien == true) && ...
                            (isfield(nprops(i),'mominertia') && ~isempty(nprops(i).mominertia) && any(nprops(i).mominertia~=0))
                        
                        warn('Moment of inertia associated with orientation-fixed node %i.',i);
                    end
                    
                    %checks in case of transfer function input or output
                    if (isfield(nprops(i),'transfer_in') && ~isempty(nprops(i).transfer_in)) || (isfield(nprops(i),'transfer_out') && ~isempty(nprops(i).transfer_out))
                        %no combination of (fix or fix_pos) and any transfer_in/transfer_out argument
                        if (isfield(nprops(i),'fix') && any(nprops(i).fix==1)) || (isfield(nprops(i),'fix_pos') && any(nprops(i).fix_pos==1))
                            err('Cannot use fixed node %u as transfer function input or output.',i);
                        end
                        
                        %no combination of fix_x and (force_x, displ_x, veloc_x)
                        if isfield(nprops(i),'fix_x') && any(nprops(i).fix_x==1) && (...
                                (isfield(nprops(i),'transfer_in') && any(strcmp('force_x',nprops(i).transfer_in))) || ...
                                (isfield(nprops(i),'transfer_out') && (any(strcmp('displ_x',nprops(i).transfer_out)) || any(strcmp('veloc_x',nprops(i).transfer_out)))) ...
                                )
                            err('Cannot use fixed node %u as transfer function input or output in x-direction',i);
                        end
                        %no combination of fix_y and (force_y, displ_y, veloc_y)
                        if isfield(nprops(i),'fix_y') && any(nprops(i).fix_y==1) && (...
                                (isfield(nprops(i),'transfer_in') && any(strcmp('force_y',nprops(i).transfer_in))) || ...
                                (isfield(nprops(i),'transfer_out') && (any(strcmp('displ_y',nprops(i).transfer_out)) || any(strcmp('veloc_y',nprops(i).transfer_out)))) ...
                                )
                            err('Cannot use fixed node %u as transfer function input or output in y-direction',i);
                        end
                        %no combination of fix_z and (force_z, displ_z, veloc_z)
                        if isfield(nprops(i),'fix_z') && any(nprops(i).fix_z==1) && (...
                                (isfield(nprops(i),'transfer_in') && any(strcmp('force_z',nprops(i).transfer_in))) || ...
                                (isfield(nprops(i),'transfer_out') && (any(strcmp('displ_z',nprops(i).transfer_out)) || any(strcmp('veloc_z',nprops(i).transfer_out)))) ...
                                )
                            err('Cannot use fixed node %u as transfer function input or output in z-direction',i);
                        end
                    end
                    
                end
                ensure(count_bcs >= 6,'The nodes seem to have insufficient (%i<6) constraints (fix, displ, or rot).',count_bcs);
            end
            
            %CHECK EPROPS INPUT VARIABLE
            if exist('eprops','var')
                allowed_eprops = {'elems','emod','smod','dens','cshape','dim','orien','nbeams','flex','color','hide','opacity','cw'};
                supplied_eprops = fieldnames(eprops);
                unknown_eprops_i = ~ismember(supplied_eprops,allowed_eprops);
                if any(unknown_eprops_i)
                    err('Unknown eprops field %s.',supplied_eprops{unknown_eprops_i});
                end
                
                ignoresets = []; %sets without elems property (e.g. empty sets) are removed after the for loop
                %note: this means that set numbers might have changed after this deletion!
                el_nr_doubles_check = []; %filling this with user defined element numbers (with .elems) to check for doubles and missing elements
                for i=1:size(eprops,2)
                    %check for elems field, the only mandatory field
                    if ~(isfield(eprops(i),'elems') && ~isempty(eprops(i).elems))
                        warn('Property elems is not defined in eprops(%u); ignoring eprops(%i).',i,i);
                        ignoresets(end+1) = i;
                    else
                        %elems exists (checked by previous for loop), so validate elems
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
                        if (isfield(eprops(i),'dim') && ~isempty(eprops(i).dim))
                            validateattributes(eprops(i).dim,{'double'},{'vector'},'',sprintf('dim property in eprops(%u)',i));
                            if ~(isfield(eprops(i),'cshape') && ~isempty(eprops(i).cshape))
                                warn('Property eprops(%u).dim is redundant without the cshape property.',i)
                            end
                        end
                        
                        %color field
                        %approach: color can be specified in various ways;
                        %always convert to rgb values between 0-1
                        colorset = struct(...
                            'name',{'blue','grey','darkblue','darkgrey','green','darkgreen'},...
                            'value',{   [162 195 214]/255,...
                            [198 198 198]/255,...
                            [0 124 176]/255,...
                            [112 112 112]/255,...
                            [128 194 143]/255,...
                            [90 137 101]/255});
                        
                        if (isfield(eprops(i),'color') && ~isempty(eprops(i).color))
                            if isfloat(eprops(i).color) %also returns true for "integers" like 2 or 2.0
                                %check for integer, or 3-vec, and validate
                                if size(eprops(i).color,2) == 1 && size(eprops(i).color,1) == 1 %integer
                                    ensure(mod(eprops(i).color,1)==0 && eprops(i).color>0,'Used as a scalar, value of eprops(%i).color should be a positive integer.',i);
                                    ensure(eprops(i).color<=size(colorset,2),'Pre-defined color with index eprops(%i).color does not exist.',i)
                                    eprops(i).color = colorset(eprops(i).color).value;
                                elseif size(eprops(i).color,2) == 3 && size(eprops(i).color,1) == 1 %1x3 vector
                                    ensure(all(eprops(i).color>=0) && all(eprops(i).color<=255),'Used as a 1x3 vector, values of eprops(%i).color should be >=0 and <=255.',i);
                                    if any(eprops(i).color>1.0) %some values are >1.0
                                        %user uses rgb between 0-255 notation
                                        ensure(all(mod(eprops(i).color,1)==0),'Used as a 1x3 vector with value(s) larger than 1.0, all values of eprops(%i).color should be integers.',i);
                                        %divide by 255
                                        eprops(i).color = eprops(i).color/255;
                                    end
                                else
                                    err('Property eprops(%i).color should be a string, an integer, or a 1x3 vector.',i);
                                end
                            elseif ischar(eprops(i).color)
                                ci = arrayfun(@(x) strcmp(x.name,eprops(i).color),colorset);
                                ensure(any(ci),'Pre-defined color with name eprops(%i).color does not exist.',i);
                                eprops(i).color = colorset(ci).value;
                            else
                                err('Property eprops(%i).color should be a string, an integer, or a 1x3 vector.',i);
                            end
                        else
                            %this is setting a default.. do we want that here?
                            ii = mod(i-1,size(colorset,2))+1; %modulo with number of predefined colors
                            %(using the -1 and +1 because modulo otherwise returns 0's)
                            eprops(i).color = colorset(ii).value;
                        end
                        
                        if (isfield(eprops(i),'hide') && ~isempty(eprops(i).hide));     validateattributes(eprops(i).hide,{'logical'},{'scalar'},'',sprintf('hide property in eprops(%u)',i)); end
                        if (isfield(eprops(i),'opacity') && ~isempty(eprops(i).opacity)); validateattributes(eprops(i).opacity,{'double'},{'scalar','>',0,'<',1},'',   sprintf('opacity property in eprops(%u)',i)); end
                        
                        if (isfield(eprops(i),'cshape') && ~isempty(eprops(i).cshape))
                            validateattributes(eprops(i).cshape,{'char'},{'nonempty'},'',sprintf('cshape property in eprops(%u)',i));
                            if ~any(strcmp(eprops(i).cshape,{'rect','circ'})), err('Element cshape should be either rect or circ.'); end
                            switch eprops(i).cshape
                                case 'rect'
                                    if ~(isfield(eprops(i),'dim') && ~isempty(eprops(i).dim));  err('Property dim is not defined in eprops(%u)',i); end
                                    validateattributes(eprops(i).dim,{'double'},{'vector','numel',2},'',sprintf('dim property in eprops(%u)',i));
                                    ensure(all(eprops(i).dim>=1e-4),sprintf('eprops(%i).dim values should be at least 1e-4 m.',i))
                                    if ~(isfield(eprops(i),'orien') && ~isempty(eprops(i).orien));  err('Property orien is not defined in eprops(%u)',i); end
                                    validateattributes(eprops(i).orien,{'double'},{'vector','numel',3},'',sprintf('orien property in eprops(%u)',i));
                                case 'circ'
                                    if ~(isfield(eprops(i),'dim') && ~isempty(eprops(i).dim));  err('Property dim is not defined in eprops(%u)',i); end
                                    validateattributes(eprops(i).dim,{'double'},{'vector','numel',1},'',sprintf('dim property in eprops(%u)',i));
                                    ensure(eprops(i).dim>=1e-4,sprintf('eprops(%i).dim value should be at least 1e-4 m.',i))
                            end
                        end
                        
                        if (isfield(eprops(i),'nbeams') && ~isempty(eprops(i).nbeams))
                            validateattributes(eprops(i).nbeams,{'double'},{'scalar','>=',1},'',sprintf('nbeams property in eprops(%u)',i));
                            ensure(mod(eprops(i).nbeams,1)==0,'Property eprops(%i).nbeams should be a positive integer.',i);
                        end
                        
                        if (isfield(eprops(i),'cw') && ~isempty(eprops(i).cw));     validateattributes(eprops(i).cw,{'logical'},{'scalar'},'',sprintf('cw property in eprops(%u)',i)); end
                        
                        %check for mandatory fields when dens field is present
                        if (isfield(eprops(i),'dens') && ~isempty(eprops(i).dens))
                            validateattributes(eprops(i).dens,{'double'},{'scalar'},'',   sprintf('dens property in eprops(%u)',i));
                            
                            %dens requires cshape
                            if ~(isfield(eprops(i),'cshape') && ~isempty(eprops(i).cshape))
                                err('Property cshape is not defined in eprops(%u)',i);
                            end
                        end
                        
                        %check for mandatory fields when flex field is present
                        if (isfield(eprops(i),'flex') && ~isempty(eprops(i).flex))
                            validateattributes(eprops(i).flex,{'double'},{'vector','positive'},'',sprintf('flex property in eprops(%u)',i));
                            
                            %check if values for flex are valid
                            validateattributes(eprops(i).flex,{'double'},{'vector'},'',  sprintf('flex property in eprops(%u)',i));
                            if any(((eprops(i).flex==1)+(eprops(i).flex==2)+(eprops(i).flex==3)+(eprops(i).flex==4)+(eprops(i).flex==5)+(eprops(i).flex==6))==0)
                                err('Invalid deformation mode in eprops(%u).flex.',i)
                            end
                            
                            %check if field exist in structure
                            if ~(isfield(eprops(i),'cshape') && ~isempty(eprops(i).cshape)); err('Property cshape is not defined in eprops(%u)',i);     end
                            if ~(isfield(eprops(i),'emod') && ~isempty(eprops(i).emod)); err('Property emod is not defined in eprops(%u)',i);     end
                            if ~(isfield(eprops(i),'smod') && ~isempty(eprops(i).smod)); err('Property smod is not defined in eprops(%u)',i);     end
                            if ~(isfield(eprops(i),'dens') && ~isempty(eprops(i).dens)); err('Property dens is not defined in eprops(%u)',i);     end
                        else
                            eprops(i).flex = [];
                            if (isfield(eprops(i),'emod') && ~isempty(eprops(i).emod)); warn('Property eprops(%u).emod is redundant without the flex property.',i);     end
                            if (isfield(eprops(i),'smod') && ~isempty(eprops(i).smod)); warn('Property eprops(%u).smod is redundant without the flex property.',i);     end
                            if (isfield(eprops(i),'cw') && ~isempty(eprops(i).cw)); warn('Property eprops(%u).cw is redundant without the flex property.',i);           end
                        end
                    end
                end
                
                %delete sets without the elems property (warning has already been issued)
                if any(ignoresets)
                    eprops(ignoresets)=[];
                end
                
                %warn user if elements are defined without adding properties to them
                el_without_prop_index = ~ismember(1:size(elements,1),el_nr_doubles_check);
                if any(el_without_prop_index)
                    el_without_prop = find(el_without_prop_index);
                    if length(el_without_prop) == 1
                        el_without_prop_str = num2str(el_without_prop);
                        warn('Element %s has no user-defined properties. Defaults (rigid massless elements) are used.',el_without_prop_str)
                    else
                        el_without_prop_str = [num2str(el_without_prop(1)) sprintf(', %i',el_without_prop(2:end))];
                        warn('Elements %s have no user-defined properties. Defaults (rigid massless elements) are used.',el_without_prop_str)
                    end
                end
                %warn user if no element has flexibility
                %if ~isfield(eprops,'flex') || cellfun(@isempty,{eprops(:).flex})
                if ~(isfield(eprops,'flex') || sum(cellfun(@isempty,{eprops(:).flex})))
                    warn('No element has the flex property. Simulation does not seem useful.')
                end
            end
            
            %start orien checks. Note: this comes *after* the dependency of cshape=rect on orien is ensured
            %loop over all elements, check if it belongs to a set with a orien property,
            %see if that property is valid. If it does not belong to a set, check if default works
            for i=1:size(elements,1)
                ii_set = arrayfun(@(x) ismember(i,x.elems),eprops);
                if ~any(ii_set) %element not in any set
                    orien_try = [0 1 0];
                    %error message if this turns out invalid:
                    orien_err = 'No orien property specified for element %i (because element not in any set). Default value [0 1 0] does not work, because (almost) parallel to element axis.';
                elseif isfield(eprops(ii_set),'orien') && ~isempty(eprops(ii_set).orien) %element in a set with orien
                    orien_try = eprops(ii_set).orien;
                    %error message if this turns out invalid:
                    orien_err = 'Orien property for element %i does not work, because (almost) parallel to element axis.';
                else %element in a set, but no orien
                    orien_try = [0 1 0];
                    %error message if this turns out invalid:
                    orien_err = 'No orien property specified for element %i. Default value [0 1 0] does not work, because (almost) parallel to element axis.';
                end
                %try to see if the planned (user-supplied or default) orien works
                xp = nodes(elements(i,1),:);
                xq = nodes(elements(i,2),:);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %this directly from SPACAR 2015 source:
                ex = xq - xp;
                ex = ex/norm(ex);
                ey_input = orien_try;
                ey = ey_input(:)/norm(ey_input);
                ex_proj = dot(ey,ex);
                noemer = sqrt(1-ex_proj^2);
                if noemer < 1e-5
                    err(orien_err,i);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            
            %CHECK OPTIONAL ARGUMENTS
            if (exist('opt','var') && ~isempty(opt))
                allowed_opts = {'filename','gravity','silent','calcbuck','showinputonly','loadsteps','rls','mode','transfer','calccompl','customvis'};
                supplied_opts = fieldnames(opt);
                unknown_opts_i = ~ismember(supplied_opts,allowed_opts);
                if any(unknown_opts_i)
                    err('Unknown opt field %s.',supplied_opts{unknown_opts_i});
                end
                
                if isfield(opt,'filename')
                    if isempty(opt.filename)
                        warn('Filename cannot be empty. Filename spacar_file is used instead.');
                    else
                        validateattributes(opt.filename,{'char'},{'vector'},'',            'filename property in opt');
                        if length(opt.filename) > 19
                            warn('Filename too long: maximum of 20 characters. Filename spacar_file is used instead.');
                            opt.filename = [];
                        elseif any(isspace(opt.filename))
                            err('No white-space is allowed in filename');
                        end
                    end
                end
                if (isfield(opt,'silent') && ~isempty(opt.silent))
                    validateattributes(opt.silent,{'logical'},{'scalar'},'',            'silent property in opt');   end
                if (isfield(opt,'calcbuck') && ~isempty(opt.calcbuck))
                    validateattributes(opt.calcbuck,{'logical'},{'scalar'},'',          'calcbuck property in opt'); end
                if (isfield(opt,'calccompl') && ~isempty(opt.calccompl))
                    validateattributes(opt.calccompl,{'logical'},{'scalar'},'',          'calccompl property in opt'); end
                if (isfield(opt,'gravity') && ~isempty(opt.gravity))
                    validateattributes(opt.gravity,{'double'},{'vector','numel',3},'',  'gravity property in opt');  end
                if isfield(opt,'loadsteps')
                    validateattributes(opt.loadsteps,{'double'},{'scalar','>=',1},'',   'loadsteps property in opt');
                    ensure(~isinf(opt.loadsteps),'Number of loadsteps (opt.loadsteps) cannot be Inf.')
                    ensure((mod(opt.loadsteps,1)==0 ),'Number of loadsteps (opt.loadsteps) has to be a positive integer.')
                    if opt.loadsteps >= 100
                        warn('A large number of load steps is specified; this might increase computation time and memory usage.');
                    end
                end
                
                %CHECK RLS INPUT VARIABLE
                if (isfield(opt,'rls') && ~isempty(opt.rls))
                    ensure(all(ismember(fieldnames(opt.rls),{'def'})),'Unknown field in rls; only def field is allowed.');
                    
                    for i=1:size(opt.rls,2)
                        if ~isempty(opt.rls(i).def)
                            validateattributes(opt.rls(i).def,{'double'},{'vector'},'',   sprintf('def property in rls(%u)',i));
                            if any(((opt.rls(i).def==1)+(opt.rls(i).def==2)+(opt.rls(i).def==3)+(opt.rls(i).def==4)+(opt.rls(i).def==5)+(opt.rls(i).def==6))==0)
                                err('Invalid deformation mode in rls(%u).',i)
                            end
                        end
                    end
                end
                
                %CHECK TRANSFER FIELD
                if (isfield(opt,'transfer') && ~isempty(opt.transfer))
                    %check whether appropriate type
                    ensure(islogical(opt.transfer) || (iscell(opt.transfer) && length(opt.transfer)==2),'opt.transfer should be logical or 2-field cell.')
                    
                    %check value for relative damping
                    if iscell(opt.transfer)
                        ensure(islogical(opt.transfer{1}),'opt.transfer{1} should be logical.');
                        ensure(isnumeric(opt.transfer{2}),'opt.transfer{2} should be scalar.');
                        ensure(length(opt.transfer{2})==1,'opt.transfer{2} should be scalar.');
                        ensure(opt.transfer{2}>=0,'Relative damping (opt.transfer{2}) should be >=0.');
                        ensure(~isinf(opt.transfer{2}),'Relative damping (opt.transfer{2}) cannot be Inf.');
                    end
                    
                    %opt.transfer can be two different types, always convert to cell
                    if islogical(opt.transfer)
                        opt.transfer = {opt.transfer};
                    end
                    
                    %CHECKS IN CASE TRANSFER == TRUE
                    if any(opt.transfer{1} == true)
                        %check 'clash' with loadstep specification
                        if isfield(opt,'loadsteps') && ~isempty(opt.loadsteps)
                            warn('Specification of load steps is ignored (default of 1 is used) since state-space equations are calculated (opt.transfer=true)');
                        end
                        
                        %check clash with buckling calculation
                        if isfield(opt,'calcbuck') && any(opt.calcbuck==true)
                            err('Simultaneous calculation of buckling load multipliers (opt.calcbuck) and state-space equations (opt.transfer) is not supported.');
                        end
                        
                        %check if at least one input and at least one output are specified
                        if isfield(nprops,'transfer_in') && isfield(nprops,'transfer_out')
                            ensure(any(~cellfun(@isempty,{nprops.transfer_in})),'Calculation of state-space equations requires at least one input (nprops.transfer_in)');
                            ensure(any(~cellfun(@isempty,{nprops.transfer_out})),'Calculation of state-space equations requires at least one output (nprops.transfer_out)');
                        else
                            err('Calculation of state-space equations requires at least one input (nprops.transfer_in) and output (nprops.transfer_out).');
                        end
                        
                        %check clash with some specified input displ, force, moment
                        clashinputs = {
                            'displ_x','displ_y','displ_z',...
                            'displ_initial_x','displ_initial_y','displ_initial_z',...
                            'force','force_initial','moment','moment_initial'...
                            };
                        for k=1:length(clashinputs)
                            if isfield(nprops,clashinputs(k))
                                ensure(all(cellfun('isempty',{nprops.(clashinputs{k})})),'Calculation of state-space equations (opt.transfer) requires absence of nprops.%s input.',clashinputs{k});
                            end
                        end
                    end
                end
            end
            
        end %END NOT-SILENT MODE BLOCK
        
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
                varargout{5} = opt;
        end
        
    end

    function print_dat(fileID,formatSpec,string)
        if ~isempty(regexp(string,'\n','once'))
            fprintf(fileID,formatSpec,string);
        end
    end

    function err(message,varargin)
        %custom error function to hide the backtrace stuff in command window
        errorstruct.stack.file = '';
        errorstruct.stack.name = 'spacarlight';
        errorstruct.stack.line = 1;
        errorstruct.identifier = '';
        message = [message '\n'];
        if nargin > 1
            errorstruct.message = sprintf(message,varargin{:});
        else
            errorstruct.message =sprintf(message);
        end
        
        fid_write=fopen([opt.filename '.log'],'w');
        fprintf(fid_write, 'Date: %s\n',date);
        fprintf(fid_write, 'Spacarlight version: %s\n',sl_version);
        fprintf(fid_write, 'Matlab version: %s\n\n\n\n',version);
        fprintf(fid_write, '------ SPACAR LIGHT LOG -----\n');
        fprintf(fid_write, 'Message displayed in command window: %s\n',errorstruct.message);
        
        try %#ok<TRYNC>
            if exist([opt.filename '.log'],'file')
                fid = fopen([opt.filename '.log']);
                log = textscan(fid,'%s','delimiter','\n');
                fclose(fid);
            end
            if exist([opt.filename '.dat'],'file')
                fid = fopen([opt.filename '.dat']);
                dat = textscan(fid,'%s','delimiter','\n');
                fclose(fid);
            end
            
            if exist('msg','var')
                fprintf(fid_write, 'Error: %s\n\n',msg.message);
                fprintf(fid_write, '\nError location:\n');
                for i=1:size(msg.stack,1)
                    fprintf(fid_write, '\nFile: %s\n',msg.stack(i).file);
                    fprintf(fid_write, 'Name: %s\n',msg.stack(i).name);
                    fprintf(fid_write, 'Line: %u\n',msg.stack(i).line);
                end
            end
            fprintf(fid_write, '\n\n\n------ SPACAR DAT INPUT -----\n');
            if exist('dat','var')
                for i=1:size(dat{1},1)
                    fprintf(fid_write, '%s\n',dat{1}{i});
                end
            else
                fprintf(fid_write, 'No datfile found');
            end
            fprintf(fid, '\n\n\n------ SPACAR LOG OUTPUT -----\n');
            if exist('log','var')
                for i=1:size(log{1},1)
                    fprintf(fid_write, '%s\n',log{1}{i});
                end
            else
                fprintf(fid_write, 'No logfile found');
            end
        end
        fclose(fid_write);
        error(errorstruct)
    end


    function warn(msg,varargin)
        %custom warning function to include sprintf syntax
        warning('on','all')
        if nargin > 1
            warning(msg,varargin{:});
        else
            warning(sprintf(msg)); %#ok<SPWRN> \n not supported
        end
        warning('off','all')
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

    function ensure_idelret(cond,msg,varargin)
        %custom assert function to hide the backtrace stuff in command window
        try
            %execute assert, but hide output
            assert(any(cond));
        catch
            idx = find(cond==0);
            nums = 'Element';
            for i=1:length(idx)
                nums = sprintf('%s %u',nums,idx(i));
            end
            %if assert failed, throw error
            err([nums ' ' msg],varargin{:});
        end
        %if assert succeeded, do nothing
    end

%% AUXILIARY FUNCTIONS
    function results = calc_results(E_list, id_inputf, id_inputx, nodes, eprops, opt, label_transfer_in, label_transfer_out, results)
        filename = opt.filename;
        nddof   = getfrsbf([filename '.sbd'],'nddof'); %number of dynamic DOFs
        t_list  =  1:getfrsbf([filename,'.sbd'],'tdef'); %list of timesteps
        lnp     = getfrsbf([filename,'.sbd'],'lnp'); %lnp data
        ln      = getfrsbf([filename,'.sbd'],'ln'); %lnp data
        
        if nddof == 0
            warn('No dynamic degrees of freedom.')
            return
        end
        
        results.ndof = getfrsbf([filename '.sbd'] ,'ndof'); %
        
        %CHECK BUCKLING SETTINGS
        calcbuck = false;
        if (isfield(opt,'calcbuck') && opt.calcbuck == 1)
            calcbuck = true;
            if id_inputx
                warn('Input displacement prescribed; buckling load multipliers are also with respect to reaction forces due to this input.');
            end
            if ~id_inputf
                warn('No external forces are prescribed. Buckling values are not calculated.');
                calcbuck = false;
            end
        end
        
        %PROCESS RESULTS PER LOADSTEP
        x       = getfrsbf([filename '.sbd'] ,'x');
        fxtot   = getfrsbf([filename '.sbd'] ,'fxt');
        M0_data = getfrsbf([filename '.sbm'] ,'m0');
        G0_data = getfrsbf([filename '.sbm'] ,'g0');
        K0_data = getfrsbf([filename '.sbm'] ,'k0');
        N0_data = getfrsbf([filename '.sbm'] ,'n0');
        nk = sqrt(size(K0_data,2));%number of elementen in K,M,G,N matrix if t_list > 1
        
        if length(t_list)==1 %#ok<*BDSCI>
            x = x';
            fxtot = fxtot';
        end
        
        %K = K0 + getfrsbf([filename '.sbm'] ,'n0', i) + G0;
        %C = getfrsbf([filename '.sbm'] ,'c0', t_list(i)) + getfrsbf([filename '.sbm'] ,'d0', t_list(i));
        
        %NODE DEPENDENT RESULTS
        for i=t_list
            for j=1:size(nodes,1)
                
                %RESTRUCT DATA
                if length(t_list) > 1
                    K0 = reshape(K0_data(i,:),nk,[]);
                    N0 = reshape(N0_data(i,:),nk,[]);
                    G0 = reshape(G0_data(i,:),nk,[]);
                    M0 = reshape(M0_data(i,:),nk,[]);
                else
                    K0 =K0_data;
                    N0 = N0_data;
                    G0 = G0_data;
                    M0 = M0_data;
                end
                
                if ~ismember((j-1)*2+1,ln) %if node not connected to an element
                    results.step(i).node(j).x = nodes(j,1:3);
                    results.node(j).x(1:3,i) = results.step(i).node(j);
                    continue;
                end
                %store results per loadstep, using "step" field
                results.step(i).node(j).p           = x(i,lnp((j-1)*2+1,1:3));
                results.step(i).node(j).r_eulzyx    = quat2eulang(x(i,lnp((j-1)*2+2,1:4)));
                results.step(i).node(j).r_axang     = quat2axang(x(i,lnp((j-1)*2+2,1:4)));
                results.step(i).node(j).r_quat      = x(i,lnp((j-1)*2+2,1:4));
                results.step(i).node(j).Freac       = fxtot(i,lnp((j-1)*2+1,1:3));
                results.step(i).node(j).Mreac       = fxtot(i,lnp((j-1)*2+2,2:4))/2;
                
                %also store results for all loadsteps combined
                results.node(j).p(1:3,i)             = results.step(i).node(j).p;
                results.node(j).r_eulzyx(1:3,i)      = results.step(i).node(j).r_eulzyx;
                results.node(j).r_axang(1:4,i)       = results.step(i).node(j).r_axang;
                results.node(j).r_quat(1:4,i)        = results.step(i).node(j).r_quat;
                results.node(j).Freac(1:3,i)         = results.step(i).node(j).Freac;
                results.node(j).Mreac(1:3,i)         = results.step(i).node(j).Mreac;
            end
            
            %EIGENFREQUENCIES
            
            if nddof>10
                [V,D]   = eigs(K0+N0+G0,M0,10,'sm');
            else
                [V,D]   = eig(K0+N0+G0,M0);
            end
            
            D       = diag(D);
            [~,o]   = sort(abs(D(:)));
            d       = D(o);
            results.step(i).freq = sqrt(d)*1/(2*pi); %per loadstep
            results.freq(1:length(d),i) = results.step(i).freq; %for all loadsteps
            
            %BUCKLING
            if calcbuck
                [~,loadmult] = eig(-K0,G0);
                results.step(i).buck = sort(abs(diag(loadmult))); %per loadstep
                results.buck(1:nddof,i) = results.step(i).buck; %for all loadsteps
            end
            
            %MAXIMUM STRESS
            [propcrossect, Sig_nums]  = calc_propcrossect(E_list,eprops);
            opt_stress.exterior = true; %only calculate exterior stresses (not possible for circ cross-section)
            [~,~,~,stressextrema] = stressbeam([filename,'.sbd'],Sig_nums,i,opt_stress,propcrossect);
            results.step(i).stressmax = stressextrema.max*1e6; %per loadstep
            results.stressmax(i) = results.step(i).stressmax; %for all loadsteps
            
            if opt.mode==9
                sys_ss = getss(filename);
                if length(opt.transfer) == 2 %relative damping has been specified
                    reldamp = opt.transfer{2};
                    nstates = size(sys_ss.a,1);
                    V = V*diag(1./sqrt(diag(V.'*M0*V))); %modeshapes mass-orthonormal
                    %             V*diag(D)*V'*M0 %M\K
                    %             V*2*reldamp*diag(sqrt(D))*V'*M0; %M\D
                    sys_ss.a(nstates/2+1:nstates,nstates/2+1:nstates) = ...
                        -V*2*reldamp*diag(sqrt(D))*V'*M0;
                end
                results.statespace = sys_ss;
            end
            for j=1:length(label_transfer_in)
                results.statespace.InputName{j} = label_transfer_in{j};
            end
            for j=1:length(label_transfer_out)
                results.statespace.OutputName{j} = label_transfer_out{j};
            end
        end
        
        if opt.calccompl
            for j=1:size(nodes,1)
                if ismember(j,elements) %Check if node in element list
                    [CMglob, CMloc] = complt(filename,(j-1)*2+1,(j-1)*2+2);
                    for i=t_list
                        results.step(i).node(j).CMglob = CMglob(:,:,i);
                        results.step(i).node(j).CMloc = CMloc(:,:,i);
                    end
                    for i=t_list
                        results.node(j).CMglob(1:6,1:6,i)    = results.step(i).node(j).CMglob;
                        results.node(j).CMloc(1:6,1:6,i)     = results.step(i).node(j).CMloc;
                    end
                else
                    %Node does not exist, do nothing
                end
            end
        end
    end

    function stiffness = calc_stiffness(eprops)
        % Compute the stiffness properties for rectangular or circular cross-section
        type    = eprops.cshape;
        dim     = eprops.dim;
        E       = eprops.emod;
        G       = eprops.smod;
        v       = E/(2*G) - 1;
        switch lower(type)
            case 'rect'
                w   = dim(1);
                t   = dim(2);
                A   = t*w;
                It 	= calc_torsStiff(t,w);
                Iy  = (1/12)*t^3*w;
                Iz  = (1/12)*t*w^3;
                k   = 10*(1+v)/(12+11*v);
            case 'circ'
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

    function Ip = calc_torsStiff(t, w)
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
        % Compute the inertia properties for rectangular or circular cross-section
        type    = eprops.cshape;
        dim     = eprops.dim;
        rho     = eprops.dens;
        switch lower(type)
            case 'rect'
                w   = dim(1);
                t   = dim(2);
                A   = t*w;
                Iy  = 1/12 * t^3*w;
                Iz  = 1/12 * t*w^3;
            case 'circ'
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

    function out = quat2axang(q)
        
        %conversion from quaternions (Euler parameters) to axis-angle representation
        %based on http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToAngle/
        %first output is angle (radians), followed by axis
        
        if q(1)>1, err('Quaternions should be normalized.'); end
        
        angle = 2*acos(q(1)); %returns in radians
        s = sqrt(1-q(1)^2);
        if s < 0.001
            %if s is close to zero, then axis is not important (just pick something normalised)
            x = 1;
            y = 0;
            z = 0;
        else
            x = q(2)/s;
            y = q(3)/s;
            z = q(4)/s;
        end
        
        out = [angle; x; y; z];
        
    end

    function eul = quat2eulang(q)
        
        %conversion from quaternions to Euler angles (radians)
        %rotation sequence is ZYX (following quat2eul from Robotics System Toolbox)
        
        %q = normc(q(:)); %normalize !requires specific toolboxes for matlab version 2017
        %versions
        q = q./norm(q);
        
        %extra check
        test = -2*(q(2)*q(4)-q(1)*q(3));
        if test>1, test = 1; end
        
        eul(1,1) = atan2(2*(q(2)*q(3)+q(1)*q(4)),q(1)^2+q(2)^2-q(3)^2-q(4)^2);
        eul(2,1) = asin(test);
        eul(3,1) = atan2(2*(q(3)*q(4)+q(1)*q(2)),q(1)^2-q(2)^2-q(3)^2+q(4)^2);
        
        if ~isreal(eul), eul = real(eul); end
        
    end

    function q = eul2quat(eul)
        %conversion from Euler angles (radians) to quaternions
        %rotation sequence is ZYX (following eul2quat from Robotics System Toolbox)
        
        % Pre-allocate output
        q = zeros(size(eul,1), 4, 'like', eul);
        
        % Compute sines and cosines of half angles
        c = cos(eul/2);
        s = sin(eul/2);
        
        % Construct quaternion
        q = [c(:,1).*c(:,2).*c(:,3)+s(:,1).*s(:,2).*s(:,3), ...
            c(:,1).*c(:,2).*s(:,3)-s(:,1).*s(:,2).*c(:,3), ...
            c(:,1).*s(:,2).*c(:,3)+s(:,1).*c(:,2).*s(:,3), ...
            s(:,1).*c(:,2).*c(:,3)-c(:,1).*s(:,2).*s(:,3)];
    end

    function q = axang2quat(axang)
        %conversion from axang to quaternions
        
        % Normalize the axis
        v = axang(1:3)./norm(axang(1:3));
        
        % Create the quaternion
        thetaHalf = axang(:,4)/2;
        sinThetaHalf = sin(thetaHalf);
        q = [cos(thetaHalf), v(1).*sinThetaHalf, v(2).*sinThetaHalf, v(3).*sinThetaHalf];
    end

    function [cw, aspect] = cw_values(L,eprops)
        w = eprops.dim(1);
        E = eprops.emod;
        G = eprops.smod;
        v = E/(2*G) - 1;
        
        aspect  = L/w;       %- aspect ratio of flexure
        c       = sqrt(24/(1+v));
        cw      = (aspect*c/(aspect*c-2*tanh(aspect*c/2)));
    end

    function [CMglob, CMloc] = complt(filename,ntr,nrot)
        % Calculate the compliance matrix in global directions and in body-fixed
        % local directions.
        % Calling syntax: [CMglob CMloc] = complt(filename,ntr,nrot)
        % Input:
        %   filename: string with the filename, without extension, of the file
        %             where the data can be found
        %   ntr:      node number of the translational node
        %   nrot:     node number of the rotational node
        % Output:
        %   CMglob: 6 x 6 x n compliance matrix in global directions
        %   CMloc:  6 x 6 x n compliance matrix in local directions
        
        if (nargin ~= 3)
            disp('complt needs 3 input arguments');
            CMglob=zeros(6,6,1);
            CMloc=zeros(6,6,1);
            return;
        end;
        
        % read total number of steps from the file
        tdef     =getfrsbf([filename '.sbd'],'tdef');
        CMglob=zeros(6,6,tdef);
        CMloc=zeros(6,6,tdef);
        lnp     =getfrsbf([filename '.sbd'],'lnp');
        nx      =getfrsbf([filename '.sbd'],'nx');
        nxp     =getfrsbf([filename '.sbd'],'nxp');
        nep     =getfrsbf([filename '.sbd'],'nep');
        DX_data      =getfrsbf([filename '.sbd'],'dx');
        K0_data      =getfrsbf([filename '.sbm'],'k0');
        G0_data      =getfrsbf([filename '.sbm'],'G0');
        X_data       =getfrsbf([filename '.sbd'],'x');
        
        
        
        % locate the place of the coordinates of the points where the compliance
        % has to be determined
        locv=[lnp(ntr,1:3), lnp(nrot,1:4)];
        % test whether the selected coordinates are feasible
        % for i=1:7
        %   if locv(i) <= 0
        %     disp('ERROR: invalid node number');
        %     return;
        %   end;
        %   if locv(i) <= nxp(1) || ...
        %      (locv(i)>(nxp(1)+nxp(2)) && locv(i) <= (nxp(1)+nxp(2)+nxp(3)))
        %     disp('WARNING: constrained node');
        %   end;
        % end;
        % search for the right degrees of freedom
        locdof=[nxp(3)+(1:nxp(4)), nxp(3)+nxp(4)+nep(3)+(1:nep(4))];
        
        if tdef > 1
            nddof = sqrt(length(K0_data(1,:)));
        end
        
        
        
        
        % start loop over the time steps
        for tstp=1:tdef
            % get data from files
            %DX      =getfrsbf([filename '.sbd'],'dx',tstp);
            %X       =getfrsbf([filename '.sbd'],'x',tstp);
            % reduce DX to the rows needed
            if tdef > 1
                DX = reshape(DX_data(tstp,:),nx,[]);
                %K0      =getfrsbf([filename '.sbm'],'k0',tstp);
                %G0      =getfrsbf([filename '.sbm'],'g0',tstp);
                K0 = reshape(K0_data(tstp,:),nddof,[]);
                G0 = reshape(G0_data(tstp,:),nddof,[]);
                X = X_data(tstp,:);
            else
                DX = DX_data;
                K0 = K0_data;
                G0 = G0_data;
                X = X_data;
            end
            DX = DX(locv,locdof);
            
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
            CMglob(:,:,tstp)=Tglob*CMlambda*(Tglob');
            CMloc(:,:,tstp)=Tloc*CMlambda*(Tloc');
        end
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
                    
                    switch eprops(id).cshape
                        case 'rect'
                            for j=1:length(Elements)
                                propcrossect(end+1).CrossSection = 'rect';
                                propcrossect(end).Dimensions = [eprops(id).dim(1),eprops(id).dim(2)];
                            end
                        case 'circ'
                            for j=1:length(Elements)
                                propcrossect(end+1).CrossSection = 'circ';
                                propcrossect(end).Dimensions = eprops(id).dim(1);
                            end
                    end
                end
            end
        end
    end

    function rls = restruct_rlse(rlse)
        for i=1:size(rlse,1)
            write = rlse(i,1:6).*[1 2 3 4 5 6];
            write(write==0) = [];
            rls(i).def = write;
        end
    end

end