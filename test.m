% EXAMPLE SCRIPT TO RUN SPACAR_LIGHT
clear
clc


%% NODE POSITIONS
%           x y z
nodes = [   0 0 0
            6e-3 0 0
            34e-3 0 0
            40e-3 0 0];

        
%% ELEMENT CONNECTIVITY
elements = [    1 2
                2 3
                3 4];

% elements = [1 4];

            
%% NODE PROPERTIES 
nprops(1).fix = true;
%nprops(4).fix = true;
nprops(4).force_initial = [-1 1e-6 0];
% nprops(4).moment = [0 0 -0.4];

%% ELEMENT PROPERTIES

eprops(1).elems    = [1 3];
eprops(1).emod     = 210E9;
eprops(1).smod     = 79E9;
eprops(1).dens     = 7800;
eprops(1).cshape   = 'circ';
eprops(1).dim      = [1e-3];
eprops(1).orien    = [0 0 1];
eprops(1).flex     = [1 2 3 4 5 6];
eprops(1).nbeams =  2;


eprops(2).elems    = [2];
eprops(2).emod     = 210E9;
eprops(2).smod     = 79E9;
eprops(2).dens     = 7800;
eprops(2).cshape   = 'circ';
eprops(2).dim      = [1e-3];
eprops(2).orien    = [0 0 1];
eprops(2).flex     = [1 2 3 4 5 6];
eprops(2).nbeams = 2;


%% RELEASES

% rls(1).def = [];                %Set of releases in element 1 to prevent overconstraints. 
rls = [];                                  %Create empty Rlse variable to check for overconstraints

%% OPTIONAL
opt.calcbuck = true;
opt.filename = 'wireflexure';

%% CALL SPACAR_LIGHT
out = spacarlight(nodes, elements, nprops, eprops, rls, opt);