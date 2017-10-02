clear

addpath ~/mlib/sp17ex

L = 80e-3;
w = 30e-3;
t = 0.5e-3;
d = 120e-3;

%% NODE POSITIONS
nodes = [0 0 0;
         0 0 L;
         d 0 L;
         d 0 0];

%% ELEMENT CONNECTIVITY
elements = [1 2;
            2 3;
            3 4];

%% NODE PROPERTIES
nprops(1).fix = true;
nprops(4).fix = true;
%nprops(2).fix_orien = true;
% nprops(2).force_initial = [0 10 0]; %initial load 1N 
% nprops(2).force = [-1e10 0 0]; %initial load 1N 

%% ELEMENT PROPERTIES
eprops(1).elems = [1 3];
eprops(1).emod = 210e9;
eprops(1).smod = 70e9;
eprops(1).dens = 7800;
eprops(1).dim = [w t];
eprops(1).cshape = 'rect';
eprops(1).flex = [1 2 3 4 5 6];
% eprops(1).orien = [0 0 1];
eprops(1).nbeams = 1;

%% RELEASES
rls(3).def = 1:6;

%% OPTIONAL
% opt.calcbuck = true;

%% DO SIMULATION
out = spacarlight(nodes,elements,nprops,eprops,rls);