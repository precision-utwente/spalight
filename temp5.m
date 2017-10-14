clear
clc

addpath ~/mlib/sp17ex
addpath ~/mlib/spalight075/

L = 0.1;    %[m]
W = 0.1;   %[m]
d = 0.01;

%% NODE POSITIONS
nodes = [0 0 0;     %node 1
         0 0 L;     %node 2
         0.3*W 0 L;     %node 3
         W 0 L;
         W 0 0;
         0 0 -d;
         W 0 -d;
         0.3*W 0 L/2];    %node 4
 
 
%% ELEMENT CONNECTIVITY
elements= [1 2;     %leafspring between node 1 and 2 
           2 3;     %intermediate body
           3 4;     %intermediate body
           4 5;     %leafspring between node 3 and 4
           1 6;
           5 7;
           3 8];     
 
 
 
%% NODE PROPERTIES
nprops(1).fix = true;       %fix node 1
nprops(5).fix = true;       %fix node 4
 
% nprops(8).displ_x = 20e-3; %start with node 2 displaced 10mm to the left
nprops(8).force = [0 0 1];

% nprops(2).displ_x =         0.02; %displace node 2 20mm to the right
 
% nprops(2).force_initial = [0 0 5]; %initial load of 5N in z-direction 
% nprops(3).force_initial = [0 0 5]; %initial load of 5N in z-direction 
 
 
%% ELEMENT PROPERTIES
 
%first element property set
eprops(1).elems = [1 4];        %assing property set 1 to element 1 and 3
eprops(1).emod = 210e9;         %E-modulus [Pa]
eprops(1).smod = 70e9;          %G-modulus [Pa]
eprops(1).dens = 7800;          %density   [kg/m^3]
eprops(1).dim = [0.03 0.0005];  %crossectional dimension [w t] in [mm]
eprops(1).cshape = 'rect';      %crossectional shape rectangular
eprops(1).flex = [2 3 4];       %torsional (2) and out-of-plane bending (3,4) deformations are flexible
eprops(1).orien = [0 1 0];      %width-direction of leafspring oriented in z-direction
eprops(1).nbeams = 4;           %Leafspring simulated with 2 spacar-beams for increased accuracy. 
eprops(1).color = [0.8549    0.8588    0.8667];

%second element property set
eprops(2).elems = [2 3 7];          %assing property set 2 to element 2
eprops(2).dens = 2700;          %density   [kg/m^3]
eprops(2).dim = [0.035 0.01];   %crossectional dimension [w t] in [mm]
eprops(2).cshape = 'rect';      %crossectional shape rectangular
eprops(2).orien = [0 1 0];      %width-direction of leafspring oriented in z-direction
eprops(2).color = [0.1686    0.3922    0.6627];
 
%third element property set
eprops(3).elems = [5 6];        %assing property set 2 to element 2
eprops(3).dens = 2700;          %density   [kg/m^3]
eprops(3).dim = [0.035 0.01];   %crossectional dimension [w t] in [mm]
eprops(3).cshape = 'rect';      %crossectional shape rectangular
eprops(3).orien = [0 1 0];      %width-direction of leafspring oriented in z-direction
eprops(3).color = [0.2980    0.3020    0.3098];

%% RELEASES
opt.rls(1).def = [1 2 3 4 5 6];
%opt.rls(4).def = [4];
% rls = []; 

opt.calcbuck = true;

%% DO SIMULATION
out = spacarlight(nodes,elements,nprops,eprops,opt);

% set(gca(1),'view',[-15.9000   19.6000])