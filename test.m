clear
clc

L = 0.1;    %[m]
W = 0.05;   %[m]

%% NODE POSITIONS
nodes = [0 0 0;     %node 1
         0 L 0;     %node 2
         W L 0;     %node 3
         W 0 0];    %node 4
     

%% ELEMENT CONNECTIVITY
elements= [1 2;     %leafspring between node 1 and 2 
           2 3;     %intermediate body
           3 4];    %leafspring between node 3 and 4 



%% NODE PROPERTIES
nprops(1).fix = true;       %fix node 1
nprops(4).fix = true;   %fix node 4

nprops(2).displ_initial_x =-0.01; %start with node 2 displaced 10mm to the left
nprops(2).displ_x =         0.02; %displace node 2 20mm to the right

nprops(2).force_initial = [0 0 5]; %initial load of 5N in z-direction 
nprops(3).force_initial = [0 0 5]; %initial load of 5N in z-direction 


%% ELEMENT PROPERTIES

%first element property set
eprops(1).elems = [1 3];        %assing property set 1 to element 1 and 3
eprops(1).emod = 210e9;         %E-modulus [Pa]
eprops(1).smod = 70e9;          %G-modulus [Pa]
eprops(1).dens = 7800;          %density   [kg/m^3]
eprops(1).dim = [0.05 0.0005];  %crossectional dimension [w t] in [mm]
eprops(1).cshape = 'rect';      %crossectional shape rectangular
eprops(1).flex = [2 3 4];       %torsional (2) and out-of-plane bending (3,4) deformations are flexible
eprops(1).orien = [0 0 1];      %width-direction of leafspring oriented in z-direction
eprops(1).nbeams = 2;           %Leafspring simulated with 2 spacar-beams for increased accuracy. 

%second element property set
eprops(2).elems = [2];          %assing property set 2 to element 2
eprops(2).dens = 2700;          %density   [kg/m^3]
eprops(2).dim = [0.05 0.05];    %crossectional dimension [w t] in [mm]
eprops(2).cshape = 'rect';      %crossectional shape rectangular
eprops(2).orien = [0 0 1];      %width-direction of leafspring oriented in z-direction




%% RELEASES
rls(1).def = [1 2 3 4 5 6];
rls(3).def = [3];


%% DO SIMULATION
out = spacarlight(nodes,elements,nprops,eprops,rls);