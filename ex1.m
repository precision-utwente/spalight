% EXAMPLE SCRIPT FOR RUNNING SPACAR LIGHT
% This example simulates a simple cross flexure rotating due to an applied moment

clear
clc
% addpath('spacar')		%Specify location of spacar folder, if not current

%% NODE POSITIONS
%           x y z
nodes = [   0 0 0;      %node 1  
            0 0 0.1;    %node 2  
            0.1 0 0.1;  %node 3
            0.1 0 0;    %node 4
            0.05 0 0.05
            ];   %node 5

        
%% ELEMENT CONNECTIVITY
%               p   q
elements = [    1   2;  %element 1
                2   3;  %element 2
                3   4; %element 3
                2   5
                ];  %element 2

            
%% NODE PROPERTIES  
%node 1
nprops(1).fix               = true;         %Fix node 1
nprops(1).fix_warp = true;

%node 3
nprops(3).force             = [150 0 0];      %Force [N] in y-direction on node 3


%node 4
nprops(4).fix               = true;         %Fix node 4
nprops(4).fix_warp = true;

%% ELEMENT PROPERTIES
%Property set 1
eprops(1).elems    = [1 3];            %Add this set of properties to elements 1 and 3
eprops(1).emod     = 210e9;            %E-modulus [Pa]
eprops(1).smod     = 69e9;             %G-modulus [Pa]
eprops(1).dens     = 7800;             %Density [kg/m^3]
eprops(1).cshape   = 'rect';           %Rectangular cross-section
eprops(1).dim      = [40e-3 1e-3];   %Width: 50 mm, thickness: 0.2 mm
eprops(1).orien    = [0 1 0];          %Orientation of the cross-section as a vector pointing along "width-direction"
eprops(1).nbeams   = 1;                %Number of beams used to model this element 
eprops(1).flex     = 1:6;        	   %Model out-of-plane bending (modes 3 and 4) as flexible
eprops(1).color    = 'grey';           %Color
eprops(1).opacity  = 0.7;              %Opacity
eprops(1).warping  = true;

%Property set 2
eprops(2).elems    = [2 4];            %Add this set of properties to element 2
eprops(2).dens     = 2700;             %Density [kg/m^3]
eprops(2).cshape   = 'rect';           %Rectangular cross-section
eprops(2).dim      = [40e-3 5e-3];     %Width: 50 mm, thickness: 10 mm
eprops(2).orien    = [0 1 0];          %Orientation of the cross-section as a vector pointing along "width-direction"
eprops(2).nbeams   = 1;                %1 beam for simulating this element (as it is rigid an no more elements are required)
eprops(2).color    = 'darkblue';
% eprops(2).hide     = true;           %Hide element (in visualization only)
% eprops(2).flex = 1:6;
eprops(2).emod = 210e9;
eprops(2).smod = 70e9;
eprops(2).warping = true;


%% OPTIONAL ARGUMENTS
% opt.filename    = 'crosshinge';     %Filename
%opt.gravity     = [0 0 -9.81];      %Gravitational acceleration [m/s^2]
% opt.calcbuck    = true;             %Enable calculation of load multipliers
%opt.calccompl   = false;            %Disable calculation of compliance matrices (can reduce computation time for large simulations)
% opt.showinputonly = true;          %Only visualize the elements and nodes that were defined (not running any simulation)
%opt.silent      = true;            %Run in silent mode

%% CALL SPACAR_LIGHT
out = spacarlight(nodes, elements, nprops, eprops);

% ah = out.fighandle.Children(1);
% axtoolbar(ah,{'export','datacursor','rotate','pan','zoomin','zoomout','restoreview'});

%out.step(i)                            Results at loadstep i
%out.step(i).freq                       List with eigenfrequencies [Hz], sorted from lowest to highest
%out.step(i).buck*                      List with load multipliers, sorted from lowest to highest
%out.step(i).stressmax                  Maximum stress [Pa]
%out.step(i).node(j)                    Results at loadstep i, node number j
%out.step(i).node(j).p                  Position [m]
%out.step(i).node(j).r_axang			Rotation in axis-angle representation
%out.step(i).node(j).r_eulzyx           Euler rotations in order z, y, x
%out.step(i).node(j).r_quat             Rotations in quaternions
%out.step(i).node(j).Freac              Reaction forces on the node
%out.step(i).node(j).Mreac              Reaction moments on the node
%out.step(i).node(j).CMglob             6x6 compliance matrix of the node in the global frame**
%out.step(i).node(j).CMloc              6x6 compliance matrix of the node in the local frame**
%
%   *load multipliers only calculated if opt.calcbuck = true.
%   **Compliance matrixes only calculated if opt.calccompl = true (if not specified, the default option "true" will be used)