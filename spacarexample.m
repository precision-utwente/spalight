% EXAMPLE SCRIPT FOR RUNNING SPACAR LIGHT
% This example simulates a simple cross-flexure rotating due to an applied moment
% For details, more examples and the full syntax list: go to spacar.nl

clear
clc
% addpath('spacar')		%Specify location of spacar folder, if not current

%% NODE POSITIONS
%           x y z
nodes = [   0 0 0;      %node 1  
            0 0.1 0;    %node 2  
            0.1 0.1 0;  %node 3
            0.1 0 0];   %node 4

        
%% ELEMENT CONNECTIVITY
%               p   q
elements = [    1   3;  %element 1
                2   3;  %element 2
                2   4]; %element 3

            
%% NODE PROPERTIES  
%node 1
nprops(1).fix               = true;         %Fix node 1

%node 3
nprops(3).force             = [0 1 0];      %Force [N] in y-direction on node 3
nprops(3).moment_initial    = [0 0 -0.025]; %Initial moment [Nm] around z-axis on node 3 
nprops(3).moment            = [0 0 0.05];   %Moment [Nm] around z-axis on node 3 (combined with moment_initial, the moment goes from -0.025 to 0.025 Nm)
nprops(3).mass              = 0.1;          %Mass [kg] of node 3

%node 4
nprops(4).fix               = true;         %Fix node 4

%% ELEMENT PROPERTIES
%Property set 1 (for the two flexible leaf springs)
eprops(1).elems    = [1 3];            %Add this set of properties to elements 1 and 3
eprops(1).emod     = 210e9;            %E-modulus [Pa]
eprops(1).smod     = 70e9;             %G-modulus [Pa]
eprops(1).dens     = 7800;             %Density [kg/m^3]
eprops(1).cshape   = 'rect';           %Rectangular cross-section
eprops(1).dim      = [50e-3 0.2e-3];   %Width: 50 mm, thickness: 0.2 mm
eprops(1).orien    = [0 0 1];          %Orientation of the cross-section as a vector pointing along "width-direction"
eprops(1).nbeams   = 3;                %Number of beams used to model elements in this set
eprops(1).flex     = 1:6;        	   %Model all deformation modes (1 to 6) as flexible
eprops(1).color    = 'grey';           %Color (visualization only)
eprops(1).opacity  = 0.7;              %Opacity (visualization only)
eprops(1).warping  = true;             %Enable the modeling of warping (including constrained warping effects)

%Property set 2 (for the rigid body between the two leaf springs)
eprops(2).elems    = 2;                %Add this set of properties to element 2
eprops(2).dens     = 3000;             %Density [kg/m^3]
eprops(2).cshape   = 'rect';           %Rectangular cross-section
eprops(2).dim      = [50e-3 10e-3];    %Width: 50 mm, thickness: 10 mm
eprops(2).orien    = [0 0 1];          %Orientation of the cross-section as a vector pointing along "width-direction"
eprops(2).nbeams   = 1;                %Number of beams used to model elements in this set (as elements in this set are rigid, 1 is enough)
eprops(2).color    = 'darkblue';       %Color (visualization only)
eprops(2).warping  = true;             %Enable the modeling of warping (including constrained warping effects)
% eprops(2).hide     = true;           %Hide element (visualization only)

%% OPTIONAL ARGUMENTS
opt.filename    = 'crosshinge';     %Filename
opt.gravity     = [0 0 -9.81];      %Gravitational acceleration [m/s^2]
% opt.transfer    = {true 0.01};    %Enable calculation of transfer function (uses modal relative damping of 0.01)
% opt.calcbuck    = true;           %Enable calculation of load multipliers
% opt.calccompl   = false;          %Disable calculation of compliance matrices (can reduce computation time for large simulations)
% opt.showinputonly = true;         %Only visualize the elements and nodes that were specified (not running any simulation)
% opt.silent      = true;           %Run in silent mode

%% CALL SPACAR_LIGHT
out = spacarlight(nodes, elements, nprops, eprops, opt);

% Some of the available output variables are listed below. For the full syntax list, go to spacar.nl

%out.step(i)                            Results at loadstep i
%out.step(i).freq                       List with eigenfrequencies [Hz], sorted from lowest to highest
%out.step(i).buck                       List with load multipliers, sorted from lowest to highest*
%out.step(i).stressmax                  Maximum stress [Pa]
%out.step(i).node(j)                    Results at loadstep i, node number j
%out.step(i).node(j).p                  Position [m]
%out.step(i).node(j).r_axang			Rotation in axis-angle representation
%out.step(i).node(j).Freac              Reaction forces on the node
%out.step(i).node(j).Mreac              Reaction moments on the node
%out.step(i).node(j).CMglob             6x6 compliance matrix of the node in the global frame**

%   *load multipliers only calculated if opt.calcbuck = true.