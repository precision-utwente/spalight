addpath ~/mlib/sp17ex/

% EXAMPLE SCRIPT TO RUN SPACAR_LIGHT
clear
clc

%this example simulates a simple cross flexure rotating from -0.25 radians rotations
%to +0.25 radians rotation

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
%Structure with index i corresponding to the node number and the opt arguments:
%
%nprops(i).fix                  value true  to constraint translational and rotational motion of node i
%nprops(i).fix_pos              value true  to constraint translational motion of node i
%nprops(i).fix_x                value true  to constraint translational motion in x-direction of node i
%nprops(i).fix_y                value true  to constraint translational motion in y-direction of node i
%nprops(i).fix_z                value true  to constraint translational motion in z-direction of node i
%nprops(i).fix_orien            value true  to constraint rotational motion of node i
%
%nprops(i).force                vector of length 3 containing a force added in 10 loadsteps at node i in Newton, [Fx Fy Fx]
%nprops(i).force_initial        vector of length 3 containing an initial force at node i in Newton, [Fx Fy Fx]
%nprops(i).moment               vector of length 3 containing the moment components applied to node i in 10 loadsteps, in Newton*meter, [Mx My Mz]
%nprops(i).moment_initial       vector of length 3 containing the initial moment components applied to node i in Newton*meter, [Mx My Mz]
%
%nprops(i).displ_x            input displacement at node i in x-direction added in 10 steps in meters
%nprops(i).displ_y            input displacement at node i in y-direction added in 10 steps in meters
%nprops(i).displ_z            input displacement at node i in z-direction added in 10 steps in meters
%nprops(i).displ_initial_x    initial input displacement at node i in x-direction in meters
%nprops(i).displ_initial_y    initial input displacement at node i in y-direction in meters
%nprops(i).displ_initial_z    initial input displacement at node i in z-direction in meters
%nprops(i).rot_x              input rotation at node i around the x-axis added in 10 steps in radians*
%nprops(i).rot_y              input rotation at node i around the y-axis added in 10 steps in radians*
%nprops(i).rot_z              input rotation at node i around the z-axis added in 10 steps in radians*
%nprops(i).rot_initial_x      initial input rotationan at node i around the x-axis in radians*
%nprops(i).rot_initial_y      initial input rotationan at node i around the x-axis in radians*
%nprops(i).rot_initial_z      initial input rotationan at node i around the x-axis in radians*
%
%   *At each node, only rotation can be prescribed around a single axis!
%
%nprops(i).mass                 point-mass at node i in kg
%nprops(i).mominertia           vector with length 6 containing the moment of inertia of node i in Nm^2, [Ixx Ixy Ixz Iyy Iyz Izz]
%
%
%
%example:
%node 1
nprops(1).fix               = true;             %Fix node 1

%node 3
% nprops(3).force           = [0 1 0]; 
% nprops(3).force_initial   = [0 0.2 0];
nprops(3).moment_initial    = [0 0 0.01];       %Initial moment [Nm] around z-axis on node 3 
nprops(3).moment            = [0 0 0.05];       %Moment [Nm] around z-axis on node 3
nprops(3).mass              = 10;               %Mass of node 3: 10kg
nprops(3).mominertia        = [1 0 0 1 0 1];    %Inertia of node 3: Ixx = 1kgm^2, Iyy = 1kgm^2, Izz = 1kgm^2

%node 4
nprops(4).fix               = true;             %Fix node 4



%% ELEMENT PROPERTIES
%Structure with index i corresponding to the element number and the opt arguments:
%
%eprops(i).elems*              vector (list) containing the element numbers to which the properties of eprops(i) apply
%eprops(i).emod*                  E-modulus in Pa
%eprops(i).smod*                  G-modulus in Pa
%eprops(i).dens*                density in kg/m^3
%eprops(i).type*               a string providing the type of element. Supported options are 'leafspring','wire' or 'rigid'
%eprops(i).dim*                dimensions of the cross-section matching the element type in meters. For leafspring or rigid body a vector with length 2, [width thickness], for wire a single value [Diameter]
%eprops(i).orien                vector with length 3 containing the orientation of the  "width direction" of the  leafspring. Can be omited for type 'wire' or 'rigid', [Ox, Oy, Oz]
%eprops(i).nbeams              number of beam-elements used in the simulation. If omited, the defaults value 1 is used. 
%eprops(i).flex                 vector containing the deformation modes to be considered flexible. 1: elongation, 2: torsion, 3 & 4: out-of-plane bending, 5 & 6, inplane bending, [2 3 4]
%eprops(i).color                vector with length 3 containing the rgb values for the color of theelement, [1 0 0]
%eprops(i).hide                 value true (1) to hide this element in the visualization (it is included in the simulation)
%
%   *required input if the element is flexible (some deformation modes are flexible)
%
%
%example:

%Property set 1
eprops = [];
eprops(1).elems    = [1 3];            %Add this set of properties to element 1 and 3
eprops(1).emod     = 210E9;            %E-modulus
eprops(1).smod     = 70E9;             %G-modulus
eprops(1).dens     = 7800;             %Density
eprops(1).type     = 'leafspring';     %Simulate leafspring
eprops(1).dim      = [50e-3 0.2e-3];     %Width: 50mm, tickness: 1mm
eprops(1).orien    = [0 0 1];          %Width-direction of the leafspring in z-direction
eprops(1).nbeams   = 4;                %4 beam elements to simulate leafspring
eprops(1).flex     = [ 3 4 ];          %Model out-of-plane bending as flexible
eprops(1).color    = [0 1 0];          %Color: green

%Property set 1
eprops(2).elems    = 2;                %Add this set of properties to element 2
eprops(2).dens     = 3000;             %Density
eprops(2).type     = 'rigid';          %Simulate leafspring
eprops(2).dim      = [50e-3 10e-3];    %Width: 50mm, tickness: 10mm
eprops(2).orien    = [0 0 1];          %Width-direction of the body in z-direction
eprops(2).nbeams   = 1;                %1 beam elements to simulate body
eprops(2).color    = [1 0 0];          %Color: red
eprops(2).hide     = false;            %Body is visible, set to true to hide body



%% RELEASES
%Structure with index i corresponding to the element number and the arguments:
%
%rls(i).def*                    vector containng the deformation modes to be considered released. 1: elongation, 2: torsion, 3 & 4: out-of-plane bending, 5 & 6, inplane bending, [2 3 4]
%
%   *if multiple beams are used (eprops(i).n_beams > 1), the releases will only be added to the last beam number
%
%
%example:

rls(1).def = [1 2 3 4 5 6];                %Set of releases in element 1 to prevent overconstraints. 
%rls = [];                                  %Create empty Rlse variable to check for overconstraints

%% OPTIONAL
%Structure with opt arguments:
%
%opt.filename                     string containing the filemane for the spacar .dat file. Length must be smaller then 20 characters
%opt.silent                       value true (1) to supress visualization after the simulation is completed. Can be usefull when calling Spacar_light.m in a for-loop during optimization
%opt.buckload*                    value true (1) to calculate the load multipliers for buckling
%opt.gravity                      vector with length 3 providing the gravitation vector in m/s^2, [0 -9.81 0]
%
%   *evaluation of the load multipliers is only allowed when external force is applied on the system. Furthermore, load multipliers are also with respect to the reaction forces caused by prescribed
%   displacements/rotations. It is recommended to only evaluate load multipliers with a single force/moment applied to the system and no prescribed rotations/displacements.
%
%example:

opt.filename    = 'spacarfile';     %Filename
opt.silent      = false;            %Run in normal mode
opt.buckload    = false;             %Disable calcuation of load multipliers (default)
opt.gravity     = [0 0 -9.81];      %Gravitation in z-direction

%% CALL SPACAR_LIGHT
results = spacarlight(nodes, elements, nprops, eprops, rls, opt);

%simulation results are stored in RESULTS structure
%
%results.step(i)                            results at loadstep i (i=1 if no additional loads or displacements are specified)
%results.step(i).Freq                       List with eigenfrequencies in Hz, sorted from lowest to highest
%results.step(i).Buck*                      List with load multipliers, sorted from lowest to highest
%results.step(i).stressmax                  Maximum stress
%results.step(i).node(j)                    results at loadstep i, node number j
%results.step(i).node(j).x                  Position
%results.step(i).node(j).rx_eulzyx          Euler rotations in order z, y, x
%results.step(i).node(j).rx_quat            Rotations in quaternions
%results.step(i).node(j).Freac              Reaction forces in the node (non-zero for fixed or input nodes)
%results.step(i).node(j).Mreac              Reaction moments in the node (non-zero for fixed or input nodes)
%results.step(i).node(j).CMglob             6x6 compliance matrix of the node in the global frame
%results.step(i).node(j).CMloc              6x6 compliance matrix of the node in the local frame
%
%   *load multipliers only calculated if opt.buck_load = true.