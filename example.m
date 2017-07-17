addpath ../sp15ex/
addpath ../spaplus
addpath ../spaprivate

% EXAMPLE SCRIPT TO RUN SPACAR_LIGHT
clear
clc

%this example simulates a simple cross flexure rotating from -0.25 radians rotations
%to +0.25 radians rotation

%% NODE POSITIONS
%           x y z
Nodes = [   0 0 0;      %node 1  
            0 0.1 0;    %node 2  
            0.1 0.1 0;  %node 3
            0.1 0 0];   %node 4

        
%% ELEMENT CONNECTIVITY
%               p   q
Elements = [    1   3;  %element 1
                2   3;  %element 2
                2   4]; %element 3

            
%% NODE PROPERTIES  
%Structure with index i corresponding to the node number and the optional arguments:
%
%NODE_PROPS(i).fix_all              value true  to constraint translational and rotational motion of node i
%NODE_PROPS(i).fix_xyz              value true  to constraint translational motion of node i
%NODE_PROPS(i).fix_x                value true  to constraint translational motion in x-direction of node i
%NODE_PROPS(i).fix_y                value true  to constraint translational motion in y-direction of node i
%NODE_PROPS(i).fix_z                value true  to constraint translational motion in z-direction of node i
%NODE_PROPS(i).fix_rxyz             value true  to constraint rotational motion of node i
%NODE_PROPS(i).fix_rx             	value true  to constraint rotational motion around the x-axis of node i
%NODE_PROPS(i).fix_ry               value true  to constraint rotational motion around the y-axis of node i
%NODE_PROPS(i).fix_rz               value true  to constraint rotational motion around the z-axis of node i
%
%NODE_PROPS(i).force                vector with length 3 containing a force added in 10 loadsteps at node i in Newtons, [Fx Fy Fx]
%NODE_PROPS(i).force_initial        vector with length 3 containing an initial force at node i in Newtons, [Fx Fy Fx]
%
%NODE_PROPS(i).disp_x               input displacement at node i in x-direction added in 10 steps in meters
%NODE_PROPS(i).disp_y               input displacement at node i in y-direction added in 10 steps in meters
%NODE_PROPS(i).disp_z               input displacement at node i in z-direction added in 10 steps in meters
%NODE_PROPS(i).disp_initial_x       initial input displacement at node i in x-direction in meters
%NODE_PROPS(i).disp_initial_y       initial input displacement at node i in y-direction in meters
%NODE_PROPS(i).disp_initial_z       initial input displacement at node i in z-direction in meters
%NODE_PROPS(i).disp_rx              input rotation at node i around the x-axis added in 10 steps in radians*
%NODE_PROPS(i).disp_ry              input rotation at node i around the y-axis added in 10 steps in radians*
%NODE_PROPS(i).disp_rz              input rotation at node i around the z-axis added in 10 steps in radians*
%NODE_PROPS(i).disp_initial_rx      initial input rotationan at node i around the x-axis in radians*
%NODE_PROPS(i).disp_initial_ry      initial input rotationan at node i around the x-axis in radians*
%NODE_PROPS(i).disp_initial_rz      initial input rotationan at node i around the x-axis in radians*
%
%   *At each node, only rotation can be prescribed around a single axis!
%
%NODE_PROPS(i).mass                 point-mass at node i in kg
%NODE_PROPS(i).inertia              vector with length 6 containing the inertia properties of node i in Nm^2, [Ixx Ixy Ixz Iyy Iyz Izz]
%
%
%
%example:
%node 1
Node_props(1).fix_all           = true;             %Fix node 1
%node 3
Node_props(3).force_initial     = [0 -0.3 0];        %Initial moment around z-axis on node 3 
Node_props(3).force             = [0 0.6 0];         %Initial moment around z-axis on node 3
Node_props(3).mass              = 1;               %Mass of node 3: 10kg
Node_props(3).inertia           = [1 0 0 1 0 1];    %Inertia of node 3: Ixx = 1kgm^2, Iyy = 1kgm^2, Izz = 1kgm^2
%node 4
Node_props(4).fix_all           = true;             %Fix node 4



%% ELEMENT PROPERTIES
%Structure with index i corresponding to the element number and the optional arguments:
%
%ELEM_PROPS(i).El_Nrs*              vector (list) containing the element numbers to which the properties of ELEM_PROPS(i) apply
%ELEM_PROPS(i).E*                  E-modulus in Pa
%ELEM_PROPS(i).G*                  G-modulus in Pa
%ELEM_PROPS(i).rho*                density in kg/m^3
%ELEM_PROPS(i).type*               a string providing the type of element. Supported options are 'leafspring','wire' or 'rigid'
%ELEM_PROPS(i).dim*                dimensions of the cross-section matching the element type in meters. For leafspring or rigid body a vector with length 2, [width thickness], for wire a single value [Diameter]
%ELEM_PROPS(i).orien                vector with length 3 containing the orientation of the  "width direction" of the  leafspring. Can be omited for type 'wire' or 'rigid', [Ox, Oy, Oz]
%ELEM_PROPS(i).n_beams              number of beam-elements used in the simulation. If omited, the defaults value 1 is used. 
%ELEM_PROPS(i).flex                 vector containing the deformation modes to be considered flexible. 1: elongation, 2: torsion, 3 & 4: out-of-plane bending, 5 & 6, inplane bending, [2 3 4]
%ELEM_PROPS(i).color                vector with length 3 containing the rgb values for the color of theelement, [1 0 0]
%ELEM_PROPS(i).hide                 value true (1) to hide this element in the visualization (it is included in the simulation)
%
%   *required input if the element is flexible (some deformation modes are flexible)
%
%
%example:

%Property set 1
Elem_props(1).El_Nrs    = [1 3];            %Add this set of properties to element 1 and 3
Elem_props(1).E         = 210E9;            %E-modulus
Elem_props(1).G         = 70E9;             %G-modulus
Elem_props(1).rho       = 7800;             %Density
Elem_props(1).type      = 'leafspring';     %Simulate leafspring
Elem_props(1).dim       = [50e-3 0.2e-3];     %Width: 50mm, tickness: 1mm
Elem_props(1).orien     = [0 0 1];          %Width-direction of the leafspring in z-direction
Elem_props(1).n_beams   = 4;                %4 beam elements to simulate leafspring
Elem_props(1).flex      = [ 3 4 ];          %Model out-of-plane bending as flexible
Elem_props(1).color     = [0 1 0];          %Color: green

%Property set 1
Elem_props(2).El_Nrs    = 2;                %Add this set of properties to element 2
Elem_props(2).rho       = 3000;             %Density
Elem_props(2).type      = 'rigid';          %Simulate leafspring
Elem_props(2).dim       = [50e-3 10e-3];    %Width: 50mm, tickness: 10mm
Elem_props(2).orien     = [0 0 1];          %Width-direction of the body in z-direction
Elem_props(2).n_beams   = 1;                %1 beam elements to simulate body
Elem_props(2).color     = [1 0 0];          %Color: red
Elem_props(2).hide      = false;            %Body is visible, set to true to hide body



%% RELEASES
%Structure with index i corresponding to the element number and the arguments:
%
%RLSE(i).def*                    vector containng the deformation modes to be considered released. 1: elongation, 2: torsion, 3 & 4: out-of-plane bending, 5 & 6, inplane bending, [2 3 4]
%
%   *if multiple beams are used (ELEM_PROPS(i).n_beams > 1), the releases will only be added to the last beam number
%
%
%example:

Rlse(1).def = [1 2 3 4 5 6];                %Set of releases in element 1 to prevent overconstraints. 
%Rlse = [];                                  %Create empty Rlse variable to check for overconstraints

%% OPTIONAL
%Structure with optional arguments:
%
%OPTIONAL.filename                          string containing the filemane for the spacar .dat file. Length must be smaller then 20 characters
%OPTIONAL.simple_mode                       value true (1) to supress visualization after the simulation is completed. Can be usefull when calling Spacar_light.m in a for-loop during optimization
%OPTIONAL.buck_load*                        value true (1) to calculate the load multipliers for buckling
%OPTIONAL.gravity                           vector with length 3 providing the gravitation vector in m/s^2, [0 -9.81 0]
%
%optional input/outputs for transferfunctions (note prediscribing displacements possibly affects input/output transfer functions as prediscribed nodes are not free to move)
%
%OPTIONAL.transfer_in(i).type               type of input for the ith input for the transferfunction. Value is a string with options:
%                                           'force_x','force_y','force_z','moment_x','moment_y','moment_z',disp_x','disp_y','disp_z','rot_x','rot_y','rot_z'
%                                           which indicates a force/dipslacement in x,y or z-direction or a moment/rotation around the x,y or z-axis
%OPTIONAL.transfer_in(i).node               node to apply the ith input for the transferfunction
%OPTIONAL.transfer_out(i).type              type of output for the ith output for the transferfunction. Optional values similar to input
%OPTIONAL.transfer_out(i).node              node to apply the ith output for the transferfunction
%
%
%   *evaluation of the load multipliers is only allowed when external force is applied on the system. Furthermore, load multipliers are also with respect to the reaction forces caused by prescribed
%   displacements/rotations. It is recommended to only evaluate load multipliers with a single force/moment applied to the system and no prescribed rotations/displacements.
%
%
%example:

Optional.filename       = 'spacarfile';     %Filename
Optional.simple_mode    = false;            %Run in normal mode
Optional.buck_load      = false;             %Disable calcuation of load multipliers (default)
Optional.gravity        = [0 0 -9.81];      %Gravitation in z-direction

Optional.transfer_in(1).type = 'force_x';
Optional.transfer_in(1).node = 2;

Optional.transfer_out(1).type = 'disp_x';
Optional.transfer_out(1).node = 2;


%% CALL SPACAR_LIGHT
Results = Spacar_light(Nodes, Elements, Node_props, Elem_props, Rlse, Optional);

%simulation results are stored in RESULTS structure
%
%Results.step(i)                            Results at loadstep i (i=1 if no additional loads or displacements are specified)
%Results.step(i).Freq                       List with eigenfrequencies in Hz, sorted from lowest to highest
%Results.step(i).Buck*                      List with load multipliers, sorted from lowest to highest
%Results.step(i).stressmax                  Maxmimum stress
%Results.step(i).node(j)                    Results at loadstep i, node number j
%Results.step(i).node(j).x                  Position
%Results.step(i).node(j).rx_eulzyx          Euler rotations in order z, y, x
%Results.step(i).node(j).rx_quat            Rotations in quaternions
%Results.step(i).node(j).Freac              Reaction forces in the node (non-zero for fixed or input nodes)
%Results.step(i).node(j).Mreac              Reaction moments in the node (non-zero for fixed or input nodes)
%Results.step(i).node(j).CMglob             6x6 compliance matrix of the node in the global frame
%Results.step(i).node(j).CMloc              6x6 compliance matrix of the node in the local frame
%
%   *load multipliers only calculated if Optional.buck_load = true.

