clear
clc

nodes = [   0 0 0
            -0.1 0.1 0
            0 0.2 0
            0 0.2 0 %note: this node coincides with the previous because they will be connected with a hinge
            0 0.2 0 %note: this node coincides with the previous because they will be connected with a hinge
            0.1 0.2 0
            0.1 0.2 0 %note: this node coincides with the previous because they will be connected with a hinge
            ];

elements = [    1   2 %flexible beam
                2   3 %rigid beam
                3   4 %hinge z-axis
                4   5 %hinge y-axis
                5   6 %rigid beam
                6   7 %hinge z-axis
                ];

%note: since elements 3 and 4 are in series, they constitute a universal
%joint in the z-y-plane

nprops(1).fix               = true;
nprops(7).fix               = true; %fixing one end of the last z-axis hinge

nprops(2).force             = [1 1 0];

%flexible beam
eprops(1).elems    = 1;
eprops(1).emod     = 210e9;
eprops(1).smod     = 70e9;
eprops(1).dens     = 7800;
eprops(1).cshape   = 'rect';
eprops(1).dim      = [50e-3 0.2e-3];
eprops(1).orien    = [0 0 1];
eprops(1).nbeams   = 1;
eprops(1).flex     = 1:6;
eprops(1).color    = 'grey';
eprops(1).opacity  = 0.7;
eprops(1).warping  = true;

%rigid beam
eprops(2).elems    = [2 5];
eprops(2).dens     = 3000;
eprops(2).cshape   = 'rect';
eprops(2).dim      = [50e-3 10e-3];
eprops(2).orien    = [0 0 1];
eprops(2).nbeams   = 1;
eprops(2).color    = 'darkblue';
eprops(2).warping  = true;

%hinge elements, z-axis
eprops(3).elems = [3 6];
eprops(3).orien = [0 0 1];
eprops(3).type = 'hinge';

%hinge element, y-axis
eprops(4).elems = 4;
eprops(4).orien = [0 1 0];
eprops(4).type = 'hinge';

opt.filename    = 'testhinge';

out = spacarlight(nodes, elements, nprops, eprops, opt);

%get rotation of hinge element 6 at step i:
i = 8;
out.step(i).element(6).e(1,1) %entry (1,1) of the e-field has the relative rotation of the hinge