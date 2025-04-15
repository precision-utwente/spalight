clear
clc

nodes = [   0 0 0
            0 0 0
            -0.1 0.1 0
            -0.1 0.1 0
            0 0.2 0
            0 0.2 0
            0.1 0.2 0
            0.1 0.2 0
            -0.1 0 0
            ];

elements = [    1   2 %hinge 1
                2   3 %rigid beam 1
                3   4 %hinge 2
                4   5 %rigid beam 2
                5   6 %hinge 3
                6   7 %rigid beam 3
                7   8 %hinge 4
                9   3
                ];

nprops(1).fix               = true;
nprops(8).fix               = true;
% nprops(3).force             = [1 0 0];
% nprops(3).displ_x = 0.15;

%rigid beams
eprops(1).elems    = [2 4 6];
eprops(1).dens     = 3000;
eprops(1).cshape   = 'rect';
eprops(1).dim      = [50e-3 10e-3];
eprops(1).orien    = [0 0 1];
eprops(1).nbeams   = 1;
eprops(1).color    = 'darkblue';
eprops(1).warping  = true;

%hinge elements, z-axis
eprops(2).elems = [1 3 5 7];
eprops(2).orien = [0 0 1];
eprops(2).type = 'hinge';

eprops(3).elems = 8;
eprops(3).orien = [0 0 1];
eprops(3).cshape = 'rect';
eprops(3).dim = [10e-3 10e-3];
eprops(3).orien = [0 0 1];
eprops(3).emod = 200e9;
eprops(3).smod = 70e9;
eprops(3).dens = 7800;

opt.filename    = 'testhingerigid';
opt.rls(2).def = 1:6;
% opt.rls(4).def = 1;

out = spacarlight(nodes, elements, nprops, eprops, opt);

%get rotation of hinge element 6 at step i:
% i = 8;
% out.step(i).element(6).e(1,1) %entry (1,1) of the e-field has the relative rotation of the hinge