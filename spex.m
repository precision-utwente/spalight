clear
clc

nodes = [   0 0 0
            0 0.1 0
            0.1 0.1 0
            0.1 0 0];

elements = [    1   3
                2   3
                2   4];

nprops(1).fix               = true;

nprops(3).force             = [0 1 0];
nprops(3).moment_initial    = [0 0 -0.025];
nprops(3).moment            = [0 0 0.05];
nprops(3).mass              = 0.1;

nprops(4).fix               = true;

eprops(1).elems    = [1 3];
eprops(1).emod     = 210e9;
eprops(1).smod     = 70e9;
eprops(1).dens     = 7800;
eprops(1).cshape   = 'rect';
eprops(1).dim      = [50e-3 0.2e-3];
eprops(1).orien    = [0 0 1];
eprops(1).nbeams   = 3;
eprops(1).flex     = 1:6;
eprops(1).color    = 'grey';
eprops(1).opacity  = 0.7;
eprops(1).warping  = true;

eprops(2).elems    = 2;
eprops(2).dens     = 3000;
eprops(2).cshape   = 'rect';
eprops(2).dim      = [50e-3 10e-3];
eprops(2).orien    = [0 0 1];
eprops(2).nbeams   = 1;
eprops(2).color    = 'darkblue';
eprops(2).warping  = true;
% eprops(2).hide     = true;

opt.filename    = 'spex';     %Filename
opt.gravity     = [0 0 -9.81];      %Gravitational acceleration [m/s^2]
% opt.transfer    = {true 0.01};    %Enable calculation of transfer function (uses modal relative damping of 0.01)
% opt.calcbuck    = true;           %Enable calculation of load multipliers
% opt.calccompl   = false;          %Disable calculation of compliance matrices (can reduce computation time for large simulations)
% opt.showinputonly = true;         %Only visualize the elements and nodes that were specified (not running any simulation)
% opt.silent      = true;           %Run in silent mode

out = spacarlight(nodes, elements, nprops, eprops, opt);