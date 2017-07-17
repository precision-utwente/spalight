addpath ../sp15ex/
addpath ../spaplus
addpath ../spaprivate

% rmpath ../sp15ex/
% rmpath ../spaplus
% rmpath ../spaprivate

clear
clc

L = 0.5;
d = 0.2;

nodes = [
    0 0 0;
    0 L 0;
    d 0 0;
    ];
    
elements = [
    1 2;
%     2 3;
    ];

nprops(1).fix_all = true;
nprops(1).mass = 1;
nprops(3).fix_rx = true;
nprops(3).fix_ry = true

elprops(1).El_Nrs = [1];
elprops(1).E = 210e9;
elprops(1).G = 70e9;
elprops(1).rho = 7800;
elprops(1).dim = [0.03 0.0005];
elprops(1).type = 'leafspring';
elprops(1).flex = 1:6;
elprops(1).orien = [0 0 1];
elprops(1).n_beams = 1;

% elprops(2).El_Nrs = [2];
% elprops(2).orien = [0 0 1];

% rls(1).def = 1:5;
rls(1).def = [4];

opt.filename = 'test';
% opt.simple_mode = 1;

Spacar_light(nodes,elements,nprops,elprops,rls);