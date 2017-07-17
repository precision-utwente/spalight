addpath ../sp15ex/
addpath ../spaplus
addpath ../spaprivate

clear
clc

L = 0.5;
d = 0.2;

nodes = [
    0 0 0;
    0 L 0;
    0 d 0;
    ];
    
elements = [
    1 2;
    1 3;
    ];

nprops(2).fix_all = true;
nprops(1).mass = 1;
nprops(3).fix_all = true;

elprops(1).E = 210e9;
elprops(1).G = 70e9;
elprops(1).rho = 7800;
elprops(1).dim = [0.03 0.0005];
elprops(1).type = 'leafspring';
elprops(1).El_Nrs = [1];
elprops(1).flex = 1:6;
elprops(1).orien = [0 0 1];

elprops(2).El_Nrs = [2];
elprops(2).orien = [0 0 1];

rls(1).def = 1:6;

Spacar_light(nodes,elements,nprops,elprops,rls);