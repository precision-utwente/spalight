addpath ../spaplus
addpath ../arrow3

close all
% clear all

nodes = [
    0 0 0;
    1 1 0;
    1 1 1;
    1 0.5 1;
    1 0 1;
    ];

elements = [
    1 2;
    1 3;
    2 3;
    3 4;
%     4 5;
    ];

nopr = [];
nopr(2).force = [1 0 0];
% nopr(4).force = [0 0 1];

showGeom(nodes,elements,nopr);

set(gcf,'Position',[4   593   845   512])