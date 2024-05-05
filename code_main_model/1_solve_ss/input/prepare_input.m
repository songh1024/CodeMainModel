clear all;

load('input.mat')
load('area.mat','d_intra1');
d=d+d_intra1';

%% write Pi
fileID = fopen('Pi.dat','w');
print1(fileID,Pi,nmkt);

%% write Z
fileID = fopen('Z.dat','w');
print1(fileID,Z,nmkt);

%% write d
fileID = fopen('d.dat','w');
print1(fileID,d,nmkt);

%% write distance
fileID = fopen('distance.dat','w');
print2(fileID,distance,nmkt,nmkt);