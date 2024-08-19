% Validation from sythetic data
clc;clear;close all
addpath(genpath([pwd '\functions']));
% put the KI, KII, KII values for stereo-DIC data as in
% Calibration_2DKIII(KI,KII,KIII);for 2D, put KIII as zerp
Maps = Calibration_2DKIII(5,1,3);
[J,KI,KII,KIII] = DIC2CAE(Maps);