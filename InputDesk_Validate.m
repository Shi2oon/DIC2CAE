% Validation from sythetic data
clc;clear;close all
addpath(genpath([pwd '\functions']));
Maps = Calibration_2DKIII(5,2,1);
[J,KI,KII,KIII] = DIC2CAE(Maps);