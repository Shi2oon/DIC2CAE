clc;clear;close all
addpath(genpath([pwd '\functions']));
DataDirect = fullfile(pwd,'2D_DIC.dat'); % file location
Maps.results = fullfile(pwd,'2D_DIC');
% Domain size (square, crack tip at centre).
Maps.units.xy     = 'mm';            % meter (m) or milmeter (mm) or micrometer(um);
Maps.units.S      = 'Pa';      
Maps.units.St     = 'Pa'; 
Maps.pixel_size   = 1;              % if DIC values are in pixel, 1 if in physical units;
Maps.Operation    = 'DIC';          % Strain, xED = xEBSD, DIC = Displacement
Maps.stressstat   = 'plane_stress'; % 'plane_stress' OR 'plane_strain'
Maps.unique       = 'Calibration';

%% INPUT MATERIAL PROPERTIES AND DATA
% 'E' for Elastic material
% Poisson's ratio,          Young's Modulus [Pa],      		Material Name     
  Maps.nu    = 0.3;         Maps.E  = 210E9;                Maps.Mat = 'Ferrite';
  Maps.type  = 'E'; 
% if 'Ramberg-Osgood' type of material input                Yield Stress [Pa] 
% Maps.Exponent = 26.67;      Maps.Yield_offset = 1.24;        Maps.yield = 4E9;
% Maps.type  = 'R';           Maps.Mat = 'Austenite';
% Maps.E  = 193E9;            Maps.nu    = 0.3;
% if 'Elastic-Anisotropic' you need to define the stifness tensor 
% Maps.type  = 'A'; 
% Maps.Stiffness = [  283  121  121   0   0   0
%                     121  283  121   0   0   0
%                     121  121  283   0   0   0
%                     0    0    0     81  0   0
%                     0    0    0     0   81  0
%                     0    0    0     0   0   81]*1e9;
% Maps.Mat = 'Ferrite';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END of USER INTERFACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
Data = importdata(DataDirect);
[~,RawData ] = reshapeData(Data.data);
Maps.X  = RawData.X1;       Maps.Y = RawData.Y1;   
Maps.Ux = RawData.Ux;       Maps.Uy = RawData.Uy;
% for stereo DIC
% Maps.Z = RawData.Z1;      Maps.Uz = RawData.Uz;
[J,KI,KII,KIII] = DIC2CAE(Maps);