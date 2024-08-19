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
% Poisson's ratio,          Young's Modulus [Pa],      		Material Name     
  Maps.nu    = 0.3;         Maps.E  = 210E9;                Maps.Mat = 'Ferrite';
% 'E' for Elastic or 'R' for Ramberg-Osgood or 'A' for Elastic-Anisotropic
% or 'P' for elasto-plastic
  Maps.type  = 'E'; 
% if 'Ramberg-Osgood' type of material input                Yield Stress [Pa] 
%   Dir.Exponent = 26.67;     Dir.Yield_offset = 1.24;        Dir.yield = 4E9;
% if 'Elastic-Anisotropic' you need to define the stifness tensor 
%     Dir.Stiffness = NaN(6,6);
% if 'Plastic' you need to define the stifness tensor 
% Dir.type  = 'P'; 
% Dir.Plastic_Strain = [0.0010,0.0023,0.0038,0.0059,0.0085,0.0120,0.0150,0.0205,0.0260,0.0336]'-0.0010;
% Dir.Plastic_Stress = [150,350,547,849,1156,1569,1934,2616,3142,3912]'*1e6;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END of USER INTERFACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
Data = importdata(DataDirect);
[~,RawData ] = reshapeData(Data.data);
Maps.X  = RawData.X1;       Maps.Y = RawData.Y1;   
Maps.Ux = RawData.Ux;       Maps.Uy = RawData.Uy;
% for stereo DIC
% Maps.Z = RawData.Z1;      Maps.Uz = RawData.Uz;
[J,KI,KII,KIII] = DIC2CAE(Maps);