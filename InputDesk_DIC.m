clc;clear;close all
addpath(genpath([pwd '\functions']));
pwdfor = 'C:\Users\ak13\Documents\240712_CT_Al_5052';
DataDirect = fullfile(pwdfor,'B0045.dat'); % file location
Maps.results = erase(DataDirect,'.dat');
% Domain size (square, crack tip at centre).
Maps.units.xy     = 'mm';            % meter (m) or milmeter (mm) or micrometer(um);
Maps.units.S      = 'Pa';      
Maps.units.St     = 'Pa'; 
Maps.pixel_size   = 1;              % if DIC values are in pixel, 1 if in physical units;
Maps.Operation    = 'DIC';          % Strain, xED = xEBSD, DIC = Displacement
Maps.stressstat   = 'plane_stress'; % 'plane_stress' OR 'plane_strain'
Maps.unique       = '240712_CT_Al_5052';

%% INPUT MATERIAL PROPERTIES AND DATA
% 'E' for Elastic material
% Poisson's ratio,          Young's Modulus [Pa],      		Material Name     
  Maps.nu    = 0.321;         Maps.E  = 70E9;                Maps.Mat = 'Al_5052';
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
% [theta,RawData.Ux,RawData.Uy,RawData.Uz,rotCentre] = ...
%     RotRemoving('true',Maps.X,Maps.Y,RawData.Ux,RawData.Uy);
[data] = Cropping10(Maps.X,Maps.Y,Maps.Ux, Maps.Uy);
Maps.X = data.xMap_px;
Maps.Y = data.yMap_px;
Maps.Ux = data.uxMap_px;
Maps.Uy = data.uyMap_px;
Maps.Uy(Maps.Uy==0)=NaN;    Maps.Ux(Maps.Ux==0)=NaN;
Maps.Uy = inpaint_nans(Maps.Uy);
Maps.Ux = inpaint_nans(Maps.Ux);
[theta,Maps.Ux,Maps.Uy,rotCentre] = ...
    RotRemoving('true',Maps.X,Maps.Y,Maps.Ux,Maps.Uy);
[J,KI,KII,KIII] = DIC2CAE(Maps);