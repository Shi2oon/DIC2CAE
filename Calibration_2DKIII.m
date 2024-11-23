function [Maps,M4,alldata] = Calibration_2DKIII(KI,KII,KIII)
% for more details see https://ora.ox.ac.uk/objects/uuid:f2ba08f3-4a27-4619-92ed-bcd3834dadf0/files/d765371972, page 261
%% Input
              close all                   
% Domain size (square, crack tip at centre).
Maps.Mat          = 'Calibration';
Maps.type         = 'E';% 'A' if u want to use anistropic matrix or 'E' for linear elastic
Maps.input_unit   = 'um';        % meter (m) or milmeter (mm) or micrometer(um); 
Maps.pixel_size   = 1;           % if DIC values are in pixel, 1 if in physical units;
Maps.Operation    = 'DIC';       % Strain, xED = xEBSD, DIC = Displacement
Maps.stressstat   = 'plane_stress'; % 'plane_stress' OR 'plane_strain'
Maps.unique       = 'Calibration';
Maps.results =fullfile(pwd,Maps.unique);
sz = 100;

switch Maps.input_unit
    case 'm'
        saf = 1;
    case 'mm'
        saf = 1e3;
    case 'um'
        saf = 1e6;
    case 'nm'
        saf = 1e9;
end

Maps.E = 210e9;             % Young's Modulus
Maps.nu = 0.3;             % Poisson ratio
G = Maps.E/(2*(1 + Maps.nu));  % Shear modulus
    switch Maps.stressstat
        case 'plane_strain'
            kappa = 3 - (4 .* Maps.nu); % [/]
        case 'plane_stress'
            kappa = (3 - Maps.nu)./(1 + Maps.nu); % [/]
    end                                                           % Bulk modulus                                                  
% SIF loading
KI = KI*1e6;        % Mode I SIF
KII = KII*1e6;      % Mode II SIF
KIII = KIII*1e6;    % Mode III SIF

%% Anayltical displacement data.
Maps.stepsize = 1/sz*2;
lin = Maps.stepsize*(ceil(-1/Maps.stepsize)+1/2):Maps.stepsize:Maps.stepsize*(floor(1/Maps.stepsize)-1/2);
[Maps.X,Maps.Y,Maps.Z] = meshgrid(lin,lin,0);
[th,r] = cart2pol(Maps.X,Maps.Y);
DataSize = [numel(lin),numel(lin),1];
Maps.X = Maps.X*saf;
Maps.Y = Maps.Y*saf;
if KIII ~= 0
    Maps.Z = Maps.Z*saf;
end
% displacement data
Maps.Ux = ( 0.5*KI/G*sqrt(r/(2*pi)).*(+cos(th/2).*(kappa-cos(th)))+...
              0.5*KII/G*sqrt(r/(2*pi)).*(+sin(th/2).*(kappa+2+cos(th))))*saf;
Maps.Uy = ( 0.5*KI/G*sqrt(r/(2*pi)).*(+sin(th/2).*(kappa-cos(th)))+...
              0.5*KII/G*sqrt(r/(2*pi)).*(-cos(th/2).*(kappa-2+cos(th))))*saf;
if KIII ~= 0
    Maps.Uz = ( 2*KIII/G*sqrt(r/(2*pi)).*sin(th/2))*saf;
end
