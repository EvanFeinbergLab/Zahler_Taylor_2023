function [ Eh, Ev, tmpEyelid, tmpPupilArea ] = extract_pupil_trace(varargin)
%extract_pupil_trace takes deepLabCut marker locations and likelihoods as 
%input and converts pupil position to degrees. Uses lens mag to convert 
%model values from mm to pixels. Outputs pupil position in degrees.
% 
% REQUIRED INPUTS
%
% dlcPupilData should be array or table with these columns:
% [1,2]: Top eyelid (x,y)
% [3,4]: bottom eyelid (x,y)
% [5,6]: Pupil left (x,y)
% [7,8]: Pupil right (x,y)
% [9,10]: Pupil top (x,y)
% [11,12]: Pupil bottom (x,y)
% [13,14]: reflection left (x,y)
% [15,16]: reflection right (x,y)
%
% dlcLikelihoodData should be array or table with these columns:
% [1]: Top eyelid
% [2]: bottom eyelid
% [3]: Pupil left
% [4]: Pupil right
% [5]: Pupil top
% [6]: Pupil bottom
% [7]: reflection left
% [8]: reflection right
%
% lens_mag: pixels/mm conversion factor
%
% OPTIONAL PARAMETERS
%
% Method: 'sakatani' (default) or 'stahl'
%
% PupilMarkers: 'leftright' (default), 'topbottom', 'all'. Determines which DLC pupil markers should be used when computing the pupil center.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse varargin inputs
p = inputParser;
addRequired(p,'dlcPupilData');
addRequired(p,'dlcLikelihoodData');
addRequired(p,'lens_mag');
addParameter(p,'Method', 'sakatani');
addParameter(p,'PupilMarkers', 'leftright');
addParameter(p,'LikelihoodThreshold', 0.9);
parse(p, varargin{:})

dlcPupilData = p.Results.dlcPupilData;
dlcLikelihoodData = p.Results.dlcLikelihoodData;
lens_mag = p.Results.lens_mag;
method = p.Results.Method;
PupilMarkers = p.Results.PupilMarkers;
LikelihoodThreshold = p.Results.LikelihoodThreshold;


% convert tables to array
if istable(dlcPupilData) && istable(dlcLikelihoodData)
    dlcPupilData = table2array(dlcPupilData);
    dlcLikelihoodData = table2array(dlcLikelihoodData);
end

% Calculate pupil center
if strcmp(PupilMarkers, 'leftright')
    Xp = dlcPupilData(:,5)/2 + dlcPupilData(:,7)/2; % Pupil center (x)
    Yp = dlcPupilData(:,6)/2 + dlcPupilData(:,8)/2; % Pupil center (y)
elseif strcmp(PupilMarkers, 'topbottom')
    Xp = (dlcPupilData(:,9) + dlcPupilData(:,11))/2; % Pupil center (x)
    Yp = (dlcPupilData(:,10) + dlcPupilData(:,12))/2; % Pupil center (y)
elseif strcmp(PupilMarkers, 'all')
    Xp = (dlcPupilData(:,5) + dlcPupilData(:,7) + dlcPupilData(:,9) + dlcPupilData(:,11))/4; % Pupil center (x)
    Yp = (dlcPupilData(:,6) + dlcPupilData(:,8) + dlcPupilData(:,10) + dlcPupilData(:,12))/4; % Pupil center (y)
else
    error("Invalid PupilMarkers")
end

% Calculate reflection center
Xc = dlcPupilData(:,13)/2 + dlcPupilData(:,15)/2;  % Reflection center (x)
Yc = dlcPupilData(:,14)/2 + dlcPupilData(:,16)/2;  % Reflection center (y)

tmpEyelid = dlcPupilData(:,4) - dlcPupilData(:,2); 
tmpPupilRadius = (dlcPupilData(:,7) - dlcPupilData(:,5))/2; % radius along the horizontal axis
tmpPupilArea = (tmpPupilRadius.^2) .* pi;

% Compute angular eye position
if strcmp(method, 'sakatani')
    [Eh, Ev] = sakatani(Xp, Yp, Xc, Yc, tmpPupilRadius, lens_mag);
elseif strcmp(method, 'stahl')
    [Eh, Ev] = stahl(Xp, Yp, Xc, Yc, lens_mag);
else
    error("Invalid Method")
end

% Set low likelihood frames to NaN
if strcmp(PupilMarkers, 'leftright')
    nanFilter = dlcLikelihoodData(:,3)<LikelihoodThreshold | dlcLikelihoodData(:,4)<LikelihoodThreshold | dlcLikelihoodData(:,7)<LikelihoodThreshold | dlcLikelihoodData(:,8)<LikelihoodThreshold;
elseif strcmp(PupilMarkers, 'topbottom')
    nanFilter = dlcLikelihoodData(:,5)<LikelihoodThreshold | dlcLikelihoodData(:,6)<LikelihoodThreshold | dlcLikelihoodData(:,7)<LikelihoodThreshold | dlcLikelihoodData(:,8)<LikelihoodThreshold;
elseif strcmp(PupilMarkers, 'all')
    nanFilter = dlcLikelihoodData(:,3)<LikelihoodThreshold | dlcLikelihoodData(:,4)<LikelihoodThreshold | dlcLikelihoodData(:,5)<LikelihoodThreshold | dlcLikelihoodData(:,6)<LikelihoodThreshold | dlcLikelihoodData(:,7)<LikelihoodThreshold | dlcLikelihoodData(:,8)<LikelihoodThreshold;
end

Eh(nanFilter) = NaN;
Ev(nanFilter) = NaN;
tmpPupilArea(nanFilter) = NaN;
tmpEyelid(nanFilter) = NaN;



    function [Eh, Ev] = sakatani(Xp, Yp, Xc, Yc, PupilRadius, lens_mag)
        % Method by Sakatani and Isa (2004) 
        
        % Constants (from Sakatani and Isa)
        % lens_mag converts mm to pixels (lens mag = pixels/mm)
        Rlens = 1.25*lens_mag; % eye lens radius
        epsilon = 0.2*lens_mag; % distance between center of eyeball and center of cornea curvature
        delta = 0.1*lens_mag;  % distance between center of lens curvature and center of the eyeball
        
        % Calculate R_effective (the effective rotation radius)
        Reff = real(sqrt(Rlens^2 - PupilRadius.^2) - delta);
        
        % Calculate rotation center in video coordinates
        Xo = (Xc-Xp).*(Reff./(Reff-epsilon))+Xp;
        Yo = (Yc-Yp).*(Reff./(Reff-epsilon))+Yp;

        % Calculate angular eye position
        Eh = -real(asind((Xp-Xo)./(Reff)));
        Ev = -real(asind((Yp-Yo)./(Reff)));
    end

    function [Eh, Ev] = stahl(Xp, Yp, Xc, Yc, lens_mag)
        % Method by Stahl (2000)
        
        Reff = 1.08*lens_mag; % measured emperically
        
        Eh = -real(asind((Xp - Xc)./(Reff)));
        Ev = -real(asind((Yp - Yc)./(Reff)));
    end
    
end
