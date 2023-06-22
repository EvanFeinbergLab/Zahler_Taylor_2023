function [ Eh, Ev, tmpEyelid, tmpPupilArea ] = extract_pupil_trace(dlcPupilData, dlcLikelihoodData, lens_mag, method)
    
%extract_pupil_trace takes deepLabCut marker locations and likelihoods as 
%input and converts pupil position to degrees. Uses lens mag to convert 
%model values from mm to pixels. Outputs pupil position in degrees.

Xp = dlcPupilData(:,5)/2 + dlcPupilData(:,7)/2; % Pupil center (x)
Xc = dlcPupilData(:,13)/2 + dlcPupilData(:,15)/2;  % Reflection center (x)
Yp = dlcPupilData(:,6)/2 + dlcPupilData(:,8)/2; % Pupil center (y)
Yc = dlcPupilData(:,14)/2 + dlcPupilData(:,16)/2;  % Reflection center (y)
tmpEyelid = dlcPupilData(:,4) - dlcPupilData(:,2); 
tmpPupilRadius = (dlcPupilData(:,7) - dlcPupilData(:,5))/2; % radius along the horizontal axis
tmpPupilArea = (tmpPupilRadius.^2) .* pi;

% Compute angular eye position
if nargin == 3
    [Eh, Ev] = sakatani(Xp, Yp, Xc, Yc, tmpPupilRadius, lens_mag);
elseif nargin == 4
    if strcmp(method, 'sakatani')
        [Eh, Ev] = sakatani(Xp, Yp, Xc, Yc, tmpPupilRadius, lens_mag);
    elseif strcmp(method, 'stahl')
        [Eh, Ev] = stahl(Xp, Yp, Xc, Yc, lens_mag);
    else
        error("method does not exist")
    end
end

% Set low likelihood frames to NaN
dlcThreshold = 0.9; % DLC likelihood threshold
nanFilter = dlcLikelihoodData(:,3)<dlcThreshold | dlcLikelihoodData(:,4)<dlcThreshold | dlcLikelihoodData(:,7)<dlcThreshold | dlcLikelihoodData(:,8)<dlcThreshold;

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
