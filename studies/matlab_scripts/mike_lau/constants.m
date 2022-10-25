classdef constants
    % Constants in cgs units unless otherwise specified
    
    properties (Constant)
        %------------------------------------------------------------------
        % Physical constants and parameters
        %------------------------------------------------------------------
        G           = 6.67408e-8;       % Gravitational constant
        c           = 2.99792458e10;    % Speed of light
        MSUN        = 1.98847e33;       % Solar mass
        RSUN        = 6.957e10;         % Solar radius
        LSUN        = 3.839e33;         % Solar luminosity
        AU          = 1.495978707e13;   % Astronomical unit
        MJ          = 1.89813e30;       % Jupiter mass
        RJ          = 7.1492e9;         % Jupiter radius
        RADCONST    = 7.5646e-15;       % Radiation constant
        STEFBOLTZ   = 5.6704e-5;        % Stefan-Boltzmann constant
        KB          = 1.38066e-16;      % Boltzmann constant
        PROTONMASS  = 1.67262158e-24;   % Proton mass
        AVOGADRO    = 6.0221408577e23;  % Avogadro's number
        
        %------------------------------------------------------------------
        % Unit conversion factors between code units (length in RSUN, mass
        % in MSUN) and cgs units
        %------------------------------------------------------------------
        TCODE             = sqrt(constants.RSUN^3 / (constants.G * constants.MSUN));    % Code units for time
        CM_S_PER_VCODE    = constants.RSUN / constants.TCODE;                           % cgs velocity unit per code unit
        G_CM3_PER_RHOCODE = constants.MSUN * constants.RSUN^-3;                         % cgs density unit per code unit
        YRS_PER_TCODE     = constants.TCODE / 60^2 / 24 / 365;                          % cgs density unit per code unit
        HR_PER_TCODE      = constants.TCODE / 60^2;
        SEC_PER_YR        = 3.15576e7;
        DYNE_PER_FCODE    = constants.MSUN * constants.RSUN * constants.TCODE^-2;
        ERG_PER_ECODE     = constants.MSUN * constants.CM_S_PER_VCODE^2;
        UPRES_PER_PRESCODE= constants.MSUN / (constants.RSUN * constants.TCODE^2);
    end
    
end

