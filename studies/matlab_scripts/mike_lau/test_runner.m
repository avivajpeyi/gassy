const = constants;

% create the doner star (center of CE)
star = stellarProfile('cored_adiabatic_fixedcomposition.dat'); 
star.mass = star.mass + 3.15 * const.MSUN; % AVI: why isnt this just star.mass?
star.stellarMass = star.stellarMass + 3.15 * const.MSUN; % AVI: why is this duplicated?
star.cs = sqrt(5/3 * star.pres ./ star.dens); 
star.calcGradients % init some vars

% place the accretor (at edge of star)
m2 = 0.5 * const.MSUN;
xFactor = 0.9; % Initial radius of companion as fraction of donor radius
rotFactor = 0.4; % Rotation of donor as a function of Keplerian angular velocity at surface (this decreases drag)

% set up Inspiral ODE equations
inspiral = analyticalInspiral(dt,endsep,xFactor,rotFactor,star, m2);
inspiral.setDragCoefficient(1); % Gamma = 5/3
inspiral.USE_BHLRAD = false;

% ODE Step settings
period = 2*pi * sqrt( (xFactor * star.stellarRadius)^3 / (const.G * star.stellarMass) );
dt = 0.001 * period;
endsep = 15 * const.RSUN;
 
% run ODE integrator and plot separation
inspiral.integrateEoM;
inspiral.plt('sep','red','solid','Separation')
saveas(gcf,'inspiral','png')