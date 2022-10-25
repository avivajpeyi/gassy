classdef stellarProfile < handle
    
    properties
        %------------------------------------------------------------------
        % Data options
        %------------------------------------------------------------------
        filePath            % (str) path to data file
        hasCompositionFlag  % (logical) Flag to indicate if data file contains H and He fractions
        hasSoundSpeedFlag   % (logical) Flag to indicate if data file contains sound speed profile
        description         % (str) Description, used as legend text
        addedCoreMass       % (logical) Flag to indicate if a core mass has been added to the mass profile
        
        %------------------------------------------------------------------
        % Stellar profile (numeric arrays, cgs units)
        %------------------------------------------------------------------
        mass    % Mass coordinate
        pres    % Pressure
        temp    % Temperature
        rad     % Radius
        dens    % Density
        ene     % Specific internal energy
        XX      % Hydrogen mass fraction
        YY      % Helium mass fraction
        cs      % Sound speed
        mu      % Mean molecular weight
        
        %------------------------------------------------------------------
        % Derived profiles (numeric arrays, cgs units)
        %------------------------------------------------------------------
        entropy             % A quantity that is equal to the gas + radiation
                            % mass-specific entropy, up to some constant
        drhodr_vs_r         % A Nx2 array with the first column containing r 
                            % (this is not the same as rad as it doens't contain
                            % repeated elements), and the second column
                            % containing -drho/dr.
        dPdr_vs_r           % A Nx2 array with the first column containing r 
                            % (this is not the same as rad as it doens't contain
                            % repeated elements), and the second column
                            % containing -dP/dr.
        dlogPdlogrho        % A Nx2 array with the first column containing r 
                            % (this is not the same as rad as it doens't contain
                            % repeated elements), and the second column
                            % containing dlogP/dlogrho.
        Hrho_vs_r           % A Nx2 array with the first column containing r 
                            % (this is not the same as rad as it doens't contain
                            % repeated elements), and the second column
                            % containing the density scale height, -dr/dlnrho.
        Ebind_grav_from_surf% Magnitude of gravitational potential energy integrated from surface
        Ebind_int_from_surf % Binding energy (including internal energy) integrated from surface
                      
        %------------------------------------------------------------------
        % Global stellar properties
        %------------------------------------------------------------------
        stellarMass
        stellarRadius
        momentOfInertia
        
    end
    
    
    
    methods
        %------------------------------------------------------------------
        % Construct an instance of this class
        %------------------------------------------------------------------
        function obj = stellarProfile(filePath,opt)
            % (assuming quantities are sorted from centre to surface)
            
            arguments
               filePath
               opt.description        = filePath
               opt.hasCompositionFlag = false
               opt.hasSoundSpeedFlag  = false
               opt.initEmptyClass     = false
            end
            
            obj.filePath = filePath;
            obj.description = opt.description;

            if (opt.initEmptyClass); return; end

            obj.addedCoreMass = false;
            obj.hasCompositionFlag = opt.hasCompositionFlag;
            obj.hasSoundSpeedFlag = opt.hasSoundSpeedFlag;
            
            starDataArray = obj.importStellarProfile;
            obj.mass  = starDataArray(:,1);
            obj.pres  = starDataArray(:,2);
            obj.temp  = starDataArray(:,3);
            obj.rad   = starDataArray(:,4);
            obj.dens  = starDataArray(:,5);
            obj.ene   = starDataArray(:,6);
            
            if obj.hasCompositionFlag
                obj.XX    = starDataArray(:,7);
                obj.YY    = starDataArray(:,8);
            end
            if obj.hasSoundSpeedFlag
                obj.cs    = starDataArray(:,9);
            end
            
            obj.stellarMass   = max(obj.mass);
            obj.stellarRadius = max(obj.rad);
            
            % Additional calculations
%             obj.calcGradients;
            obj.calcMOI;
        end
        
        %------------------------------------------------------------------
        % If star is a softened star, add back core mass onto mass
        % coordinate and total mass
        %------------------------------------------------------------------
        function obj = addCoreMass(obj,mcore)
            if ~obj.addedCoreMass
                obj.stellarMass = obj.stellarMass + mcore;
                obj.mass = obj.mass + mcore;
                obj.addedCoreMass = true;
            end
        end
            
        %------------------------------------------------------------------
        % Plot quantity against radius
        %------------------------------------------------------------------
        function myplot = plt(obj,quantity,opt)
            
            % Specify default value for optional arguments
            arguments
               obj
               quantity
               opt.leg = obj.description;
               opt.imass = false;
               opt.sunUnits = false; % If true, plot m in solar masses and r in solar radii instead
               opt.logx = true;
               opt.logy = true;
               opt.yscale = 1;
            end
            
            u = constants;
            switch quantity
                case 'mass'
                    if opt.sunUnits
                        yplot = obj.mass / u.MSUN;
                        myylab = 'm / M${}_\odot$';
                    else
                        yplot = obj.mass;
                        myylab = 'm / g';
                    end
                case 'pres'
                    yplot = obj.pres;
                    myylab = 'P / dyn cm${}^{-2}$';
                case 'ene'
                    yplot = obj.ene;
                    myylab = 'u / erg g${}^{-1}$';
                case 'dens'
                    yplot = obj.dens;
                    myylab = '$\rho$ / g cm${}^{-3}$';
                case 'temp'
                    yplot = obj.temp;
                    myylab = 'T / K';
                case 'XX'
                    yplot = obj.XX;
                    myylab = 'X';
                case 'YY'
                    yplot = obj.YY;
                    myylab = 'Y';
                case 'mu'
                    if isempty(obj.mu); error("Mean molecular weight not calculated"); end
                    yplot = obj.mu;
                    myylab = '$\mu$';
                case 'cs'
                    yplot = obj.cs * 1.e-5;
                    myylab = '$c_s$ / km s${}^{-1}$';
                case 'entropy'
                    if isempty(obj.entropy); error("Entropy not calculated"); end
                    yplot = obj.entropy;
                    myylab = 'S / erg K${}^{-1}$ g${}^{-1}$';
                case 'Ebind'
                    yplot = abs(obj.ene - u.G * obj.mass ./ obj.rad);
                    myylab = '$|E_\mathrm{bind}|$ / erg g${}^{-1}$';
                case 'Ebind_grav_from_surf'
                    yplot = obj.Ebind_grav_from_surf / 1e48;
                    myylab = '$|E_\mathrm{bind}(>r)|$ / $10^{48}$ erg';
                case 'Ebind_int_from_surf'
                    yplot = obj.Ebind_int_from_surf / 1e48;
                    myylab = '$|E_\mathrm{bind}(>r)|$ / $10^{48}$ erg';
                otherwise
                    warning('Unexpected plot type.')
            end
            
            if opt.imass
                if opt.sunUnits
                        xplot = obj.mass / u.MSUN;
                        myxlab = '$m$ / M${}_\odot$';
                else
                    xplot = obj.mass;
                    myxlab = '$m$ / g';
                end
            else
                if opt.sunUnits
                    xplot = obj.rad / u.RSUN;
                    myxlab = '$r$ / R${}_\odot$';
                else
                    xplot = obj.rad;
                    myxlab = '$r$ / cm';
                end
            end
            
            if quantity == "Ebind_int"
               xplot = xplot(2:end); 
            end
            
            myplot = plot(xplot, yplot * opt.yscale,'DisplayName',opt.leg);
            if opt.logx; set(gca,'Xscale','log'); end
            if opt.logy; set(gca,'Yscale','log'); end
            xlabel(myxlab, 'interpreter', 'latex');
            ylabel(myylab, 'interpreter', 'latex');
                      
        end
       
        %------------------------------------------------------------------
        % Print summary about binding energy. Provide mcore in g
        %------------------------------------------------------------------
        function printBindingEneSummary(obj,mcore)
            [bindingEneG, bindingEneB,lambdaG, lambdaB] = obj.calcEnvBindingEne(mcore);
            intEne = obj.calcTotalIntEne;
            fprintf('\nE_bind,G = %e erg\n', bindingEneG);
            fprintf('E_bind,B = %e erg\n', bindingEneB);
            fprintf('E_int = %e erg\n', intEne);
            fprintf('lambdaG = %e\n', lambdaG);
            fprintf('lambdaB = %e\n', lambdaB);
        end
        
        %------------------------------------------------------------------
        % Calculate alpha inferred from simulation given final separation
        %------------------------------------------------------------------
        function [alphaB, alphaG] = calcAlpha(obj,sep_final,coreMass,m2)
            u = constants;
            Delta_Eorb = 0.5 * u.G * m2 * abs(obj.stellarMass / obj.stellarRadius - coreMass / sep_final);
            [bindingEneG,bindingEneB,~,~] = calcEnvBindingEne(obj,coreMass);
            alphaB = bindingEneB / Delta_Eorb;
            alphaG = bindingEneG / Delta_Eorb;
            
            fprintf('alphaG = %e\n', alphaG);
            fprintf('alphaB = %e\n', alphaB);
        end
        
        %------------------------------------------------------------------
        % Calculate total binding energy integrated from the surface
        % Returns cell array Ebind containing
        % * Ebind(1): Purely gravitational binding energy
        % * Ebind(2): Binding energy including internal energy
        %------------------------------------------------------------------
        function calcIntegratedEbind(obj)
            u = constants;
            
            % Initialise arrays
            Ebind_grav_from_centre = zeros(size(obj.mass));
            Ebind_int_from_centre  = zeros(size(obj.mass));
            eint_from_centre       = zeros(size(obj.mass));
            
            integrand = - u.G * obj.mass ./ obj.rad; % -Gm/r
            Ebind_grav_from_centre(2:end) = cumtrapz(obj.mass(2:end), integrand(2:end)); % Central Ebind, Ebind(1), is infinity due to division by r=0
            eint_from_centre(2:end)       = cumtrapz(obj.mass(2:end), obj.ene(2:end));
            Ebind_int_from_centre(2:end)  = Ebind_grav_from_centre(2:end) + eint_from_centre(2:end);
                        
            Ebind_grav_from_centre = abs(Ebind_grav_from_centre);
            Ebind_int_from_centre  = abs(Ebind_int_from_centre);
            
            % Binding energy interated from surface
            obj.Ebind_grav_from_surf = Ebind_grav_from_centre(end) - Ebind_grav_from_centre(1:end);
            obj.Ebind_int_from_surf = Ebind_int_from_centre(end) - Ebind_int_from_centre(1:end);
        end
        
        %------------------------------------------------------------------
        % Calculate entropy (see properties for description) given mean
        % molecular weight (mu)
        %------------------------------------------------------------------
        function calcEntropy(obj,mu,iExcludeRadiation)
            arguments
               obj
               mu                        % Can be either a scalar (if mu = constant in star)
                                         % or the full mean molecular
                                         % weight profile
               iExcludeRadiation = false % True if only gas entropy wanted
            end
            
            if length(mu) == 1
                mu = mu * ones( size(obj.rad) );
            end
            
            u = constants;
            Sgas = u.KB / u.PROTONMASS ./ mu .* log( obj.temp.^1.5 ./ obj.dens );
            Srad = 4/3 * u.RADCONST * obj.temp.^3 ./ obj.dens;

            if iExcludeRadiation
                obj.entropy = Sgas;
            else
                obj.entropy = Sgas + Srad;
            end
                      
        end
        
        %------------------------------------------------------------------
        % Calculate mean molecular weight profile using pressure,
        % temperature, and density
        %------------------------------------------------------------------
        function calcMu(obj,iExcludeRadiation)
            arguments
               obj
               iExcludeRadiation = false % True if excluding radiation pressure
            end
            
            u = constants;
            
            if iExcludeRadiation
                pgas = obj.pres;
            else
                pgas = obj.pres - 1/3 * u.RADCONST * obj.temp.^4;
            end
            
            obj.mu = obj.dens * u.KB .* obj.temp ./ ( u.PROTONMASS * pgas );
                      
        end
        
        %------------------------------------------------------------------
        % Calculate envelope gravitational binding energy and lambda parameter
        %------------------------------------------------------------------
        function [bindingEneG,bindingEneB,lambdaG,lambdaB] = calcEnvBindingEne(obj,coreMass,coreEnvBoundary)
            arguments
                obj
                coreMass
                coreEnvBoundary = nan       % Supply core-envelope boundary.
            end
            
            envMass = obj.stellarMass - coreMass;
            
            if isnan(coreEnvBoundary) % Integrate through entire star
                if (obj.rad(1) == 0.) 
                    startidx = 2; % If r_0 = 0, integrate from r_1 instead to avoid singularity
                else
                    startidx = 1;
                end
            else
                [~, startidx] = min( abs( obj.rad - coreEnvBoundary ) );
            end
            
            integrand = obj.mass ./ obj.rad;
            bindingEneG = abs(constants.G * trapz(obj.mass(startidx:end), integrand(startidx:end)));
            bindingEneB = bindingEneG - calcTotalIntEne(obj,coreEnvBoundary);
            lambdaG = constants.G * coreMass * envMass / (obj.stellarRadius * bindingEneG);
            lambdaB = constants.G * coreMass * envMass / (obj.stellarRadius * bindingEneB);
        end
        
        %------------------------------------------------------------------
        % Integrate total internal energy (of envelope)
        %------------------------------------------------------------------
        function intEne = calcTotalIntEne(obj,coreEnvBoundary)
            arguments
                obj
                coreEnvBoundary = nan       % Supply core-envelope boundary.
            end
                        
            if isnan(coreEnvBoundary) % Integrate through entire star
                startidx = length(obj.rad);
            else
                [~, startidx] = min( abs( obj.rad - coreEnvBoundary ) );
            end
            
            intEne = abs(trapz(obj.mass(startidx:end), obj.ene(startidx:end))); % Abs value in case quantities are sorted from surface to centre.
        end
        
        
        %------------------------------------------------------------------
        % Calculate CE final separation from alpha-lambda prescription
        %------------------------------------------------------------------
        function [finalSep_G, finalSep_B] = alphaLambdaFinalSep(obj,Macc,coreMass,coreEnvBoundary,alpha,initSep)
            % finalSep_G: Assuming only gravitational binding energy
            % finalSep_B: Assuming both gravitational binding energy and
            %             internal energy
            arguments
                obj
                Macc                        % Mass of accretor
                coreMass                    % Donor core mass
                coreEnvBoundary = nan       % Supply core-envelope boundary. By
                                            % default, integrate entire star
                alpha           = 1.        % Common-envelope efficiency
                initSep = obj.stellarRadius % Initial separation
            end
            
            % Energy to unbind
            [E_G, E_B, ~, ~] = obj.calcEnvBindingEne(coreMass,coreEnvBoundary);
            
            finalSep_G = initSep / (2*E_G * initSep / (alpha...
                       * constants.G * Macc * coreMass) + obj.stellarMass / coreMass);
            finalSep_B = initSep / (2*E_B * initSep / (alpha...
                       * constants.G * Macc * coreMass) + obj.stellarMass / coreMass);
            
            fprintf('\nFinal separation is %e Rsun\n', finalSep_G / constants.RSUN);
            fprintf('Final separation with internal energy is %e Rsun\n', finalSep_B / constants.RSUN);
        end
        
        
        %------------------------------------------------------------------
        % Calculate density and pressure gradients, assumming hydrostatic
        % equilibrium
        %------------------------------------------------------------------
        function calcGradients(obj)
            
            r    = obj.rad(2:end);   % Avoid r = 0 point
            rho  = obj.dens(2:end);
            P = obj.pres(2:end);
            
            % Density gradient
            [xData, yData] = prepareCurveData( log10(r), log10(rho) );
            ft = fittype( 'smoothingspline' );
            opts = fitoptions( 'Method', 'SmoothingSpline' );
            opts.SmoothingParam = 0.999999999999687;
            [logrho_vs_logr_fit, gof] = fit(xData, yData, ft, opts);
            dlogrho_dlogr = differentiate( logrho_vs_logr_fit, log10(r) );
            obj.drhodr_vs_r(:,1) = r;
            obj.drhodr_vs_r(:,2) = rho ./ r .* dlogrho_dlogr;
            
            % Density scale height
            obj.Hrho_vs_r(:,1) = r;
            obj.Hrho_vs_r(:,2) = -r ./ dlogrho_dlogr;
                                                           
            % Pressure gradient
            [xData, yData] = prepareCurveData( log10(r), log10(P) );
            ft = fittype( 'smoothingspline' );
            [logP_vs_logr_fit, gof] = fit(xData, yData, ft);
            dlogP_dlogr = differentiate( logP_vs_logr_fit, log10(r) );
            obj.dPdr_vs_r(:,1) = r;
            obj.dPdr_vs_r(:,2) = P ./ r .* dlogP_dlogr;
                 
            % Polytropic index of profile
            obj.dlogPdlogrho(:,1) = r;
            obj.dlogPdlogrho(:,2) = dlogP_dlogr ./ dlogrho_dlogr;
                            
        end
  
        
        %------------------------------------------------------------------
        % Calculate moment of inertia of envelope
        %------------------------------------------------------------------
        function MoI = calcMOI(obj, opt)
            arguments
                obj
                opt.rmin = 0. % Radial coordinate to integrate from (default is zero, integrating the MoI for the entire star)
            end
            
            if opt.rmin == 0.
                ind1 = 1;
            else
                [~,ind1] = min(abs(obj.rad - opt.rmin));
            end

            % Check whether quantities are sorted from surface to centre
            if obj.rad(1) < obj.rad(2)
                dom = ind1:size(obj.rad);  % Integration domain
            else
                dom = 1:ind1;
            end

            MoI = 8*pi/3 * trapz(obj.rad(dom), obj.dens(dom) .* obj.rad(dom).^4);
            obj.momentOfInertia = MoI;
            
        end
        
        %------------------------------------------------------------------
        % Read data file
        %------------------------------------------------------------------
        function starDataArray = importStellarProfile(obj)
            filePath = obj.filePath;
            startRow = 2;

            if (obj.hasSoundSpeedFlag && obj.hasCompositionFlag)
                formatSpec = '%13f%15f%15f%15f%15f%15f%15f%15f%f%[^\n\r]';
            elseif (obj.hasCompositionFlag && ~obj.hasSoundSpeedFlag)
                formatSpec = '%13f%15f%15f%15f%15f%15f%15f%f%[^\n\r]';
            else
                formatSpec = '%13f%15f%15f%15f%15f%f%[^\n\r]';
            end
            fileID = fopen(filePath,'r');

            dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', ...
                '', 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', ...
                false, 'EndOfLine', '\r\n');
            
            fclose(fileID);
            starDataArray = [dataArray{1:end-1}];
        end
    end
    
    methods(Static)
        %------------------------------------------------------------------
        % Generate plot labels
        %------------------------------------------------------------------        
        function label = makeLabel(quantity,ax)

            % Specify default value for optional arguments
            arguments
               quantity
               ax = 'y'             % Either 'x' or 'y' 
            end

            switch quantity
                case 'rad'
                    labstr = '$r$ / cm';
                case 'mass'
                    labstr = '$m$ / g';
                case 'pres'
                    labstr = '$P$ / dyn cm${}^{-2}$';
                case 'ene'
                    labstr = 'u / erg g${}^{-1}$';
                case 'dens'
                    labstr = '$\rho$ / g cm${}^{-3}$';
                case 'temp'
                    labstr = 'T / K';
                case 'XX'
                    labstr = '$X$ mass fraction';
                case 'YY'
                    labstr = '$Y$ mass fraction';
                case 'cs'
                    labstr = '$c_s$ / cms${}^{-1}$';
                case 'entropy'
                    labstr = '$S$ / erg K${}^{-1}$ g${}^{-1}$';
                otherwise
                    warning('Unexpected plot type.')
            end

            if strcmp(ax,'x')
                label = xlabel(labstr,'Interpreter','latex');
            elseif strcmp(ax,'y')
                label = ylabel(labstr,'Interpreter','latex');
            else
                warning('"ax" can only take value "x" or "y"')
            end

        end 


    end
end

