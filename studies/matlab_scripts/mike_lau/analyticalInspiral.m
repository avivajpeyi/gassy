classdef analyticalInspiral < handle
    % Inspiral of companion in donor envelope

    properties
        %------------------------------------------------------------------
        % Inspiral settings
        %------------------------------------------------------------------
        DT                % Size of initial timestep
        ENDSEP            % Separation below which simulation is stopped
        MAXSTEPNO = 10000 % Maximum number of timesteps
        XFACTOR           % Initial radius of companion as fraction of donor radius
        ROTFACTOR         % Rotation of donor as r function of Keplerian angular velocity at surface (this decreases drag)
        MACC              % Mass of accretor
        DONOR             % Stellar profile that the companion is inspiralling through
        iDRAGCOEF         % 0: Constant drag coefficient; 1: De+20 gamma = 5/3; 2: De+20 gamma = 4/3
                          % 3,4: Same as 1,2, but using density scale height for interpolation
                          % instead of mass ratio.
        USE_BHLRAD = true % True: Use BHL radius. False: Use HL radius.

        %------------------------------------------------------------------
        % Quantities at every step
        %------------------------------------------------------------------
        currentStep       % Current time step number
        time
        dt                % (Adaptive) time step
        pos               % Array containing x and y coordinates
        vel
        accel
        dragAccel         % Acceleration due to BHL drag
        totene            % Sum of companion KE and GPE
        mach              % Upstream mach number
        cs                % Sound speed
        dragCoff          % Multiplicative factor to Bondi-Hoyle-Lyttleton drag
        velCon            % Velocity contrast (vx, vy)
        Racc              % Accretion radius
        Hrho              % Density scale height, -dr/dlogrho
        epsilonrho        % Ratio of accretion radius to density scale height
        massIn            % Mass interior to the companion
        vKep              % Keplerian velocity, (GM/r)^1/2
        dlogP_dlogrho


        %------------------------------------------------------------------
        % Overall quantities
        %------------------------------------------------------------------
        omegaCrit         % Critical angular velocity: Angular velocity at surface of (rotating donor)
        initPeriod        % Keplerian period at initial position of companion
        pos0              % Initial position (x,y) of companion
        vel0              % Initial velocity (vx,vy) of companion

        %------------------------------------------------------------------
        % Other global variables
        %------------------------------------------------------------------
        err               % 1: err < maxerr in current timestep; 0: err > minerr in current timestep
        dtfactor = 2      % Adaptive timestepping factor
    end

    methods
        %------------------------------------------------------------------
        % Construct an instance of this class
        %------------------------------------------------------------------
        function obj = analyticalInspiral(dt,endsep,xFactor,rotFactor,donorProfile,macc)
            obj.DT           = dt;
            obj.ENDSEP       = endsep;
            obj.currentStep  = 1;
            obj.MACC         = macc;

            % Initialise arrays
            npts = obj.MAXSTEPNO; % Choose r large number as array length
            [obj.time,...
             obj.dt,...
             obj.totene,...
             obj.dragCoff,...
             obj.mach,...
             obj.cs,...
             obj.Racc,...
             obj.Hrho,...
             obj.epsilonrho,...
             obj.vKep,...
             obj.dlogP_dlogrho,...
             obj.massIn] = deal(zeros(npts,1));

            [obj.pos,...
             obj.vel,...
             obj.accel,...
             obj.velCon] = deal(zeros(npts,2));

            obj.dragAccel = zeros(npts,3); % Third columns holds sign of drag (> 0 for drag)

            obj.XFACTOR = xFactor;
            obj.ROTFACTOR = rotFactor;
            obj.DONOR = donorProfile;

            setInitPosVel(obj);
            obj.pos(1,:) = obj.pos0;
            obj.vel(1,:) = obj.vel0;
            obj.dt(1)    = obj.DT;

            getTotEne(obj);

            const = constants; % Load constants
            obj.omegaCrit = rotFactor * sqrt( const.G * donorProfile.stellarMass *...
                            donorProfile.stellarRadius^-3 );
            obj.initPeriod = 2*pi * sqrt( (xFactor * donorProfile.stellarRadius)^3 /...
                             (const.G * donorProfile.stellarMass) );
        end

        %------------------------------------------------------------------
        % Sets default initial position and velocity of companion
        %------------------------------------------------------------------
        function setInitPosVel(obj)

            [unique_rad, index] = unique(obj.DONOR.rad);
            x0 = obj.XFACTOR * obj.DONOR.stellarRadius;

            % Get mass interior to initial companion position
            init_m = interp1(unique_rad, obj.DONOR.mass(index), x0);

            % Keplerian velocity at initial companion position
            v0 = sqrt(constants.G * init_m / (obj.XFACTOR * obj.DONOR.stellarRadius));

            obj.pos0 = [x0, 0];
            obj.vel0 = [0, v0];
        end

        %------------------------------------------------------------------
        % Sets options for drag coefficient
        %------------------------------------------------------------------
        function setDragCoefficient(obj,iDRAGCOEF,Cd)
            obj.iDRAGCOEF = iDRAGCOEF;
            if iDRAGCOEF == 0
                obj.dragCoff(:) = Cd;
            end
        end

        %------------------------------------------------------------------
        % Integrate equation of motion
        %------------------------------------------------------------------
        function integrateEoM(obj)

            obj.currentStep = 1;
            while ( (norm(obj.pos(obj.currentStep,1:2)) > obj.ENDSEP) &...
                    (norm(obj.pos(obj.currentStep,1:2)) < obj.DONOR.stellarRadius) &...  % Stop if companion flies out of star
                    (obj.currentStep <= obj.MAXSTEPNO) )

                obj.eulerStep;

                %obj.RKF45step; % increment i in this call
                %disp(obj.currentStep);
                %disp(obj.dt(obj.currentStep));
            end
            obj.currentStep = obj.currentStep - 1;
            i = obj.currentStep;

            % Trim arrays
            obj.time        = obj.time(1:i);
            obj.dt          = obj.dt(1:i);
            obj.pos         = obj.pos(1:i,:);
            obj.vel         = obj.vel(1:i,:);
            obj.accel       = obj.accel(1:i,:);
            obj.dragAccel   = obj.dragAccel(1:i,:);
            obj.totene      = obj.totene(1:i);
            obj.mach        = obj.mach(1:i);
            obj.dragCoff    = obj.dragCoff(1:i);
            obj.velCon      = obj.velCon(1:i,:);
            obj.cs          = obj.cs(1:i,:);
            obj.Racc        = obj.Racc(1:i,:);
            obj.Hrho        = obj.Hrho(1:i,:);
            obj.epsilonrho  = obj.epsilonrho(1:i,:);
            obj.vKep        = obj.vKep(1:i,:);
            obj.dlogP_dlogrho = obj.dlogP_dlogrho(1:i,:);

            obj.getTotEne;      % Calculate total energy at every timestep
            obj.getvKep;        % Get Keplerian velocity at every timestep
            obj.getdlogPdlogrho; % Get dlogP/dlogrho at every timestep

        end

        %-------------------------------------------------------------------
        % Calculates dlogP/dlogrho at every timestep
        %-------------------------------------------------------------------
        function getdlogPdlogrho(obj)
           companionRadiusArray = sqrt(sum(obj.pos.^2, 2));
           [rArray, index] = unique(obj.DONOR.dlogPdlogrho(:,1));
           obj.dlogP_dlogrho = interp1(rArray, obj.DONOR.dlogPdlogrho(index,2), companionRadiusArray);
        end

        %-------------------------------------------------------------------
        % Calculates local Keplerian velocity at every timestep
        %-------------------------------------------------------------------
        function getvKep(obj)
           companionRadiusArray = sqrt(sum(obj.pos.^2, 2));
           [rArray, index] = unique(obj.DONOR.rad);
           interiorMass = interp1(rArray, obj.DONOR.mass(index), companionRadiusArray);
           interiorMass(isnan(interiorMass)) = obj.DONOR.stellarMass; % Interpolation fails if companion is outside donor

           obj.vKep = sqrt(constants.G * (interiorMass + obj.DONOR.stellarMass)...
                      ./ companionRadiusArray);
        end

        %-------------------------------------------------------------------
        % Calculates total specific energy of companion. Note: Energy
        % will be infinite for uncalculated timesteps.
        %-------------------------------------------------------------------
        function getTotEne(obj)
           kin = 0.5 * sum(obj.vel.^2, 2);

           companionRadiusArray = sqrt(sum(obj.pos.^2, 2));
           [rArray, index] = unique(obj.DONOR.rad);
           interiorMass = interp1(rArray, obj.DONOR.mass(index), companionRadiusArray);
           interiorMass(isnan(interiorMass)) = obj.DONOR.stellarMass; % Interpolation fails if companion is outside donor

           pot = - constants.G * interiorMass ./ companionRadiusArray;

           obj.totene = kin + pot;
        end


        %------------------------------------------------------------------
        % Returns total acceleration (gravity plus BHL drag) without
        % writing properties at current timestep. Used for intermediate
        % steps in RKF45step
        %------------------------------------------------------------------
        function derivs = getAccel(obj,coords,iwrite)
            % derivs: 1x4 array containing time derivative of coords
            arguments
                obj
                coords              % 1x4 array containing x, y, vx, vy
                iwrite = false      % true to write quantities to object properties
            end

            pos = coords(1:2);
            vel = coords(3:4);
            accel = analyticalInspiral.getGravForce(obj.DONOR,pos) +...
                    getBHLdrag(obj,pos,vel,obj.DONOR,obj.omegaCrit,obj.MACC,obj.iDRAGCOEF,obj.USE_BHLRAD,iwrite);
            derivs = [vel(1), vel(2), accel(1), accel(2)];

        end


        %------------------------------------------------------------------
        % Calls calcBHLdrag to calculate acceleration [Fx,Fy] due to
        % Bondi-Hoyle-Lyttleton drag as r function of current position
        % [x,y] and velocity [vx,vy], and the stellar profile
        %------------------------------------------------------------------
        function accel = getBHLdrag(obj,pos,vel,star,omega,macc,iDRAGCOEF,USE_BHLRAD,iwrite)
            % omega: Angular velocity of star
            % macc: Companion mass
            % iwrite: True if writing quantities to object property arrays

            companionRad = sqrt(sum(pos.^2,2));

            if companionRad >= star.stellarRadius
                Fx = 0;
                Fy = 0;
                signdrag = 0;
            else
                % Interpolate density and sound speed
                [rArray, index] = unique(star.rad);
                rho = interp1(rArray, star.dens(index), companionRad);
                cs = interp1(rArray, star.cs(index), companionRad);

                % Calculate velocity contrast
                vcon(1) = vel(1) + omega * pos(2);
                vcon(2) = vel(2) - omega * pos(1);
                vcon2 = sum(vcon.^2);

                % Get Mach number
                interiorMass = interp1(rArray, star.mass(index), companionRad);
                mach = sqrt(vcon2) / cs;

                Racc = 2 * constants.G * macc / vcon2; % Should really be /(cs^2 + vcon2)
                Hrho = interp1(star.Hrho_vs_r(:,1), star.Hrho_vs_r(:,2), companionRad);
                epsilonrho = Racc / Hrho;

                % Calculate drag coefficient if using De+20 interpolation
                if (iDRAGCOEF == 1 || iDRAGCOEF == 2)
                    qr = macc / interiorMass;
                    dragCoef = CEwindTunnel.getDragFromQandMach(qr, mach, iDRAGCOEF);
                elseif (iDRAGCOEF == 3 || iDRAGCOEF == 4)
                    if iDRAGCOEF == 3; igamma = 1; end
                    if iDRAGCOEF == 4; igamma = 2; end
                    dragCoef = CEwindTunnel.getDragFromEpsandMach(mach, epsilonrho, igamma);
                elseif iDRAGCOEF == 0
                    dragCoef = obj.dragCoff(obj.currentStep);
                end


                if USE_BHLRAD % Use 1/(v^2 + c_s^2)^3/2
                    Fx = dragCoef * 4*pi * constants.G^2 * macc * rho *...
                        ( vcon2 + cs^2 )^-1.5 * -vcon(1);
                    Fy = dragCoef * 4*pi * constants.G^2 * macc * rho *...
                        ( vcon2 + cs^2 )^-1.5 * -vcon(2);
                else % Use 1/v^3
                    Fx = dragCoef * 4*pi * constants.G^2 * macc * rho *...
                        vcon2^-1.5 * -vcon(1);
                    Fy = dragCoef * 4*pi * constants.G^2 * macc * rho *...
                        vcon2^-1.5 * -vcon(2);
                end

                signdrag = sign( sum( vcon .* vel ) ); % > 0 for r drag
            end
            accel = [Fx, Fy];

            if iwrite
                i = obj.currentStep;
                obj.velCon(i,:)     = vcon;
                obj.cs(i)           = cs;
                obj.massIn(i)       = interiorMass;
                obj.mach(i)         = mach;
                obj.accel(i,:)      = obj.accel(i,:) + [Fx, Fy];
                obj.dragAccel(i,:)  = [Fx, Fy, signdrag];
                obj.Racc(i)         = Racc;
                obj.Hrho(i)         = Hrho;
                obj.epsilonrho(i)   = epsilonrho;
                obj.dragCoff(i)     = dragCoef;
            end

        end


        %------------------------------------------------------------------
        % Step with semi-implicit Euler method and increment step
        %------------------------------------------------------------------
        function eulerStep(obj)
            i = obj.currentStep;
            iwrite = true;
            accel = analyticalInspiral.getGravForce(obj.DONOR,obj.pos(i,:)) +...
                    getBHLdrag(obj,obj.pos(i,:),obj.vel(i,:),obj.DONOR,...
                    obj.omegaCrit,obj.MACC,obj.iDRAGCOEF,obj.USE_BHLRAD,iwrite);

            obj.time(i+1)  = obj.time(i)  + obj.DT;
            obj.vel(i+1,:) = obj.vel(i,:) + accel * obj.DT;
            obj.pos(i+1,:) = obj.pos(i,:) + obj.vel(i+1,:) * obj.DT;


            obj.currentStep = i + 1;
        end


        %------------------------------------------------------------------
        % Step with adaptive Runge-Kutta-Fehlberg 45 method
        %------------------------------------------------------------------
        function RKF45step(obj)
            i = obj.currentStep;
            r = obj.pos(i,:);      % (xi, yi)
            v = obj.vel(i,:);      % (vxi, vyi)
            dt_guess = obj.dt(i);  % Guess for timestep

            % Coefficients in RKF45 method
            b21                   = 0.25;
            [b31,b32]             = deal(3/32,9/32);
            [b41,b42,b43]         = deal(1932/2197,-7200/2197,7296/2197);
            [b51,b52,b53,b54]     = deal(439/216,-8,3680/513,-845/4104);
            [b61,b62,b63,b64,b65] = deal(-8/27,2,-3544/2565,1859/4104,-11/40);
            [c1,c3,c4,c5]         = deal(25/216,1408/2565,2197/4104,-0.20);
            [d1,d3,d4,d5,d6]      = deal(16/135,6656/12825,28561/56430,-0.18,2/55);

            k1 = obj.getAccel([r,v], true); % write quantities for current timestep
            k2 = obj.getAccel([r,v] + dt_guess * b21*k1);
            k3 = obj.getAccel([r,v] + dt_guess * (b31*k1 + b32*k2));
            k4 = obj.getAccel([r,v] + dt_guess * (b41*k1 + b42*k2 + b43*k3));
            k5 = obj.getAccel([r,v] + dt_guess * (b51*k1 + b52*k2 + b53*k3 + b54*k4));
            k6 = obj.getAccel([r,v] + dt_guess * (b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5));

            poorStep = dt_guess * (c1*k1 + c3*k3 + c4*k4 + c5*k5);         % Fourth order step
            goodStep = dt_guess * (d1*k1 + d3*k3 + d4*k4 + d5*k5 + d6*k6); % Fifth order step
            err = max(abs(poorStep - goodStep) ./ poorStep);

            % Check tolerance
            maxerr = 0.001;  % Criterion for halving timestep
            minerr = 0.0005; % Criterion for doubling timestep
            dtmin  = 1e-15 * obj.initPeriod;
            dtmax  = 1e-3 * obj.initPeriod;

            if ( (err > maxerr) && (dt_guess > dtmin) )
                % Error is too large; decrease step size and try again
                if obj.err == 0
                % If in the previous timestep, error was too large, but in
                % current timestep, error is too small, obj.dtfactor is too
                % large and so we halve it.
                    obj.dtfactor = obj.dtfactor / 2;
                end
                obj.dt(i) = dt_guess / (1 + obj.dtfactor);
                obj.err = 1;
            elseif ( (err < minerr) && (dt_guess < dtmax) )
                % Error is small enough so increase timestep and try again
                if obj.err == 1
                    obj.dtfactor = obj.dtfactor / 2;
                end
                obj.dt(i) = dt_guess * (1 + obj.dtfactor);
                obj.err = 0;
%             elseif any(abs(goodStep(3:4)./obj.vel(i,:)) > 0.01)
%                 obj.dt(i) = dt_guess / (1 + obj.dtfactor);
%                 obj.err = 1;
            else % Accept step size and step forward
                % Write acceleration and other quantities
                obj.dtfactor = 2; % Reset default dtfactor to 2
                obj.currentStep = i + 1;
                obj.dt(i+1)     = dt_guess;
                obj.time(i+1)   = obj.time(i)  + dt_guess;
                obj.pos(i+1,:)  = obj.pos(i,:) + goodStep(1:2);
                obj.vel(i+1,:)  = obj.vel(i,:) + goodStep(3:4);
            end

        end

        %------------------------------------------------------------------
        % Plot against time / yrs
        %------------------------------------------------------------------
        function myplot = plt(obj,quantity,color,linestyle,legendText)
            % Possible strings that quantity can take:
            %'sep' : Radius of companion
            %'drag': Acceleration due to drag force
            %'mach': Mach number
            %'Cd'  : Drag coefficient
            %'kin' : Specific kinetic energy of companion
            %'gpe' : Magnitude of specific gravitational potential energy of companion
            %'ene' : Total energy of companion (kinetic + potential)

            % First calculate the plot range. Plot companion down to small
            % fraction of inspiral radius

            sep = sqrt(obj.pos(:,1).^2 + obj.pos(:,2).^2);

            switch quantity
                case 'sep'
                    yplot = sep / constants.RSUN;
                case 'drag'
                    drag = obj.dragAccel(:,3) .* sqrt(sum( obj.dragAccel(:,1:2).^2, 2));
                    yplot = drag;
                case 'mach'
                    yplot = obj.mach;
                case 'Cd'
                    yplot = obj.dragCoff;
                case 'kin'
                    kin = 0.5 * sqrt(obj.vel(:,1).^2 + obj.vel(:,2).^2);
                    yplot = kin;
                case 'gpe'
                    companionRadiusArray = sqrt(sum(obj.pos.^2, 2));
                    [rArray, index] = unique(obj.DONOR.rad);
                    interiorMass = interp1(rArray, donorProfile.mass(index), companionRadiusArray);
                    interiorMass(isnan(interiorMass)) = donorProfile.stellarMass; % Interpolation fails if companion is outside donor

                    gpe = - constants.G * interiorMass ./ companionRadiusArray;
                    yplot = gpe;
                case 'ene'
                    yplot = obj.totene;
                case 'velCon'
                    vcon_mag = sqrt(sum(transpose(obj.velCon.^2), 2));
                    vcon_sign = obj.dragAccel(:,3);
                    yplot = vcon_mag
                    %yplot = vcon_sign .* vcon_mag;
                case 'velCon/vKep'
                    companionRadiusArray = sqrt(sum(obj.pos.^2, 2));
                    [rArray, index] = unique(obj.DONOR.rad);
                    interiorMass = interp1(rArray, obj.DONOR.mass(index), companionRadiusArray);
                    interiorMass(isnan(interiorMass)) = obj.DONOR.stellarMass; % Interpolation fails if companion is outside donor

                    vKep = sqrt(constants.G * (interiorMass + obj.MACC) ./ companionRadiusArray);
                    vcon = sqrt(sum(obj.velCon.^2, 2));

                    vcon_sign = obj.dragAccel(:,3);
                    yplot = vcon_sign .* vcon ./ vKep ;
                case 'vel/vKep'
                    companionRadiusArray = sqrt(sum(obj.pos.^2, 2));
                    [rArray, index] = unique(obj.DONOR.rad);
                    interiorMass = interp1(rArray, obj.DONOR.mass(index), companionRadiusArray);
                    interiorMass(isnan(interiorMass)) = obj.DONOR.stellarMass; % Interpolation fails if companion is outside donor

                    vKep = sqrt(constants.G * (interiorMass + obj.MACC) ./ companionRadiusArray);
                    vel = sqrt(sum((obj.velCon.^2), 2));

                    yplot = vel ./ vKep ;
                case 'racc/sep'
                    companionRadiusArray = sqrt(sum(obj.pos.^2, 2));
                    [rArray, index] = unique(obj.DONOR.rad);
                    csArray = interp1(rArray, obj.DONOR.cs(index), companionRadiusArray);

                    vcon2 = sum(obj.velCon.^2, 2);

                    yplot = 2 * constants.G * obj.MACC ./ (csArray.^2 + vcon2);
                    yplot = yplot ./ companionRadiusArray;

                case 'bondirad'
                    companionRadiusArray = sqrt(sum(obj.pos.^2, 2));
                    [rArray, index] = unique(obj.DONOR.rad);
                    csArray = interp1(rArray, obj.DONOR.cs(index), companionRadiusArray);

                    yplot = 2 * constants.G * obj.MACC ./ csArray.^2;
                    yplot = yplot ./ companionRadiusArray;

                otherwise
                    warning('Unexpected plot type.')
            end

            myplot = plot(obj.time / constants.SEC_PER_YR,...
                          yplot,linestyle,'Color',color,...
                          'DisplayName',legendText);

        end

        %------------------------------------------------------------------
        % Plot trajectory
        %------------------------------------------------------------------
        function myplot = plttraj(obj,opt)
            % First calculate the plot range. Plot companion down to small
            % fraction of inspiral radius
            arguments
               obj
               opt.color = 'b'
               opt.linestyle = '-'
               opt.legendText = 'Trajectory'
            end

            myplot = plot(obj.pos(:,1), obj.pos(:,2),...
                          opt.linestyle,'Color',opt.color,...
                          'DisplayName',opt.legendText); hold on
            plotcircle(0,0,obj.DONOR.stellarRadius);
            axis equal tight
        end

        %------------------------------------------------------------------
        % Plot trajectory in epsilon-Mach space
        %------------------------------------------------------------------
        function pltEpsMachTraj(obj)
           arguments
               obj
           end

           plot(obj.epsilonrho(1:obj.currentStep-1), obj.mach(1:obj.currentStep-1)); hold on
           startMarker = plot(obj.epsilonrho(1), obj.mach(1), '*');
           set(get(get(startMarker,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Do not create legend entry

           epsilonRhoGrid = linspace(0,10,300);
           machGrid = sqrt(2*epsilonRhoGrid);
           q1_Line = plot(epsilonRhoGrid,machGrid,'--k','LineWidth',1);
           text(0.75*epsilonRhoGrid(end), 0.95*machGrid(end), "$q=1$", 'interpreter','latex','fontsize',12);
           set(get(get(q1_Line,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Do not create legend entry

           ylabel("$\mathcal{M}_\infty$", "interpreter", "latex");
           xlabel("$\epsilon_\rho$", "interpreter", "latex");
        end
    end

    methods (Static)
        %------------------------------------------------------------------
        % Returns gravitational acceleration (ax, ay) given r
        % stellarProfile object and an array of positions, pos = [x,y]
        %------------------------------------------------------------------
        function accel = getGravForce(star, pos)
            % pos: nx2 array
            % star: stellarProfile object
            % Fx, Fy: nx2 array of gravitational acceleration components

            companionRad = sqrt(sum(pos.^2,2));

            if companionRad >= star.stellarRadius
                m = star.stellarMass;
            else
                % Interpolate enclosed mass
                [rArray, index] = unique(star.rad);
                m = interp1(rArray, star.mass(index), companionRad); % m has same dimensions as companionRad
            end

            Fx = - constants.G * m .* pos(:,1) ./ companionRad^3;
            Fy = - constants.G * m .* pos(:,2) ./ companionRad^3;
            accel = [Fx, Fy];

        end

    end
end

