classdef CEwindTunnel
    % Bunch of methods related to calculation of interpolated drag
    % coefficients from Common-Envelope Wind Tunnel simulations
    
    methods(Static)
        %------------------------------------------------------------------
        % Calculate drag coefficient according to the fits by De+2020 to
        % their common-envelope wind tunnel simulations, given mass ratio q
        % and mach number.
        %
        % NOTE: Upstream Mach number and epsilon_rho are the ?true? input
        % parameters that directly specify a wind tunnel simulation, and
        % the interpolation formulae should really be a function of these
        % two parameters. De+20 wanted to parametrise their results in
        % terms of q instead of epsilon_rho, and so used the hydrostatic
        % equilibrium (Eq. 10) expression to recast (epsilon_rho, mach)
        % into (q, mach). 
        %------------------------------------------------------------------
        function Cd = getDragFromQandMach(q,mach,igamma)
            % q:      Ratio of companion mass to interior mass
            % mach:   Upstream mach number
            % igamma: 1: gamma = 5/3, 2: gamma = 4/3

            switch igamma
                case 2
                % gamma = 4/3
                d(1)  = 0.551005528;
                d(2)  = 0.450225376;
                d(3)  = -0.674076156;
                d(4)  = 0.594634757;
                d(5)  = -14.9500309;
                d(6)  = 0.315907987;
                d(7)  = -0.0203304265;
                d(8)  = -1.70433247;
                d(9)  = 30.4494171;
                d(10) = -0.0309267276;

                logCd = d(1) + d(2)*q + d(3)*mach + d(4)*q.*mach + d(5)*q.^2 + d(6)*mach.^2 +...
                        d(7)*q.*mach.^2 + d(8)*q.^2.*mach + d(9)*q.^3 + d(10)*mach.^3;

                case 1
                % gamma = 5/3
                d(1) = -0.15515258;
                d(2) = -3.03230504;
                d(3) = 0.27564942;
                d(4) = 0.19765314;
                d(5) = 1.41864156;
                d(6) = -0.0092128;

                logCd = d(1) + d(2)*q + d(3)*mach + d(4)*q.*mach + d(5)*q.^2 + d(6)*mach.^2;

                otherwise
                    warning('Unknown gamma. Only gamma = 4/3 or gamma = 5/3 is allowed.')
            end

            Cd = 10.^logCd;
        end
        
        %------------------------------------------------------------------
        % Calculate drag coefficient according to the fits by De+2020 to
        % their common-envelope wind tunnel simulations, given density
        % scale height and Mach number
        %------------------------------------------------------------------
        function Cd = getDragFromEpsandMach(mach,epsilon,igamma)
            % Epsilon: Ratio of density scale height, - d ln(rho) / dr, to
            %          the gravitational focusing radius
            % mach:    Upstream mach number
            % igamma:  1: gamma = 5/3, 2: gamma = 4/3
            
            % Convert epsilon to a mass ratio using the Mach number, using
            % Eq (10) of De+20
            q = CEwindTunnel.getQFromEpsilonAndMach(epsilon,mach);
            
            Cd = CEwindTunnel.getDragFromQandMach(q,mach,igamma);
            
        end
        
        %------------------------------------------------------------------
        % Calculate the local mass ratio (q_r) from the density scale
        % height parameter (epsilon_rho) and the upstream Mach number using
        % Eq. 10 of De+20
        %------------------------------------------------------------------
        function q = getQFromEpsilonAndMach(epsilon,mach)
            % We assume f_k = 1 and gamma/Gamma_s = 1
            
            blob = mach^2 / epsilon;
            q = blob - 1 - sqrt(blob * (blob - 2));
            if (~(q>0) || ~isreal(q))
                error('CEwindTunnel: q < 0')
            end
            
        end
    end
end

