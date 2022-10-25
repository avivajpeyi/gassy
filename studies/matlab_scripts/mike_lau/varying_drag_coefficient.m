% addpath('~/Documents/MATLAB/utils');
const = constants; % Create instance of constants class

%%
% Create instance of stellar profile class, importing the MESA profile
filepath = './cored_adiabatic_fixedcomposition.dat';
starDescription = 'Softened profile, $\gamma = 5/3$';
star = stellarProfile(filepath, starDescription);

% Add back core mass for adiabatic profile
star.mass = star.mass + 3.15 * const.MSUN;
star.stellarMass = star.stellarMass + 3.15 * const.MSUN;
star.cs = sqrt(5/3 * star.pres ./ star.dens);

% Create instance of hydro simulation sink data
sink1path = './binarySink0001N01.ev'; % Data of primary core
sink2path = './binarySink0002N01.ev'; % Data of secondary
sink1 = ptmass(sink1path);
sink2 = ptmass(sink2path);
sink2.centreOnPrimaryCore(sink1); % Centre companion pos, vel on primary core


%% Create instances of fixed inspiral
xFactor = 0.95; % Initial radius of companion as fraction of donor radius
rotFactor = 0.4; % Rotation of donor as a function of Keplerian angular velocity at surface (this decreases drag)

% Fixed drag coefficients
nCdbins = 5;
CdArray = linspace(0.9,2.1,nCdbins); % Multiplicative factor to Bondi-Hoyle-Lyttleton drag

% Step settings
period = 2*pi * sqrt( (xFactor * star.stellarRadius)^3 / (const.G * star.stellarMass) );
dt = 0.001 * period;
endsep = 15 * const.RSUN;

% Create instances of CE inspiral
for i = 1:nCdbins
   CEcell{i} = analyticalInspiral(dt,endsep,xFactor,rotFactor,star,sink2.mass(1) * const.MSUN);
   CEcell{i}.setDragCoefficient(0,CdArray(i));
end

% Create instances of inspiral with De+20 drag coefficient
CEcell{nCdbins+1} = analyticalInspiral(dt,endsep,xFactor,rotFactor,star,sink2.mass(1) * const.MSUN); 
CEcell{nCdbins+1}.setDragCoefficient(1); % Gamma = 5/3
CEcell{nCdbins+1}.USE_BHLRAD = false;
CEcell{nCdbins+2} = analyticalInspiral(dt,endsep,xFactor,rotFactor,star,sink2.mass(1) * const.MSUN); 
CEcell{nCdbins+2}.setDragCoefficient(2); % Gamma = 4/3
CEcell{nCdbins+2}.USE_BHLRAD = false;


%% Match initial conditions for inspiral to full hydro simulation
sinkradius = sqrt( sink2.x.^2 + sink2.y.^2 + sink2.z.^2 ) * const.RSUN;
startidx = find( sinkradius < xFactor * star.stellarRadius, 1);
vx0 = sink2.vx(startidx);
vy0 = sink2.vy(startidx);
x0 = sink2.x(startidx) * const.RSUN;
y0 = sink2.y(startidx) * const.RSUN;

for i = 1:nCdbins+2
   CEcell{i}.vel(1,1:2) = [vx0, vy0];
   CEcell{i}.pos(1,1:2) = [x0, y0];
end

%% Integrate EOM using semi-implicit Euler method
for i = 1:nCdbins+2
    CEcell{i}.integrateEoM;
end

%% Plotting

% Colours
mycmap = hsv(nCdbins+2);

% Sink plot settings
sink_timeshift_yrs = 10.57; 
sinkplotrange = 97000:size(sink2.x(:))/2;
sinktplot = sink2.time(sinkplotrange) - sink_timeshift_yrs;



%%
%-------------------------------------------------------------------------
% Plot separation, drag
%--------------------------------------------------------------------------
sinksep = sqrt(sink2.x.^2 + sink2.y.^2);
clearvars leftPlotArray rightPlotArray

yyaxis left
hold on
for i = 1:nCdbins
   legendText = sprintf('Sep, $C_d = %.2f$',CEcell{i}.dragCoff(1));
   leftPlotArray(i) = CEcell{i}.plt('sep',mycmap(i,:),'-',legendText);
end
leftPlotArray(nCdbins+1) = CEcell{nCdbins+1}.plt('sep',mycmap(nCdbins+1,:),'-','Sep, interp. $C_d$, $\gamma = 5/3$');
leftPlotArray(nCdbins+2) = CEcell{nCdbins+2}.plt('sep',mycmap(nCdbins+2,:),'-','Sep, interp. $C_d$, $\gamma = 4/3$');
leftPlotArray(nCdbins+3) = plot(sinktplot, sinksep(sinkplotrange),'-k', 'DisplayName', 'Sep, from hydro');
hold off
yline(star.stellarRadius / constants.RSUN,':m');
ylabel('Separation / $R_\odot$','interpreter','latex');


yyaxis right
hold on
for i = 1:nCdbins
   legendText = sprintf('Drag, $C_d = %.2f$',CEcell{i}.dragCoff(1));
   rightPlotArray(i) = CEcell{i}.plt('drag',mycmap(i,:),'--',legendText); 
end
rightPlotArray(nCdbins+1) = CEcell{nCdbins+1}.plt('drag',mycmap(nCdbins+1,:),'--','Drag, interp. $C_d$, $\gamma = 5/3$');
rightPlotArray(nCdbins+2) = CEcell{nCdbins+2}.plt('drag',mycmap(nCdbins+2,:),'--','Drag, interp. $C_d$, $\gamma = 4/3$');
ylabel('$F_\textrm{drag} / M_2$ / cms${}^{-2}$','Interpreter','latex')

xlabel('$t$ / yr','Interpreter','latex')
legend([leftPlotArray,rightPlotArray],'Interpreter','latex')

%%
%--------------------------------------------------------------------------
% Plot separation, Mach number
%--------------------------------------------------------------------------
clearvars leftPlotArray rightPlotArray
yyaxis left
hold on
for i = 1:nCdbins
   legendText = sprintf('Sep, $C_d = %.2f$',CEcell{i}.dragCoff(1));
   leftPlotArray(i) = CEcell{i}.plt('sep',mycmap(i,:),'-',legendText);
end
leftPlotArray(nCdbins+1) = CEcell{nCdbins+1}.plt('sep',mycmap(nCdbins+1,:),'-','Sep, interp. $C_d$, $\gamma = 5/3$');
leftPlotArray(nCdbins+2) = CEcell{nCdbins+2}.plt('sep',mycmap(nCdbins+2,:),'-','Sep, interp. $C_d$, $\gamma = 4/3$');
leftPlotArray(nCdbins+3) = plot(sinktplot, sinksep(sinkplotrange),'-k', 'DisplayName', 'Sep, from hydro');
hold off
yline(star.stellarRadius,':m');
ylabel('Separation / cm')


yyaxis right
hold on
for i = 1:nCdbins
   legendText = sprintf('$\\mathcal{M}$, $C_d = %.2f$',CEcell{i}.dragCoff(1));
   rightPlotArray(i) = CEcell{i}.plt('mach',mycmap(i,:),'--',legendText); 
end
rightPlotArray(nCdbins+1) = CEcell{nCdbins+1}.plt('mach',mycmap(nCdbins+1,:),'--','$$\mathcal{M}$$, interp. $C_d$, $\gamma = 5/3$');
rightPlotArray(nCdbins+2) = CEcell{nCdbins+2}.plt('mach',mycmap(nCdbins+2,:),'--','$$\mathcal{M}$$, interp. $C_d$, $\gamma = 4/3$');
ylabel('$\mathcal{M}$','Interpreter','latex')
legend([leftPlotArray,rightPlotArray],'Interpreter','latex')


%%
%--------------------------------------------------------------------------
% Drag coefficient
%--------------------------------------------------------------------------
clearvars plotArray
hold on
for i = 1:nCdbins
   legendText = sprintf('$C_d = %.2f$',CEcell{i}.dragCoff(1));
   plotArray(i) = CEcell{i}.plt('Cd',mycmap(i,:),'-',legendText); 
end
plotArray(nCdbins+1) = CEcell{nCdbins+1}.plt('Cd',mycmap(nCdbins+1,:),'-','Interpolated $C_d$, $\gamma = 5/3$');
plotArray(nCdbins+2) = CEcell{nCdbins+2}.plt('Cd',mycmap(nCdbins+2,:),'-','Interpolated $C_d$, $\gamma = 4/3$');

ylabel('$C_d$','Interpreter','latex')
xlabel('$t$ / yr','Interpreter','latex')
legend([plotArray],'Interpreter','latex')

%%
%--------------------------------------------------------------------------
% Total energy
%--------------------------------------------------------------------------
hold on
for i = 1:nCdbins
   legendText = sprintf('$C_d = %.2f$',CEcell{i}.dragCoff(1));
   plotArray(i) = CEcell{i}.plt('ene',mycmap(i,:),'-',legendText); 
end
plotArray(nCdbins+1) = CEcell{nCdbins+1}.plt('ene',mycmap(nCdbins+1,:),'-','Interpolated $C_d$');
xlabel('$t$ / yr','Interpreter','latex')
ylabel('Total energy')
set(gca,'YScale','log')
legend([plotArray],'Interpreter','latex')

%%
%--------------------------------------------------------------------------
% Energy error
%--------------------------------------------------------------------------
% nexttile 
% hold on
% plot(tplot,(CE.totene(plotrange) - CE.totene(1)) / CE.totene(1));
% set(gca,'YScale','log')
% xlabel('$t$ / yr','Interpreter','latex')
% ylabel('$\Delta E / E_0$','Interpreter','latex')


   
%% Plot 3-d sink trajectory in globe-like star
h = figure;
filename = 'testAnimated.gif';

% Plot sink trajectories
curve2 = animatedline('LineWidth',0.5,'Color','b');
curve1 = animatedline('LineWidth',0.5,'Color','r');
sink2x = sink2.x;
sink2y = sink2.y;
sink2z = sink2.z;
sink1x = sink1.x;
sink1y = sink1.y;
sink1z = sink1.z;
for k = 1:8
    sink2x(1:2:end) = [];
    sink2y(1:2:end) = [];
    sink2z(1:2:end) = [];
    sink1x(1:2:end) = [];
    sink1y(1:2:end) = [];
    sink1z(1:2:end) = [];
end

view(43,24);
lim = 6e13;
set(gca,'XLim',[-lim lim],'YLim',[-lim lim],'ZLim',[-lim lim]);
axis square
for i = 1:length(sink2x)
    addpoints(curve2, sink2x(i), sink2y(i), sink2z(i));
    addpoints(curve1, sink1x(i), sink1y(i), sink1z(i));
    drawnow
    
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    
    % Write to the GIF File 
    if i == 1 
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.1); 
    else 
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1); 
    end 
end


%% Plot globe-like donor star

plot3(sink2.x, sink2.y, sink2.z, 'b'); hold on
plot3(sink1.x, sink1.y, sink1.z, 'r');
  
latspacing = 20; 
lonspacing = 20; 
plotGlobe(star.stellarRadius, latspacing, lonspacing); hold off

xlabel('x / cm')
ylabel('y / cm')
zlabel('z / cm')
axis equal


%% Plot animation of 2-d analytical inspiral
h = figure;
idx = 1;
filename = 'analyticalInspiralAnimated.gif';

% Plot sink trajectories
curve = animatedline('LineWidth',0.5,'Color','b');

%view(43,24);
lim = 1.1 * star.stellarRadius / constants.RSUN;
set(gca,'XLim',[-lim lim],'YLim',[-lim lim]);
axis square


hold on
plotcircle(0.,0.,star.stellarRadius / constants.RSUN);
xlabel('$x$ / $R_\odot$', 'Interpreter', 'latex');
ylabel('$y$ / $R_\odot$', 'Interpreter', 'latex');
writecount = 1;
for i = 1:length(CEcell{idx}.time)
    addpoints(curve, CEcell{idx}.pos(i,1) / constants.RSUN, CEcell{idx}.pos(i,2) / constants.RSUN);
    if ( (writecount > 25) || (i == 1) )
        delete(findall(gcf,'Tag','stream'));
        annstring = sprintf('t = %0.2f yr', CEcell{idx}.time(i) / 3.15e7);
        annotation('textbox',[0.65 0.65 0.20 0.25],'String',annstring,...
            'FitBoxToText','on','EdgeColor','none','Tag','stream');
        drawnow
        % Capture the plot as an image 
        frame = getframe(h); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 

        % Write to the GIF File 
        if i == 1 
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.00001); 
        else 
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.00001); 
        end 
        writecount = 0;
    else
        writecount = writecount + 1;
    end
end


%%
%--------------------------------------------------------------------------
% Plot trajectories
%--------------------------------------------------------------------------
plotrange = 1:1000;
clearvars plotTrajectories

%tiledlayout(1,3);
%nexttile
plotcircle(0.,0.,star.stellarRadius); hold on
plotTrajectories(1) = CEcell{6}.plttraj('b','-','Fixed $C_d$');
plotTrajectories(2) = CEcell{7}.plttraj('r','-','Fixed $C_d$'); hold off
xlabel('x / cm');
ylabel('y / cm');
legend(plotTrajectories,'Interpreter','latex')
axis equal tight

%nexttile
%plotcircle(0.,0.,star.stellarRadius); hold on
%plt_trajectory2 = plot(CEcell{6}.pos(plotrange,1), CEcell{6}.pos(plotrange,2), 'r', 'DisplayName', 'Interpolated $C_d$'); hold off
%xlabel('x / cm');
%ylabel('y / cm');
%legend(plt_trajectory2,'Interpreter','latex')
%axis equal tight

%nexttile
%plotcircle(0.,0.,star.stellarRadius); hold on
%plt_trajectoryhydro = plot(sink2.x, sink2.y, 'k', 'DisplayName', 'From hydro'); hold off
%xlabel('x / cm');
%ylabel('y / cm');
%legend(plt_trajectoryhydro,'Interpreter','latex')
%axis equal tight

%% 
%--------------------------------------------------------------------------
% Plot helical trajectory
%--------------------------------------------------------------------------
plt_helix = plot3(CEcell{1}.pos(:,1), CEcell{1}.pos(:,2),...
                        CEcell{1}.time(:) / constants.SEC_PER_YR, 'b');
xlabel('x / cm');
ylabel('y / cm');
zlabel('t / yr');



   

