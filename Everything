% ========== Graph templates ========== %
fontSize = 45;
set(groot, defaultTextInterpreter='latex', ...
    defaultAxesTickLabelInterpreter='latex', ...
    defaultLegendInterpreter='latex', ...
    defaultAxesFontSize=fontSize, ...
    defaultTextFontSize=fontSize, ...
    defaultTextFontSizeMode='manual', ...
    defaultLegendFontSize=fontSize, ...
    defaultLegendFontSizeMode='manual')

% ========== Cooling with forward difference ========== %
function [phi_before_wave_dt, phi_after_wave_dt] = coolingFD(General, Heat, Wave)

    % === Simulation setup === %
    dx = General.dx;
    dt = Heat.dt;
    L = General.L;

    tpoints = 0:dt:Heat.runtime+Wave.dt;
    tsize = length(tpoints);
    xpoints = -L:dx:L;
    xsize = length(xpoints);
    
    initial = NaN(xsize+2, 1);
    initial(2:end-1) = 2 .* General.ranbound .* rand(size(xpoints)) - General.ranbound;
    initial(1) = initial(2); initial(end) = initial(end-1);

    phi = NaN(xsize+2, 2); phi(:, 1) = initial;

    % fig = figure(Name = 'Initial plot', NumberTitle = 'off', ...
    %     Units = 'pixels', Position = [100, 100, General.figwidth, General.figheight]);
    % plot(xpoints, phi(2:end-1, 1), LineWidth=General.linewidth, Color=General.linecolor)
    % hold on;
    % plot(xpoints,     ones(size(xpoints)), '--', LineWidth=2, Color='#FFAD2A')
    % plot(xpoints, -1.*ones(size(xpoints)), '--', LineWidth=2,Color='#FFAD2A')
    % axis xy; ax = gca; ax.TickLength = [0 0];
    % xlim([-L, L]); ylim([-1.5 1.5]); xlabel('$x$'); ylabel('$\phi$')
    % title(['Heat\quad\quad\quad Seed ' num2str(General.seed) '\quad\quad\quad$r=' num2str(General.ranbound) '$\quad\quad\quad$t=$' num2str(Heat.runtime)])
    % fig.Color = General.backcolor; set(gca, color = General.backcolor);
    
    % === Video setup === %
    tstep = Heat.frametime / dt;
    fig = figure(Name = 'Heat Movie', NumberTitle = 'off', ...
        Units = 'pixels', Position = [100, 100, General.figwidth, General.figheight]);
    axis xy; ax = gca; ax.NextPlot = 'replaceChildren';
    xlim([-L, L]); ylim([-1.5 1.5]); ax.TickLength = [0 0];
    xlabel('$x$'); ylabel('$\phi$')
    v = VideoWriter([General.folder '/' Heat.moviename], 'MPEG-4');
    open(v)
    title(['Heat\quad\quad\quad Seed ' num2str(General.seed) '\quad\quad\quad$r=' num2str(General.ranbound) '$\quad\quad\quad$t=0$'])
    plot(xpoints, phi(2:end-1, 1), LineWidth=General.linewidth, Color=General.linecolor)
    fig.Color = General.backcolor; set(gca, color = General.backcolor);
    frame = getframe(fig);
    writeVideo(v, frame)
    
    % === Simulation and video === %
    for t = 1:tsize-1
        phi(2:end-1, 2) = dt.*((phi(3:end, 1) - 2.*phi(2:end-1, 1) + phi(1:end-2, 1))/dx^2 + 2.*phi(2:end-1, 1).*(1 - phi(2:end-1, 1).^2)) + phi(2:end-1, 1);
        phi(  1, 2) = phi(  1, 1) + (phi(    3, 2) - phi(    2, 2))*(dt/dx);
        phi(end, 2) = phi(end, 1) - (phi(end-1, 2) - phi(end-2, 2))*(dt/dx);
        phi(:, 1) = phi(:, 2);
        if rem(t, tstep) == 0
            title(['Heat\quad\quad\quad Seed ' num2str(General.seed) '\quad\quad\quad$r=' num2str(General.ranbound) '$\quad\quad\quad$t=$' num2str(tpoints(t+1), '%.1f')])
            % title(['$t=$' num2str(round(tpoints(t+1)))])
            plot(xpoints, phi(2:end-1, 1), LineWidth=General.linewidth, Color=General.linecolor)
            frame = getframe(fig);
            writeVideo(v, frame)
        end
        if tpoints(t+1) == Heat.runtime; phi_before_wave_dt = phi(:, 1); end
    end
    close(v)
    phi_after_wave_dt = phi(:, 1);

end

% ========== Cooling with centred difference ========== % NOT GOOD
function [phi, initial] = coolingCD(endtime, dt, frametime, dx, L, ranBound, seed)

    % === Simulation setup === %
    tpoints = 0:dt:endtime;
    tsize = length(tpoints);
    xpoints = -L:dx:L;
    xsize = length(xpoints);
    
    initial = ranBound .* rand(size(xpoints)) - ranBound/2;
    phi = NaN(xsize+2, 3); phi(2:end-1, 1:2) = [initial' initial'];
    phi(1, 1:2) = phi(2, 1:2); phi(end, 1:2) = phi(end-1, 1:2);
    
    % === Video setup === %
    tstep = frametime / dt;
    fig = figure(Name = 'Kink Dynamics Movie', NumberTitle = 'off');
    axis xy; ax = gca; ax.NextPlot = 'replaceChildren';
    xlim([-L, L]); ylim([-1.5 1.5]);
    v = VideoWriter(['Movies/cooling CD ' num2str(seed)], 'MPEG-4');
    open(v)
    title(['$t=0$'])
    plot(xpoints, phi(2:end-1, 1), LineWidth=1)
    frame = getframe(fig);
    writeVideo(v, frame)
    
    % === Simulation and video === %
    for t = 1:tsize-1
        phi(2:end-1, 3) = dt.*((phi(3:end, 2) - 2.*phi(2:end-1, 2) + phi(1:end-2, 2))/dx^2 + 2.*phi(2:end-1, 2).*(1 - phi(2:end-1, 2).^2)) + 2.*phi(2:end-1, 2) - phi(2:end-1, 1);
        phi(  1, 3) = phi(  1, 2) + (phi(    3, 3) - phi(    2, 3))*(dt/dx);
        phi(end, 3) = phi(end, 2) - (phi(end-1, 3) - phi(end-2, 3))*(dt/dx);
        phi(:, 1) = phi(:, 2); phi(:, 2) = phi(:, 3);
        if rem(t, tstep) == 0
            title(['$t=$' num2str(tpoints(t+1), '%.2f')])
            plot(xpoints, phi(2:end-1, 1), LineWidth=1)
            frame = getframe(fig);
            writeVideo(v, frame)
        end
    end
    close(v)
end

% ========== Expansion with centred difference ========== %
function phi = expansionCD(General, Heat, Wave, initial)

    % === Simulation setup === %
    dx = General.dx;
    dt = Wave.dt;
    L = General.L;

    tpoints = 0:dt:Wave.runtime;
    tsize = length(tpoints);
    xpoints = -L:dx:L;
    xsize = length(xpoints);

    phi = NaN(xsize+2, 3); phi(:, 1:2) = [initial initial];
    
    % === Video setup === %
    tstep = Wave.frametime / dt;
    fig = figure(Name = 'Kink Dynamics Movie', NumberTitle = 'off');
    axis xy; ax = gca; ax.NextPlot = 'replaceChildren';
    xlim([-L, L]); ylim([-1.5 1.5]);
    v = VideoWriter(['Movies/seed ' num2str(General.seed) ' expansion CD'], 'MPEG-4');
    open(v)
    title(['Seed ' num2str(General.seed) '\quad\quad\quad$t=$' num2str(Heat.runtime)])
    plot(xpoints, phi(2:end-1, 1), LineWidth=1)
    frame = getframe(fig);
    writeVideo(v, frame)
    
    % === Simulation and video === %
    for t = 1:tsize-1
        phi(2:end-1, 3) = dt^2.*((phi(3:end, 2) - 2.*phi(2:end-1, 2) + phi(1:end-2, 2))/dx^2 + 2.*phi(2:end-1, 2).*(1 - phi(2:end-1, 2).^2)) + 2.*phi(2:end-1, 2) - phi(2:end-1, 1);
        phi(  1, 3) = phi(  1, 2) + (phi(    3, 3) - phi(    2, 3))*(dt/dx);
        phi(end, 3) = phi(end, 2) - (phi(end-1, 3) - phi(end-2, 3))*(dt/dx);
        phi(:, 1) = phi(:, 2); phi(:, 2) = phi(:, 3);
        if rem(t, tstep) == 0
            % title(['$t=$' num2str(tpoints(t+1)+Heat.runtime, '%.2f')])
            title(['Seed ' num2str(General.seed) '\quad\quad\quad$t=$' num2str(round(tpoints(t+1)+Heat.runtime))])
            plot(xpoints, phi(2:end-1, 1), LineWidth=1)
            frame = getframe(fig);
            writeVideo(v, frame)
        end
    end
    close(v)

end

% ========== Expansion with centred difference with double initial ========== %
 function phi = expansion2IC(General, Heat, Wave, initial1, initial2)

    % === Simulation setup === %
    dx = General.dx;
    dt = Wave.dt;
    L = General.L;

    tpoints = 0:dt:Wave.runtime;
    tsize = length(tpoints);
    xpoints = -L:dx:L;
    xsize = length(xpoints);

    phi = NaN(xsize+2, 3);
    phi(:, 1) = initial1; phi(:, 2) = initial2;

    plotPhi = NaN(101); plotPhi(:, 1) = initial1(2:20:end-1);
    
    % === Video setup === %
    tstep = Wave.frametime / dt;
    fig = figure(Name = 'Wave Movie', NumberTitle = 'off', ...
        Units = 'pixels', Position = [100, 100, General.figwidth, General.figheight]);
    axis xy; ax = gca; ax.NextPlot = 'replaceChildren';
    xlim([-L, L]); ylim([-1.5 1.5]); ax.TickLength = [0 0];
    xlabel('$x$'); ylabel('$\phi$')
    v = VideoWriter([ General.folder '/' Wave.moviename], 'MPEG-4');
    open(v)
    title(['Wave\quad\quad\quad Seed ' num2str(General.seed) '\quad\quad\quad$r=' num2str(General.ranbound) '$\quad\quad\quad$t=$' num2str(Heat.runtime)])
    plot(xpoints, phi(2:end-1, 1), LineWidth=General.linewidth, Color=General.linecolor)
    fig.Color = General.backcolor; set(gca, color = General.backcolor);
    frame = getframe(fig);
    writeVideo(v, frame)

    % check50 = true;
    % check100 = true;
    
    % === Simulation and video === %
    for t = 2:tsize-1
        phi(2:end-1, 3) = dt^2.*((phi(3:end, 2) - 2.*phi(2:end-1, 2) + phi(1:end-2, 2))/dx^2 + 2.*phi(2:end-1, 2).*(1 - phi(2:end-1, 2).^2)) + 2.*phi(2:end-1, 2) - phi(2:end-1, 1);
        phi(  1, 3) = phi(  1, 2) + (phi(    3, 3) - phi(    2, 3))*(dt/dx);
        phi(end, 3) = phi(end, 2) - (phi(end-1, 3) - phi(end-2, 3))*(dt/dx);
        phi(:, 1) = phi(:, 2); phi(:, 2) = phi(:, 3);
        if rem(t, tstep) == 0
            % title(['$t=$' num2str(tpoints(t+1)+Heat.runtime, '%.2f')])
            title(['Wave\quad\quad\quad Seed ' num2str(General.seed) '\quad\quad\quad$r=' num2str(General.ranbound) '$\quad\quad\quad$t=$' num2str(round(tpoints(t+1)+Heat.runtime))])
            plot(xpoints, phi(2:end-1, 1), LineWidth=General.linewidth, Color=General.linecolor)
            frame = getframe(fig);
            writeVideo(v, frame)
            % if check50 && tpoints(t+1)+Heat.runtime == 50
            %     check50 = false;
            %     phi50 = phi(:, 1);
            % end
            % if check100 && tpoints(t+1)+Heat.runtime == 100
            %     check100 = false;
            %     phi100 = phi(:, 1);
            % end
        end
        if ismember(tpoints(t+1), tpoints(1:400:end))
            plotPhi(:, find(tpoints(1:400:end)==tpoints(t+1))) = phi(2:20:end-1, 1);
        end
    end
    close(v)

    if Wave.colormap
        fig = figure(Name = 'Colormap', NumberTitle = 'off', ...
            Units = 'pixels', Position = [100, 100, General.figwidth, General.figheight]);
        % fig.Color = General.backcolor; set(gca, color = General.backcolor);
        pcolor(xpoints(1:20:end), tpoints(1:400:end)+Heat.runtime, plotPhi');
        xlabel('$x$'); ylabel('time'); shading interp
        title(['Field $\phi$ through time for seed ' num2str(General.seed) ' and $r=$' num2str(General.ranbound)])
        colormap hsv; colorbar(gca, TickLabelInterpreter="latex"); clim([-1.5 1.5]); %ylim([10 70])
        exportgraphics(gca,[General.folder '/' Wave.moviename ' colormap' '.eps'])
    end
    
    if Wave.plot3D
        fig = figure(Name = '3D graph', NumberTitle = 'off', ...
        Units = 'pixels', Position = [100, 100, General.figwidth, General.figheight]);
        tiledlayout(Padding="loose"); nexttile
        % fig.Color = General.backcolor; set(gca, color = General.backcolor);
        surf(xpoints(1:20:end), tpoints(1:400:end)+Heat.runtime, plotPhi', ...
            FaceAlpha='0.5', FaceColor='interp', EdgeAlpha=0.1, EdgeColor='black');
        xlabel('$x$'); ylabel('time'); zlabel('$\phi$'); ylim([0 Wave.runtime]+Heat.runtime) %ylim([10 70])
        title(['Field $\phi$ through time for seed ' num2str(General.seed) ' and $r=$' num2str(General.ranbound)])
        colormap hsv; c = colorbar(gca, TickLabelInterpreter="latex"); clim([-1.5 1.5]);
        view(25, 20); c.Layout.Tile = 'east';
        exportgraphics(gca,[General.folder '/' Wave.moviename ' 3D plot' '.eps'])
    end

    % figure(Name = 'Noise', NumberTitle = 'off'); tlo = tiledlayout(1, 2, TileSpacing="tight");
    % 
    % tile1 = nexttile;
    % title('$t=50$'); plot(xpoints, phi50(2:end-1), LineWidth=1)
    % ylim([-1.5 1.5]); xlabel('$x$'); ylabel('$\phi$');
    % 
    % tile2 = nexttile;
    % title('$t=100$'); plot(xpoints, phi100(2:end-1), LineWidth=1)
    % ylim([-1.5 1.5]); xlabel('$x$'); yticklabels({})
    % 
    % title(tlo, 'Demonstration of energy buildup at boundaries over time', interpreter='Latex', fontsize=General.fontsize)

end

% ========== Expansion with cooling continuing with centred difference ========== % NOT GOOD
function phi = expandcool(endtime, dt, frametime, dx, L, initial, seed)

    % === Simulation setup === %
    tpoints = 0:dt:endtime;
    tsize = length(tpoints);
    xpoints = -L:dx:L;
    xsize = length(xpoints);
    divfact = 1/dt^2 + 5/dt;

    phi = NaN(xsize+2, 3); phi(:, 1:2) = initial;
    
    % === Video setup === %
    tstep = frametime / dt;
    fig = figure(Name = 'Kink Dynamics Movie', NumberTitle = 'off');
    axis xy; ax = gca; ax.NextPlot = 'replaceChildren';
    xlim([-L, L]); ylim([-1.5 1.5]);
    v = VideoWriter(['Movies/expandcool CD' num2str(seed)], 'MPEG-4');
    open(v)
    title('$t=10$')
    plot(xpoints, phi(2:end-1, 1), LineWidth=1)
    frame = getframe(fig);
    writeVideo(v, frame)
    
    % === Simulation and video === %
    for t = 1:tsize-1
        phi(2:end-1, 3) = ((phi(3:end, 2) - 2.*phi(2:end-1, 2) + phi(1:end-2, 2))/dx^2 + 2.*phi(2:end-1, 2).*(1 - phi(2:end-1, 2).^2)) / divfact + 2.*phi(2:end-1, 2) - phi(2:end-1, 1);
        phi(  1, 3) = phi(  1, 2) + (phi(    3, 3) - phi(    2, 3))*(dt/dx);
        phi(end, 3) = phi(end, 2) - (phi(end-1, 3) - phi(end-2, 3))*(dt/dx);
        phi(:, 1) = phi(:, 2); phi(:, 2) = phi(:, 3);
        if rem(t, tstep) == 0
            title(['$t=$' num2str(round(tpoints(t+1)))])
            plot(xpoints, phi(2:end-1, 1), LineWidth=1)
            frame = getframe(fig);
            writeVideo(v, frame)
        end
    end
    close(v)
end

% ========== Expansion with cooling initial with forward difference ========== % NOT GOOD
function phi = excool(endtime, dt, frametime, dx, L, initial, seed)

    % === Simulation setup === %
    tpoints = 0:dt:endtime;
    tsize = length(tpoints);
    xpoints = -L:dx:L;
    xsize = length(xpoints);
    divfact = 1/dt^2 + 1/dt;

    phi = NaN(xsize+2, 2); phi(2:end-1, 1) = initial;
    phi(1, 1) = phi(2, 1); phi(end, 1) = phi(end-1, 1);
    
    % === Video setup === %
    tstep = frametime / dt;
    fig = figure(Name = 'Kink Dynamics Movie', NumberTitle = 'off');
    axis xy; ax = gca; ax.NextPlot = 'replaceChildren';
    xlim([-L, L]); ylim([-1.5 1.5]);
    v = VideoWriter(['Movies/excool FD' num2str(seed)], 'MPEG-4');
    open(v)
    title('$t=0$')
    plot(xpoints, phi(2:end-1, 1), LineWidth=1)
    frame = getframe(fig);
    writeVideo(v, frame)
    
    % === Simulation and video === %
    for t = 1:tsize-1
        phi(2:end-1, 2) = ((phi(3:end, 1) - 2.*phi(2:end-1, 1) + phi(1:end-2, 1))/dx^2 + 2.*phi(2:end-1, 1).*(1 - phi(2:end-1, 1).^2)) / divfact + phi(2:end-1, 1);
        phi(  1, 2) = phi(  1, 1) + (phi(    3, 2) - phi(    2, 2))*(dt/dx);
        phi(end, 2) = phi(end, 1) - (phi(end-1, 2) - phi(end-2, 2))*(dt/dx);
        phi(:, 1) = phi(:, 2);
        if rem(t, tstep) == 0
            title(['$t=$' num2str(round(tpoints(t+1)))])
            plot(xpoints, phi(2:end-1, 1), LineWidth=1)
            frame = getframe(fig);
            writeVideo(v, frame)
        end
    end
    close(v)
end

% ========== Simulation settings ========== %
for seed = 6
r = 1;

% for r = 1 %[0.001 0.01 0.1 1 10 100]
% seed = 6;

rng(seed)

General.seed      = seed;
General.ranbound  = r;
General.L         = 10;
General.dx        = 1e-2;
General.fontsize  = fontSize;
General.folder    = 'Misc movies';
General.backcolor = '#181818';
General.linecolor = '#6DCBC9';
General.linewidth = 1;
General.figwidth  = 1280;
General.figheight = 720;

Heat.dt           = 1e-5;
Heat.runtime      = 10;
Heat.frametime    = 0.05;
Heat.moviename    = ['D r=' num2str(General.ranbound) ' cool seed 6']; 

Wave.dt           = 1e-3;
Wave.runtime      = 90;
Wave.frametime    = 0.2;
Wave.moviename    = ['D r=' num2str(General.ranbound) ' expa seed 6'];
Wave.colormap     = true;
Wave.plot3D       = true;

% ========== Simulations ========== %
[initial1, initial2] = coolingFD(General, Heat, Wave);
                    expansion2IC(General, Heat, Wave, initial1, initial2);

end
