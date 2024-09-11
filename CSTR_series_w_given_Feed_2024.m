% CSTR in Series with given feed rates. Dilutions varied.
 
clc
clear
close all

% Plot Configuration for labeling of vessel 1 figure
plot_cfg = struct( xlabel = "Dilution, D1 (1/h)", ...
    ylabel = ["Concentrations (g/L)", "Concentrations (g/L)", "Concentrations (g/L)", ...
    "Productivities (g/h)", "Productivities (g/h)"], ...
    legend = ["Substrate", "Cell", "Product", "Cell Productivity (DX)", "Productivity (PD)"], ...
    lineStyle  = ["-", "-", "-", "-", "-"], lineColor = ["r", "b", "g", "b", "g"], ...
    layout = [2 1], subplotIndex = [1, 1, 1, 2, 2], yAxis = ["L", "L", "L", "L", "L"], ...
    y_max = [inf, inf]);

% Cell and Product Model Constants
mumax  = 0.48;      % max growth rate constant; 1/h
YXS    = 0.20;      % gX/gS, g cell per g substrate
YPS    = 0.5;       % gP/gx, g product per g dry cell
alpha  = 3.3;       % productivity; gP per gX*h
beta   = 0.05;      % productivity; 1/h
KS     = 1.5;       % gS/L, g substrate per L
kd     = 0.01;      % death rate constant; 1/h

% Inlet Concentrations & Maximum Dilution
SF = 75.0;          % g/L, inlet substrate concentration
D_max = 0.99*mumax; % maximum dilution factor (g/l) -- do not go higher than mumax

% Model Equations for Cell and Product for Vessel #1
n_steps   = 100;                         % linspace observation
D1 = linspace(0.15, D_max, n_steps)';    % dilution rate column vector
[S1, X1, P1, DX1, PX1] = vessel_1(D1, kd, KS, mumax, SF, alpha, beta, YXS, YPS);

% Plot results for Vessel #1
y1 = [S1, X1, P1, DX1, PX1];
plot_results(D1, y1, plot_cfg)

% Additional parameters for 2nd vessel in series
F        = 20;       % Feed rate into first vessel; L/h
F_prime  = 60;       % Feed rate into second vessel; L/h

alpha2   = 3.3;      % vessel #2 productivity; gP per gX*h
beta2    = 0.05;     % vessel #2 productivity; 1/h
SF_prime = 75;       % g/L, inlet substrate concentration for the second vessel
D2plot_mult = 1.5;   % range of D2 relative to D1
D2 = linspace(0.15, D2plot_mult * D_max, n_steps+20)';  % Dilution rate for Vessel #2

% Matrix of dilutions for which calculations will be run
[m.D2, m.D1] = meshgrid(D2, D1);

% Volume calculations for V1 and V2 based on the flow rates
m.V1 = F ./ m.D1;
m.V2 = (F + F_prime) ./ m.D2;

% Set any V1 and V2 == Inf to NaN 
if any(m.V1 == Inf | m.V2 == Inf, "all")
    inv_idx = (m.V1 == Inf | m.V2 == Inf);
    [m.V1(inv_idx), m.V2(inv_idx)] = deal(NaN);
end

% initialize matrices for vessel 2 iterative calculations
[m.mu2, m.S2, m.X2, m.P2, m.count] = deal(NaN(size(m.D1)));

fprintf('Computation Progress: %3d %3d\n',0,0)

max_count = 0;      % variable to track maximum number of iterations for convergence
max_iter  = 250;    % maximum iterations allowed
toler     = 1e-4;   % convergence tolerance (1e-4 == 0.01%)

for k1 = 1:length(D1)
    for k2 = 1:length(D2)

        fprintf('\b\b\b\b\b\b\b\b %3d %3d', k1, k2)

        D1_k1 = D1(k1);   % also stored as m.D1(k1,k2)
        P1_k1 = P1(k1);
        X1_k1 = X1(k1);
        S1_k1 = S1(k1);
        D2_k2 = D2(k2);  % also stored as m.D2(k1,k2)
        V2_k2 = m.V2(k1,k2);

        % if any variable extraction is nan, then go to next iteration
        if any(isnan([D1_k1, P1_k1, X1_k1, S1_k1, D2_k2, V2_k2]))
            continue
        end

        % initialize variables before iterative solution loop
        X2 = X1_k1;   % Initial guess of X2
        delta = 1e9;  % Arbitrary big number greater than tolerance
        count = 0;    % Set iteration counter

        while delta > toler && count < max_iter
            mu2 = D2_k2 - ((F * X1_k1) / (V2_k2 * X2)) + kd;
            S2 = KS * mu2 / (mumax - mu2);
            qp2 = alpha2 * mu2 + beta2;
            X2_num = (F / V2_k2 * S1_k1) + (F_prime / V2_k2 * SF_prime) - D2_k2 * S2;
            X2_new = X2_num / (mu2/YXS + qp2/YPS);
            delta = abs((X2_new - X2) / X2);
            X2 = X2_new;
            count = count + 1;
        end

        % skip results where mu2 is outside physical reality
        if mu2 <=0 || mu2 >= mumax
            continue
        end

        if count ~= max_iter && count > max_count
            max_count = count;
        end

        % store good solutions (i.e., convergence achieved)
        if count ~= max_iter
            P2 = (F * P1_k1) / (V2_k2*D2_k2) + (qp2 * X2 / D2_k2);
            [m.P2(k1,k2), m.mu2(k1,k2), m.S2(k1,k2), m.X2(k1,k2), ...
                m.count(k1,k2)] = deal(P2, mu2, S2, X2, count);
        end

    end  % ~ vessel #2 dilution domain
end      % ~ vessel #1 dilution domain

fprintf('\nMax Count: %3d\n', max_count)

% Remove invalid values (<0) for S2, X2, P2
inv_idx = (m.S2 < 0) | (m.X2 < 0) | (m.P2 < 0) ;
if any(inv_idx, "all")
    % warning('Invalid values found: %d indices set to NaN.', sum(inv_idx, 'all'));
    [m.X2(inv_idx), m.S2(inv_idx), m.P2(inv_idx)] = deal(NaN);
end

% Calculate Productivities
m.XD2 = m.X2 .* m.D2;
m.PD2 = m.P2 .* m.D2;

% Plot vessel 2 results
surf_conf = struct( ...
    'axis_1', "Dilution, D2 (g/L)", 'axis_2', "Dilution, D1 (g/L)", ...  
    'subplot_labels', ["mu2", "S2", "V1", "V2"], ...
    'colormap', jet, 'view', [0, 90] );

% first figure set (F_prime, mu2, S2, count)
plot_surf_wrap(m.D1, m.D2, m, surf_conf);

% Update surf_conf for the second set of plots
surf_conf.subplot_labels = ["X2", "P2", "XD2", "PD2"];

% second figure set (X2, P2, XD2, PD2)
plot_surf_wrap(m.D1, m.D2, m, surf_conf);

%% Vessel 1 Calculations
function [S, X, P, XD, PD] = vessel_1(D, kd, KS, mumax, SF, alpha, beta, YXS, YPS)

    % Model equations for Cell and Product
    mu   = D + kd;                     % growth rate
    S    = (mu * KS) ./ (mumax - mu);  % substrate concentration
    qp   = alpha * mu + beta;          % product formation rate
    
    if any(qp == 0, "all")
        warning('Product formation rate, qp, equal to 0!')
    end
    
    % Cell and Product Concentrations
    X_num = D .* (SF - S);             % cell concentration numerator
    X_den = ((mu / YXS) + (qp / YPS)); % cell concentration denominator
    X =  X_num ./ X_den;               % cell concentration
    P    = qp .* X ./ D;               % product concentration
    XD   = D .* X;                     % cell productivity
    PD   = D .* P;                     % product productivity
    
    % Remove invalid values for S, X, P, XD, and PD
    inv_idx = (X <= 0) | (S <= 0) | (P < 0);
    if any(inv_idx, "all")
        % warning('vessel_1:inv_idx', '%d invalid indices found. Values set to NaN.', sum(inv_idx, 'all'));
        [S(inv_idx), X(inv_idx), P(inv_idx), XD(inv_idx), PD(inv_idx)] = deal(NaN);
    end

end

%% Surface Plot wrapper function
function plot_surf_wrap(D1, D2, m_struct, surf_conf)
    % Inputs:
    % - D1, D2: Dilution rates for Vessel 1 and 2
    % - m_struct: Structure containing the matrices to plot (e.g., F_prime, mu2, etc.)
    % - surf_conf: Configuration structure for the surface plots
    
    % Initialize an array to store subplot handles
    h = gobjects(1, length(surf_conf.subplot_labels));
    
    % Loop over subplots (2x2 grid)
    figure(Color="w");
    for sp = 1:length(surf_conf.subplot_labels)
        subplot_label = strrep(surf_conf.subplot_labels{sp}, '_', '\_');
        h(sp) = subplot(2, 2, sp);
        plot_surf(D1, D2, m_struct.(surf_conf.subplot_labels{sp}), ...
            [0 inf], [0 inf], [-inf inf], ...
            surf_conf.axis_1, surf_conf.axis_2, ...
            subplot_label, surf_conf.colormap, surf_conf.view);
    end
    
    % Enable rotation and link the camera properties for all subplots
    for sp = 1:length(h)
        rotate3d(h(sp), 'on');
        % Set each subplot to manual camera mode
        set(h(sp), 'CameraPositionMode', 'manual', ...
            'CameraTargetMode', 'manual', ...
            'CameraUpVectorMode', 'manual', ...
            'CameraViewAngleMode', 'manual');
    end
    
    % Add listeners to all subplots to synchronize rotation
    for sp = 1:length(h)
        addlistener(h(sp), 'View', 'PostSet', @(src, evt) syncRotations(h, sp));
    end
end

%% Synchronize Rotations 
function syncRotations(h, sourceIndex)
    % Synchronize the rotations of all subplots to the one that was rotated
    currentView = get(h(sourceIndex), 'View');
    for sp = 1:length(h)
        if sp ~= sourceIndex
            set(h(sp), 'View', currentView);
        end
    end
end

%% Surface Plot Function
function plot_surf(D1, D2, Z, xlim_values, ylim_values, zlim_values, ...
    axis_1, axis_2, zlabel_str, colormap_option, view_cfg)

    % Screen Z for values outside limits
    idx = Z < zlim_values(1) | Z > zlim_values(2);
    Z(idx) = NaN;
    
    % Plot the surface
    surf(D2, D1, Z, 'EdgeColor', 'none');
    set(gca, 'XDir', 'normal', 'YDir', 'reverse');  % Adjust axis directions
    xlabel(axis_1);
    ylabel(axis_2);
    zlabel(zlabel_str);
    title(zlabel_str);  % Title matches the z-axis label
    axis([xlim_values ylim_values zlim_values]);
    view(view_cfg(1), view_cfg(2));  % Set view
    c = colorbar;
    colormap(colormap_option);
    c.Label.String = zlabel_str;
end

%% Plotting Function
function plot_results(x, y, cfg)

lwd = 2; % Line width
fsize = 14; % Font size for labels
n_subplots = length(unique(cfg.subplotIndex)); % Number of unique subplots

% Determine whether overwriting an existing figure
no_figure = isempty(findobj('Type', 'figure'));
if no_figure
    HV = "on";
else
    HV = "off";
end

for sp = 1:n_subplots

    currentPlotIndices = find(cfg.subplotIndex == sp); % current subplot index
    subplot(cfg.layout(1), cfg.layout(2), sp)
    hold on

    yyaxis left % Start with left axis by default
    ax = gca; % Get current axes
    ax.YAxis(2).Visible = 'off'; % Explicitly turn off right axis

    for k = currentPlotIndices
        if strcmpi(cfg.yAxis(k), "L")
            yyaxis left
            plot(x, y(:, k), Color=cfg.lineColor(k), LineStyle=cfg.lineStyle(k), ...
                LineWidth=lwd, Marker="none", HandleVisibility=HV)
            ylabel(cfg.ylabel(k), Color="k", FontSize=fsize)
            ax.YAxis(1).Color = "k";
        else
            yyaxis right
            plot(x, y(:, k), Color=cfg.lineColor(k), LineStyle=cfg.lineStyle(k), ...
                LineWidth=lwd, Marker="none", HandleVisibility=HV)
            ylabel(cfg.ylabel(k), Color=cfg.lineColor(k), FontSize=fsize)
            ax.YAxis(2).Color = cfg.lineColor(k);
            ax.YAxis(2).Visible = 'on';
        end
    end

    axis([0 inf 0 cfg.y_max(sp)])
    xlabel(cfg.xlabel, FontSize=fsize)
    grid on

    % Add legend if no figure exists (i.e., not overwriting)
    if  no_figure
        legend(cfg.legend(currentPlotIndices), Location="best", TextColor="k")
    end

end

% Set figure to white background
set(gcf, Color="w");

% Link axes
linkaxes(findall(gcf,'type','axes'), 'x');

end
