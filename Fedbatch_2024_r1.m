%% Fed-Batch Monod Growth Model with Substrate Inhibition
% Simulates and plots a batch followed by fed-batch vessel with Monod growth
% without product formation and substrate inhibition in the growth model.

clc
clear
close all   % Comment to overlay results on existing figure

%% Main Program

% Plot Configuration for Labeling
plot_cfg = struct( xlabel = "Time (h)", ...
    ylabel = ["Substrate (g/L)", "Cell Density (g/L)", "Growth Rate (1/h)", "Volume (L)"], ...
    legend = ["Substrate", "Cell Density", "Growth Rate", "Volume"], ...
    lineStyle  = ["-", "-", "-", "-"], lineColor = ["r", "b", "m", "k"], ...
    layout = [4 1], subplotIndex = [1, 2, 3, 4], yAxis = ["L", "L", "L", "L"], ...
    y_max = [inf, inf, inf, inf]);  %y_max = [150, 150, 1, 200]);  % 
    
% Experiment Length & Initial Conditions
tmax  = 25;           % length of experiment/simulation run time (hours)
tspan = 0:0.01:tmax;  % time steps (h), 0 to tmax with step size 0.01
S0    = 4.0;          % initial conc of substrate (g/L)
X0    = 0.0052;       % initial cell density (g/L)
V0    = 100;          % initial volume (L)

% Model Constants
mumax = 0.95;  % max growth rate constant (1/h)
KS    = 0.65;  % half-saturation constant, g substrate per L (gS/L)
KI    = 0.02;  % substrate inhibition constant, (L/gS)^2
YXS   = 0.48;  % g cell per g substrate (gX/gS)

% Fed-Batch Constants
muset = 0.32;  % exponential growth rate setpoint
SF    = 481;   % feed concentration (g/L)
TB    = 7;     % starting feed time (hours)
Vmax  = 200;   % maximum volume (L)

Fi    = (muset * S0 * V0 / SF);   % Exponential feed rate term 

params = struct("mumax", mumax, "KS", KS, "KI", KI, "YXS", YXS, ...
                "muset", muset, "SF", SF, "TB", TB, "Vmax", Vmax, ...
                "Fi", Fi);

% ODE solver (passes initial conditions to function)
[t, y] = ode23(@(t, y) fedbatch_model(t, y, params), tspan, [X0; S0; V0]);

% Extract results to state variable names
[X, S, V] = deal(y(:, 1), y(:, 2), y(:, 3));

% Calculate Growth Rate
mu = mumax * S ./ (KS + S + KI .* S .* S);

% Plot results
plot_results(t, [S, X, mu, V], plot_cfg)

%% Fed-Batch Model with Substrate Inhibition
function dydt = fedbatch_model(t, y, params)

    % Unpack the state variables X, S, and V from y
    [X, S, V] = deal(y(1), y(2), y(3));
    
    % Unpack the parameters into individual variables
    [mumax, KS, KI, YXS, muset, SF, TB, Vmax, Fi] = deal(params.mumax, ...
            params.KS, params.KI, params.YXS, params.muset, params.SF, ...
            params.TB, params.Vmax, params.Fi);
    
    % Calculate feed rate F based on time t
    if t < TB
        F = 0.0;
    elseif V < Vmax
        % F = 0;         % BATCH  
        % F = 6.6667;  % 100/15 =6.6667; constant feed rate
        F = Fi * exp(muset * (t - TB));
    else
        F = 0.0;
    end

    % Growth rate with substrate inhibition
    mu = mumax * S / (KS + S + (S * S * KI));
    
    % Design Equations
    dXdt = mu * X - F * X / V;
    dSdt = F * (SF - S) / V - mu * X / YXS;
    dVdt = F;

    dydt = [dXdt; dSdt; dVdt];

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
