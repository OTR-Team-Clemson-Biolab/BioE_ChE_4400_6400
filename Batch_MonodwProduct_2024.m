% Monod Model with Product
% Simulates and plots a batch vessel with monod growth and
% mixed growth and non-associated product formation

clc
clear
close all   % comment to overlay results on existing figure

%% Main Program

% Plot Configuration for Labeling
plot_cfg = struct( xlabel = "Time (h)", ...
    ylabel = ["Cell Density (g/L)", "Glucose (g/L)", "Product (g/L)"], ...
    legend = ["Cells", "Glucose", "Product"], ...
    lineStyle  = ["-", "--", "-"], lineColor = ["k", "r", "b"], ...
    layout = [2 1], subplotIndex = [1, 2, 2], yAxis = ["L", "L", "R"], ...
    y_max = [inf, inf]);

% Experiment Length & Initial Conditions
t_max = 40;            % length of experiment/simulation run time (hours)
tspan = 0:0.01:t_max;  % time steps (h), 0 to t_max with step size 0.01
S0    = 5.0;           % initial conc of substrate (g/L)
X0    = 0.1;           % innoculation density of cell (g/L)
P0    = 0.01;          % initial product, (g/L)

% Cell and Product Model Constants
mumax = 0.25;     % max growth rate constant (1/h)
KS    = 1.0;     % half-saturation constant, g substrate per L (gS/L)
YXS   = 0.5;     % g cell per g substrate (gX/gS)
kd    = 0.0;      % death rate constant (1/h) 
KI    = 0.0;      % substrate inhibition (L/g S)^2
YPS   = 0.7;     % g product per g substrate (gP/gS)
alpha = 0.5;     % productivity, g product per g cell (gP/gX)
beta  = 0.01;    % productivity, g product per g cell-hr (gP/gX*h)
params = struct("mumax", mumax, "KS", KS, "YXS", YXS, ...
    "kd", kd, "KI", KI, "YPS", YPS, "alpha", alpha, "beta", beta);

% ODE solver (passes initial conditions to function)
options = odeset("NonNegative", [2,3]); % sets to zero when calc'ed <0
[t,y] = ode45(@(t,y) monod_w_product(t, y, params), tspan, [X0;S0;P0]);

% extract results to state variable names
[X, S, P] = deal(y(:,1), y(:,2), y(:,3));

% plot results
plot_results(t, y, plot_cfg)

%% Batch Vessel -- Monod Growth w/ Product
function [dydt] = monod_w_product(~, y, params)

% Unpack the state variables X, S, and P from y
[X, S, P] = deal(y(1), y(2), y(3));

% Unpack the parameters into individual variables
[mumax, KS, YXS, YPS, kd, KI, alpha, beta] = deal(params.mumax, params.KS, ...
    params.YXS, params.YPS, params.kd, params.KI,  params.alpha, params.beta);

%Cell Characterisitic Equations
mu    = mumax * S / (KS + S + KI * S^2); % substrate inhibition
qp = alpha * mu + beta; % mixed model

% Evaluate rate expressions -- Design (Batch) Equations
dXdt  = mu * X - kd * X;
dSdt = -mu * X / YXS - qp * X / YPS;
dPdt =  qp * X ;
dydt =  [dXdt; dSdt; dPdt];

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
