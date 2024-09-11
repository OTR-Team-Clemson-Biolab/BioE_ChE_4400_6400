% Monod Model
% Simulates and plots a batch vessel with monod growth 

clc
clear 
close all   % comment to overlay results on existing figure

%% Main Program

% Plot Configuration for Labeling
plot_cfg = struct( xlabel = "Time (h)", ...
    ylabel = ["Concentration (g/L)", "Concentration (g/L)"], ...
    legend = ["Cells", "Glucose"], ...
    lineStyle  = ["-", "--"], lineColor = ["k", "r"], ...
    layout = [2 1], subplotIndex = [1, 2], yAxis = ["L", "L"], ...
    y_max = [inf, inf]);

% Experiment Length & Initial Conditions
t_max = 50;            % length of experiment/simulation run time (hours)
tspan = 0:0.01:t_max;  % time steps (h), 0 to t_max with step size 0.01
S0    = 50.0;          % g/L, initial conc of substrate
X0    = 0.1;           % g/L, innoculation density of cell

% Cell and Product Model Constants
mumax = 0.25;     % max growth rate constant (1/h)
KS    = 0.5;      % half-saturation constant (gS/L) g substrate per L
YXS   = 0.48;     % (gX/gS), g cell per g substrate 
kd    = 0.0;      % death rate constant (1/h) 
KI    = 0.0;      % substrate inhibition (L/g S)^2
params = struct("mumax", mumax, "KS", KS, "YXS", YXS, "kd", kd, "KI", KI);

% ODE solver (passes initial conditions to function)
options = odeset("NonNegative", [2,3]); % sets to zero when calc'ed <0
[t,y] = ode45(@(t,y) monod(t,y,params), tspan, [X0; S0]); 

% extract results to state variable names
[X, S] = deal(y(:,1), y(:,2));

% plot results
plot_results(t, y, plot_cfg);

%% Batch Vessel -- Monod Growth 
function [dydt] = monod(~, y, params)

% Unpack the state variables X and S from y
[X, S] = deal(y(1), y(2)); 

% Unpack the parameters into individual variables
[mumax, KS, YXS, kd, KI] = deal(params.mumax, params.KS, params.YXS, params.kd, params.KI);

% Model Equations - Cell Characteristics and Product Equations
mu    = mumax * S / (KS + S + KI * S^2); % substrate inhibition

% Design (Batch) Equations
dXdt  = mu * X - kd * X;
dSdt  = -mu * X / YXS;
dydt  = [dXdt; dSdt];

end
%% Notes on ODE Solver: 
%
% options = odeset("NonNegative", [2,3]); 
%
% [t, y] = ode45(@(t,y) monod(t,y,params), tspan, [X0; S0]); 
% 
% - ode45 (4th and 5th order Runge-Kutta) is suitable for many problems.
%
% - ode23 (2nd and 3rd order Runge-Kutta) & ode23t (trapezoidal rule) are
%   better for specific cases, such as moderately stiff or difficult equations.
%
% NonNegative settings: options = odeset('NonNegative', [2,3]) is used to
% specify that certain components of the solution should be constrained to
% be non-negative, e.g., substrate concentrations cannot be negative. The
% [2,3] syntax refers to the 2nd and 3rd components of the solution vector y.
% However, in this example, there are only two components in y ([X0; S0]),
% corresponding to cell density and substrate concentration, respectively.
% Since there isn't a 3rd component, the 3 is superfluous. Therefore,
% using odeset('NonNegative', 2) would be
% appropriate and would give the same result.
%
% Function Handle: @(t,y) monod(t,y,params) creates an anonymous function
% that passes the current time t, the current state vector y (containing
% the current values of X(t) and S(t)), and the constants (stored in params)
% to the monod function.
%
% tspan: Defines the time interval over which to integrate the ODEs. For
% example, tspan = 0:0.01:t_max means the solver will integrate from time
% 0 to t_max, with time steps of 0.01. Note that the solver may not
% strictly adhere to these steps, as it uses variable step sizes internally.
%
% [X0; S0]: Initial values of X(t) and S(t) at time t = 0, provided as a
% column vector.
%
% [t, y]: The returned `t` is a vector of time points chosen by the ODE
% solver. Although `t` often aligns with `tspan`, the solver uses a
% variable step size, so `t` may contain additional or fewer points. Always
% use `t` for plotting, as it reflects the actual time points at which the
% solution was computed.
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
