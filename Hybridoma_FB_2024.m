%% Hybridoma Fed-Batch Simulation based on Tatiraju et al., 1999.

% Tatiraju S, Soroush M, Mutharasan R. 1999. 
% Multi-rate nonlinear state and parameter estimation in a bioreactor. 
% Biotechnology and Bioengineering 63(1):22-32.

% Notes:
% Orginal model required the use of on-line OTR measurments as a third
% substrate.  Reworked to fit batch data. Higher kd needed and higher
% mumax. Production rate dependence on glucose was also added.

% converted from "qp=gamma1 - gamma2*mu" to "qp=alpha*mu + beta", 
% since negative gamma2 not physically meaningful

% RECALL, product changes after substrate reaches zero or goes negative,
% are not real; however, production is still possible if fed on and and S
% low Initial set up is a batch. To get FB, provide a value and edit code
% for PUMPON and flowrate (F) (line 70 and lines 107-109). Will also need
% to provide a value for muset (line 55).

clc
clear
close all   % Comment to overlay results on existing figure

%% Plot Configuration
plot_cfg = struct( xlabel = "Time (h)", ...
    ylabel = ["Concentrations", "Concentrations", "Concentrations", "Concentrations", "Concentrations", "Volume"], ...
    legend = ["Viable Cells (cells/mL)", "Total Cells (cells/mL)", ...
                "Glutamine (mM)", "Glucose (mM)", "MAb (mg/L)", "Volume (L)"], ...
    lineStyle = ["--", "-", "-", "--", "--", "--"], ...
    lineColor = ["k",    "g",  "m",   "b",   "r",   "c"], ...
    layout = [2, 2], subplotIndex = [1, 1, 2, 2, 3, 4], ...
    yAxis = ["L", "L", "L", "L", "L", "L"], ...
    y_max = [inf, inf, inf, inf, inf, inf]); 

%% Experiment Length & Initial Conditions
tmax  = 200;           % hours
tspan = 0:0.01:tmax;   % time steps (hours)
X0    = 1.0e5;         % initial cells/mL
S10   = 3.0;           % initial glutamine conc (mM)
S20   = 20;            % initial glucose conc (mM)
P0    = 1.0;           % initial product conc (mg/L)
XT0   = 1.0e5;         % initial total cells (cells/mL)
V0    = 100.0;         % initial volume (L)

% Parameters - Reaction-rate constants

% Growth yield coefficients
YXS1   = 6.16e5;   % cell/mL per mM of glutamine
YXS2   = 1.5e5;    % cell/mL per mM of glucose

% Growth and death rate constants
mumax  = 0.044;    % maximum growth rate (1/h)
kd     = 0.008;    % death rate (1/h)
muset  = NaN ;     % fraction of mumax for muset

% Degradation and saturation constants
k1     = 4.8e-3;   % natural degradation rate of glutamine (1/h)
KS     = 0.10;     % half-saturation constant for glucose (mM)

% Product formation rate constants
beta   = 5.0e-8;   % basal production rate of MAb (mg*mL per cell*h)
alpha  = 3.4e-5;   % production rate coefficient dependent on growth (mg*mL/cell)

% Feed concentrations and initial setup values
SF1    = 50;       % feed concentration of glutamine (mM)
SF2    = 250;      % feed concentration of glucose (mM)

% Pump activation time
PUMPON = 20;      % time to start pump (h)

% Creating the params structure
params = struct(...
    "YXS1",   YXS1,   "YXS2",   YXS2,   "mumax",  mumax, ...
    "kd",     kd,     "muset",  muset,  "k1",     k1,     "KS",     KS, ...
    "beta",   beta,   "alpha",  alpha,  "SF1",    SF1, ...
    "SF2",    SF2,    "S20",    S20,    "V0",     V0, ...
    "PUMPON", PUMPON);

%% ODE solver
[t, y] = ode15s(@(t, y) fedbatch_model(t, y, params), tspan, [X0; S10; S20; P0; XT0; V0]);

%% Extract results
[X, S1, S2, P, XT, V] = deal(y(:, 1), y(:, 2), y(:, 3), y(:, 4), y(:, 5), y(:, 6));

% Plot results
plot_results(t, [X, XT, S1, S2, P, V], plot_cfg);
cfg=plot_cfg;
x=t;
y=[X, XT, S1, S2, P, V];

%% Hybridoma Fed-Batch Model Function
function dydt = fedbatch_model(t, y, params)
    % Unpack state variables
    [X, S1, S2, P, XT, V] = deal(y(1), y(2), y(3), y(4), y(5), y(6));
    
    % Unpack parameters
    [YXS1, YXS2, mumax, kd, muset, k1, KS, beta, alpha, SF1, SF2, S20, V0, PUMPON] = ...
        deal(params.YXS1, params.YXS2, params.mumax, params.kd, params.muset, params.k1, params.KS, ...
             params.beta, params.alpha, params.SF1, params.SF2, params.S20, ...
             params.V0, params.PUMPON);

    % Set feed rate F based on time t
    if t < PUMPON
        F = 0.0;
    elseif t >= PUMPON
        F = 0; % batch, change to turn on as constant pump
        %  F = 0.0; % constant flow rate pump option
        %  F = (muset*S20*V0/SF2)*exp(muset*(t - PUMPON)); % exp feed pump
    end

    % Growth and production rates
    mu = mumax * S2 / (KS + S2);
    qp = beta + alpha * mu;

    % Model equations
    dXdt = mu * X - kd * X - X * F / V;
    dXTdt = mu * X - XT * F / V;
    dS1dt = (SF1 - S1) * F / V - mu * X / YXS1 - k1 * S1; % natural degradation included
    dS2dt = (SF2 - S2) * F / V - mu * X / YXS2;  % neglecting maintainance
    dPdt = qp * X - P * F / V;
    dVdt = F;

    dydt = [dXdt; dS1dt; dS2dt; dPdt; dXTdt; dVdt];
    
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

%% Notes

