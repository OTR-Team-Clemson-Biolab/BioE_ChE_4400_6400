% Continuous Stirred-Tank Reactor with Recycle -- Monod Model with Product.
% Simulates and plots a CSTR vessel with monod growth and product formation
% with a recycle stream for various dilutions.

clc
clear
close all   % comment to overlay results on existing figure

% Main Program

% Plot Configuration for Labeling
plot_cfg = struct( xlabel = "Dilution, D (1/h)", ...
    ylabel = ["Concentrations (g/L)", "Concentrations (g/L)", "Concentrations (g/L)", ...
               "Concentrations (g/L)", "Productivities (g/L*h)", "Productivities (g/L*h)"], ...
    legend = ["Substrate", "Cell (Xinside)", "Cell (Xexit)", "Product", ...
               "Productivity (DXexit)", "Productivity (PD)"], ...
    lineStyle  = ["-", "-", "-", "-", "-", "-"], lineColor = ["r", "m", "b", "g", "b", "g"], ...
    layout = [2 1], subplotIndex = [1, 1, 1, 1, 2, 2], yAxis = ["L", "L", "L", "L", "L", "L"], ...
    y_max = [inf, inf]);

% Cell and Product Model Constants
mumax = 0.5;        % max growth rate constant; 1/h
YXS = 0.5;          % gX/gS, g cell per g substrate 
YPS = 0.70;         % gP/gS, g product per g substrate
alpha = 0.1;        % productivity; gP per gX*h
beta = 0.01;        % productivity; 1/h
KS = 0.1;           % gS/L, g substrate per L
kd = 0.0;           % death rate constant; 1/h 

% Inlet Concentrations & Maximum Dilution
SF = 10.0;           % g/L, inlet substrate concentration
D_max = 0.99*mumax;  % maximum dilution factor (g/l) -- do not go higher than mumax

% Recycle Stream Parameters
alphaR = 0.1;       % recycle ratio; (= 0 for standard CSTR)
Conc = 1.1;         % concentration factor; (= 1 for standard CSTR)

% Dilution Rate and Substrate Concentration Calculations
D = linspace(0.0, D_max, 300)';        % Dilution rate column vector
mu = D * (1 + alphaR - alphaR * Conc); % Effective specific growth rate (adjusted for recycle)
S = (mu * KS) ./ (mumax - mu);         % Substrate concentration (Monod model)
qp = mu * alpha + beta;                % Specific product formation rate

% Cell and Product Concentrations
X_num = D .* (SF - S);                     % cell concentration numerator
X_den = ((mu / YXS) + (qp / YPS));         % cell concentration denominator 
X =  X_num ./ X_den;                       % cell concentration
Xexit = X .* (1 + alphaR - alphaR * Conc); % Cell concentration at the reactor exit
P = qp .* X ./ D;                          % Product concentration
DX = D .* X;                               % Cell productivity inside the reactor
DXexit = D .* Xexit;                       % Cell productivity at the reactor exit
PD = D .* P;                               % Product productivity

% Plot Results
y = [S, X, Xexit, P, DXexit, PD];
plot_results(D, y, plot_cfg)

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
