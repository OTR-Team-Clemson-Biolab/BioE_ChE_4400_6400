% Continuous Stirred-Tank Reactor -- Monod Model with Product.
% Simulates and plots a CSTR vessel with monod growth and
% product formation for various dilutions.

clc
clear
close all   % comment to overlay results on existing figure

% Main Program

% Plot Configuration for Labeling
plot_cfg = struct( xlabel = "Dilution, D (1/h)", ...
    ylabel = ["Concentrations (g/L)", "Concentrations (g/L)", "Concentrations (g/L)", ...
              "Productivities (g/h)", "Productivities (g/h)"], ...
    legend = ["Substrate", "Cell", "Product", "Cell Productivity (DX)", "Productivity (PD)"], ...
    lineStyle  = ["-", "-", "-", "-", "-"], lineColor = ["r", "b", "g", "b", "g"], ...
    layout = [2 1], subplotIndex = [1, 1, 1, 2, 2], yAxis = ["L", "L", "L", "L", "L"], ...
    y_max = [inf, inf]);

% Cell and Product Model Constants
mumax = 0.50;       % max growth rate constant; 1/h
YXS = 0.40;         % gX/gS, g cell per g substrate 
YPS = 0.15;         % gP/gx, g product per g dry cell
alpha = 0.1;        % productivity; gP per gX*h
beta = 0.01;        % productivity; 1/h
KS = 0.10;          % gS/L, g substrate per L
kd = 0;             % death rate constant; 1/h

% Inlet Concentrations & Maximum Dilution
SF = 5.0;            % g/L, inlet substrate concentration
D_max = 0.99*mumax;  % maximum dilution factor (g/l) -- do not go higher than mumax

% Model Equations for Cell and Product
D = linspace(0.00, D_max, 400)';   % dilution rate column vector
mu = D + kd;                       % growth rate
S = (mu * KS) ./ (mumax - mu);     % substrate concentration
qp = alpha * mu + beta;            % product formation rate

% Cell and Product Concentrations
X_num = D .* (SF - S);             % cell concentration numerator
X_den = ((mu / YXS) + (qp / YPS)); % cell concentration denominator 
X =  X_num ./ X_den;               % cell concentration
P = qp .* X ./ D;                  % product concentration
DX = D .* X;                       % cell productivity
PX = D .* P;                       % product productivity

% Plot results
y = [S, X, P, DX, PX];
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
