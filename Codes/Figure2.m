%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project : How noise shapes decision-making dynamics in fish schools: 
% Insights from a stochastic model with state-dependent diffusion
% Author  : Deze Liu, Daniel Burbano Lombana
% Lab     : The Swarm Intelligence Lab
% Date    : 08/01/2026
% Description :
% This script generates the results for Fig. 2a-d.
% It computes deterministic equilibrium-region maps for symmetric and
% asymmetric well configurations and plots representative noise-free
% trajectories of the reflected decision model on [0,xmax].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;
clear;
close all;

% Initial settings
RUN.a_b = true;  % Figure 2a,b: equilibrium-region maps
RUN.c_d = true;  % Figure 2c,d: deterministic trajectories

% Figure 2a,b parameters.
AB.xmax = 100;
AB.xlFixed = 0.01;
AB.dMin = 1e-3;
AB.dMaxPadding = 0.1;
AB.nD = 160;
AB.nR = 240;
AB.nBoundaryPoints = 1500;
AB.rootTolerance = 1e-8;
AB.stabilityTolerance = 1e-10;
AB.asymmetricRLimit = 10.0e5;
AB.symmetricPr = [-4.2e5, -1.7e5, 0.0e5, 4.2e5];
AB.symmetricPd = [72, 72, 72, 72];
AB.symmetricDxText = [0, 0, 0, 0];
AB.symmetricDyText = [3, 3, 3, 3];
AB.asymmetricPr = [-7.0e5, -4.0e5, -0.55e5, 4.5e5];
AB.asymmetricPd = [90, 65, 90, 65];
AB.asymmetricDxText = [0, 0, 0.6e4, 0];
AB.asymmetricDyText = [3, 3, 3, 3];

% Figure 2c,d parameters.
CD.xmax = 100;
CD.xlFixed = 0.01;
CD.beta = (0.1/125/2) * 0.3;
CD.dt = 1e-3;               % seconds
CD.simulationTime = 20;     % seconds
CD.initialState = 40;
CD.d = 70;
CD.plotStride = 1;

% Figure style
STYLE.figureSizeIn = [7.0 5.25];
STYLE.tickFontSize = 28;
STYLE.labelFontSize = 32;
STYLE.annotationFontSize = 22;

STYLE.axisLineWidth = 2.6;
STYLE.plotLineWidth = 3.2;
STYLE.outputDir = fullfile(fileparts(mfilename('fullpath')), ...
    'publication_figures_Figure2_new');
if ~exist(STYLE.outputDir, 'dir')
    mkdir(STYLE.outputDir);
end

% Symmetric equilibrium map
if RUN.a_b
    cfg = make_eq_config(true, AB, STYLE);
    plot_equilibrium_map(cfg);
    cfg = make_eq_config(false, AB, STYLE);
    plot_equilibrium_map(cfg);
end

% Representative time-series
if RUN.c_d
    plot_time_series_examples(true, CD, STYLE);
    plot_time_series_examples(false, CD, STYLE);
end





function cfg = make_eq_config(use_centered_well, ab, style)

    % Style parameters
    cfg.FS = style.tickFontSize;
    cfg.LabelFS = style.labelFontSize;
    cfg.AnnotationFS = style.annotationFontSize;
    cfg.LW = style.plotLineWidth;
    cfg.AxisLW = style.axisLineWidth;
    cfg.FigureSizeIn = style.figureSizeIn;
    cfg.OutputDir = style.outputDir;

    % Fixed parameters
    cfg.xmax = ab.xmax;
    cfg.xl_fixed = ab.xlFixed;
    cfg.rootTolerance = ab.rootTolerance;
    cfg.stabilityTolerance = ab.stabilityTolerance;
    cfg.use_centered_well = use_centered_well;
    if use_centered_well
        cfg.PanelName = 'Figure2a_symmetric_map';
    else
        cfg.PanelName = 'Figure2b_asymmetric_map';
    end

    % Parameter grid
    cfg.d_min = ab.dMin;
    cfg.d_max = ab.xmax - ab.dMaxPadding;

    cfg.Nd = ab.nD;
    cfg.Nr = ab.nR;

    cfg.d_vals = linspace(cfg.d_min, cfg.d_max, cfg.Nd);

    % Analytical boundary curves
    cfg.d_curve = linspace(cfg.d_min, cfg.d_max, ab.nBoundaryPoints);

    % Interior-equilibrium saddle-node thresholds
    cfg.r_SN_pos =  cfg.d_curve.^3 ./ (3*sqrt(3));
    cfg.r_SN_neg = -cfg.d_curve.^3 ./ (3*sqrt(3));

    % Well geometry
    if use_centered_well
        cfg.xl_curve = (ab.xmax - cfg.d_curve)/2;
        cfg.xh_curve = cfg.xl_curve + cfg.d_curve;
    else
        cfg.xl_curve = ab.xlFixed + 0*cfg.d_curve;
        cfg.xh_curve = cfg.xl_curve + cfg.d_curve;
    end

    % Boundary threshold for x = 0
    % alpha/beta = 2*x_l*x_h*(x_l+x_h)
    cfg.r_bdy_0 = 2 .* cfg.xl_curve .* cfg.xh_curve .* ...
        (cfg.xl_curve + cfg.xh_curve);

    % Boundary threshold for x = x_max
    % alpha/beta = -2*(xmax-xl)*(xmax-xh)*(2*xmax-xl-xh)
    cfg.r_bdy_xmax = -2 .* (ab.xmax - cfg.xl_curve) .* ...
        (ab.xmax - cfg.xh_curve) .* ...
        (2*ab.xmax - cfg.xl_curve - cfg.xh_curve);

    % Alpha/beta range
    if use_centered_well
        r_abs_max = max(abs([cfg.r_SN_pos, cfg.r_SN_neg, ...
            cfg.r_bdy_0, cfg.r_bdy_xmax]));
        cfg.r_lim = 1.05 * r_abs_max;
        cfg.use_axis_tight = true;
    else
        cfg.r_lim = ab.asymmetricRLimit;
        cfg.use_axis_tight = false;
    end

    cfg.r_vals = linspace(-cfg.r_lim, cfg.r_lim, cfg.Nr);

    % Representative points P1-P4
    if use_centered_well
        cfg.P_r = ab.symmetricPr;
        cfg.P_d = ab.symmetricPd;
        cfg.dy_text = ab.symmetricDyText;
        cfg.dx_text = ab.symmetricDxText;
    else
        cfg.P_r = ab.asymmetricPr;
        cfg.P_d = ab.asymmetricPd;
        cfg.dy_text = ab.asymmetricDyText;
        cfg.dx_text = ab.asymmetricDxText;
    end

    cfg.P_labels = {'P1','P2','P3','P4'};
end


function plot_equilibrium_map(cfg)

    xmax = cfg.xmax;
    xl_fixed = cfg.xl_fixed;
    use_centered_well = cfg.use_centered_well;

    % Theoretical classification from the deterministic drift polynomial.
    nIntEq = zeros(cfg.Nd, cfg.Nr);
    bdy0_num = false(cfg.Nd, cfg.Nr);
    bdyX_num = false(cfg.Nd, cfg.Nr);

    for ii = 1:cfg.Nd
        d = cfg.d_vals(ii);
        if use_centered_well
            xl = (xmax - d)/2;
            xh = xl + d;
        else
            xl = xl_fixed;
            xh = xl + d;
        end

        % F(x) = -r - 2(x-xl)(x-xh)(2x-xl-xh).
        driftWithoutBias = -2 * conv( ...
            conv([1, -xl], [1, -xh]), [2, -(xl+xh)]);

        for jj = 1:cfg.Nr
            driftPolynomial = driftWithoutBias;
            driftPolynomial(end) = driftPolynomial(end) - cfg.r_vals(jj);

            equilibriumRoots = roots(driftPolynomial);
            equilibriumRoots = real(equilibriumRoots( ...
                abs(imag(equilibriumRoots)) < cfg.rootTolerance));
            equilibriumRoots = equilibriumRoots( ...
                equilibriumRoots > cfg.rootTolerance & ...
                equilibriumRoots < xmax-cfg.rootTolerance);

            driftDerivative = polyder(driftPolynomial);
            stableInterior = polyval( ...
                driftDerivative, equilibriumRoots) < ...
                -cfg.stabilityTolerance;
            nIntEq(ii,jj) = sum(stableInterior);

            bdy0_num(ii,jj) = polyval(driftPolynomial, 0) < ...
                -cfg.stabilityTolerance;
            bdyX_num(ii,jj) = polyval(driftPolynomial, xmax) > ...
                cfg.stabilityTolerance;
        end
    end

    % class_map values:
    %   0 = none
    %   1 = 1 stable interior equilibrium only
    %   2 = 2 stable interior equilibria only
    %   3 = x = 0 stable boundary equilibrium only
    %   4 = x = x_max stable boundary equilibrium only
    %   5 = 1 stable interior equilibrium + x = 0 boundary equilibrium
    %   6 = 1 stable interior equilibrium + x = x_max boundary equilibrium

    class_map = zeros(cfg.Nd, cfg.Nr);

    % Interior-only regions
    class_map((nIntEq == 1) & ~bdy0_num & ~bdyX_num) = 1;
    class_map((nIntEq >= 2) & ~bdy0_num & ~bdyX_num) = 2;

    % Boundary-only regions
    class_map((nIntEq == 0) & bdy0_num & ~bdyX_num) = 3;
    class_map((nIntEq == 0) & ~bdy0_num & bdyX_num) = 4;

    % Coexistence regions:
    % one stable interior equilibrium + one stable boundary equilibrium
    class_map((nIntEq == 1) & bdy0_num & ~bdyX_num) = 5;
    class_map((nIntEq == 1) & ~bdy0_num & bdyX_num) = 6;

    % If two interior equilibria and a boundary are detected, keep the interior class.
    class_map((nIntEq >= 2) & (bdy0_num | bdyX_num)) = 2;

    % If both boundaries are detected and no interior equilibrium is detected,
    % choose x = 0 for display. This should be rare.
    class_map((nIntEq == 0) & bdy0_num & bdyX_num) = 3;

    % Custom colormap
    c_none = [1.0000 1.0000 1.0000];     % class 0, white

    c_int1 = [202 200 201]/255;          % class 1, gray
    c_int2 = [238 184 100]/255;          % class 2, yellow

    c_x0   = [143 181 216] / 255;        % class 3, light blue
    c_xmax = [201 116 104] / 255;        % class 4, light red

    % Coexistence colors
    c_int1_x0   = [0.1 0.35 0.9];        % class 5, deep blue
    c_int1_xmax = [0.55 0.10 0.10];      % class 6, deep red

    myColors = [
        c_none;        % 0
        c_int1;        % 1
        c_int2;        % 2
        c_x0;          % 3
        c_xmax;        % 4
        c_int1_x0;     % 5
        c_int1_xmax    % 6
    ];

    % Equilibrium-region map
    fig = publication_figure(cfg.FigureSizeIn);

    imagesc(cfg.r_vals, cfg.d_vals, class_map);
    set(gca,'YDir','normal');
    colormap(myColors);
    caxis([-0.5 6.5]);

    hold on;
    
    % Class 5: deep-blue region, vertical hatching
    overlay_class_hatching(gca, cfg.r_vals, cfg.d_vals, class_map, ...
        5, 'vertical', 8, [0 0 0], (2/3)*cfg.LW);
    
    % Class 6: deep-red region, horizontal hatching
    overlay_class_hatching(gca, cfg.r_vals, cfg.d_vals, class_map, ...
        6, 'horizontal', 6, [0 0 0], (2/3)*cfg.LW);
    
    % Theoretical boundaries
    add_theoretical_boundaries(cfg.r_SN_pos, cfg.r_SN_neg, ...
        cfg.r_bdy_0, cfg.r_bdy_xmax, cfg.d_curve, cfg.LW);

    % Mark representative points P1-P4
    plot(cfg.P_r, cfg.P_d, 'k*', ...
        'MarkerSize', 10, ...
        'LineWidth', 1.6, ...
        'HandleVisibility','off');

    for pp = 1:4
        text(cfg.P_r(pp) + cfg.dx_text(pp), ...
            cfg.P_d(pp) - cfg.dy_text(pp), cfg.P_labels{pp}, ...
            'Color','k', ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','top', ...
            'FontSize', cfg.AnnotationFS, ...
            'FontWeight','normal', ...
            'Interpreter','none');
    end

    hold off;

    ax = gca;
    ax.FontSize = cfg.FS;
    ax.FontWeight = 'normal';
    ax.LineWidth = cfg.AxisLW;
    ax.TickDir = 'out';

    xlabel('$\alpha/\beta$','Interpreter','latex', ...
        'FontSize',cfg.LabelFS,'FontWeight','normal');
    ylabel('$d=x_h-x_l$','Interpreter','latex', ...
        'FontSize',cfg.LabelFS,'FontWeight','normal');

    if cfg.use_axis_tight
        axis tight;
    else
        xlim([-cfg.r_lim, cfg.r_lim]);
        ylim([cfg.d_min, cfg.d_max]);
    end

    box on;
    grid on;

    export_publication_panel(fig, cfg.OutputDir, cfg.PanelName);
end


function add_theoretical_boundaries(r_SN_pos, r_SN_neg, ...
    r_bdy_0, r_bdy_xmax, d_curve, LW)

    % Black dashed: theoretical saddle-node boundary
    plot(r_SN_pos, d_curve, 'k--', 'LineWidth', LW+0.4);
    plot(r_SN_neg, d_curve, 'k--', 'LineWidth', LW+0.4);

    % White dashed: theoretical boundary-equilibrium thresholds
    plot(r_bdy_0, d_curve, 'w--', 'LineWidth', LW+0.4);
    plot(r_bdy_xmax, d_curve, 'w--', 'LineWidth', LW+0.4);
end


function plot_time_series_examples(use_centered_well, cd, style)

    FS = style.tickFontSize;
    LW = style.plotLineWidth;
    xmax = cd.xmax;
    xl_fixed = cd.xlFixed;
    beta = cd.beta;
    dt = cd.dt;
    d0 = cd.d;

    if use_centered_well
        xl = (xmax - d0)/2;
        xh = xl + d0;
    else
        xl = xl_fixed;
        xh = xl + d0;
    end

    % Thresholds
    r_SN = d0^3 / (3*sqrt(3));

    % Boundary threshold for x = 0
    r_bdy_0 = 2 * xl * xh * (xl + xh);

    % Boundary threshold for x = x_max
    r_bdy_xmax = -2 * (xmax - xl) * (xmax - xh) * ...
                      (2*xmax - xl - xh);

    % Four representative alpha/beta values

    % Case 1: two stable interior equilibria
    r_two_int = 0;

    % Case 2: one stable interior equilibrium, no boundary equilibrium
    gap_pos = r_bdy_0 - r_SN;          % interval: (r_SN, r_bdy_0)
    gap_neg = -r_SN - r_bdy_xmax;      % interval: (r_bdy_xmax, -r_SN)

    if gap_pos > gap_neg && gap_pos > 0
        % positive side has a valid one-interior region
        r_one_int = 0.5 * (r_SN + r_bdy_0);
    else
        % negative side has a valid one-interior region
        r_one_int = 0.5 * (r_bdy_xmax - r_SN);
    end

    % Case 3: boundary-only region for x = 0
    r_x0 = 1.15 * max(r_bdy_0, r_SN);

    % Case 4: boundary-only region for x = x_max
    r_xmax = 1.15 * min(r_bdy_xmax, -r_SN);

    % P1: light red, P2: light blue, P3: gray, P4: yellow
    r_cases = [r_xmax, r_one_int, r_two_int, r_x0];

    labels = {
        'P1', ...
        'P2', ...
        'P3', ...
        'P4'
    };

    % Colors match the equilibrium maps.
    c_gray   = [202 200 201]/255;      % 1 stable interior equilibrium
    c_yellow = [238 184 100]/255;      % 2 stable interior equilibria
    c_blue   = [143 181 216]/255;      % x = 0 boundary equilibrium
    c_red    = [201 116 104]/255;      % x = x_max boundary equilibrium

    % P1 red, P2 blue, P3 gray, P4 yellow
    case_colors = [
        c_red;
        c_gray;
        c_yellow;
        c_blue
    ];

    % Physical time settings are defined together at the top of the script.
    Tseries = cd.simulationTime;
    Nt = round(Tseries / dt);
    t_sec = (0:Nt) * dt;

    x0 = cd.initialState;

    Xtraj = zeros(Nt+1, length(r_cases));
    Xtraj(1,:) = x0;

    % Euler integration
    for m = 1:length(r_cases)

        r = r_cases(m);
        x = x0;

        for k = 1:Nt

            F = -r - 2 * (x - xl) * (x - xh) * (2*x - xl - xh);

            x = x + beta * F * dt;

            x = min(max(x, 0), xmax);

            Xtraj(k+1,m) = x;
        end
    end

    % Plot four time series
    fig = publication_figure(style.figureSizeIn);
    hold on;

    plot_stride = cd.plotStride;
    idx_plot = 1:plot_stride:length(t_sec);

    for m = 1:length(r_cases)
        plot(t_sec(idx_plot), Xtraj(idx_plot,m), ...
            'LineWidth', LW, ...
            'Color', case_colors(m,:));
    end

    % Reference lines for boundaries
    yline(0, 'k:', 'LineWidth', style.axisLineWidth);
    yline(xmax, 'k:', 'LineWidth', style.axisLineWidth);
    ylim([-5 105]);
    hold off;

    ax = gca;
    ax.FontSize = FS;
    ax.FontWeight = 'normal';
    ax.LineWidth = style.axisLineWidth;
    ax.TickDir = 'out';

    xlabel('Time (s)','Interpreter','latex', ...
        'FontSize',style.labelFontSize,'FontWeight','normal');
    ylabel('$x(t)$','Interpreter','latex', ...
        'FontSize',style.labelFontSize,'FontWeight','normal');

    legend(labels, ...
        'Interpreter','none', ...
        'Location','best', ...
        'FontSize',style.annotationFontSize);

    box on;
    grid on;

    if use_centered_well
        panelName = 'Figure2c_symmetric_time_series';
    else
        panelName = 'Figure2d_asymmetric_time_series';
    end
    export_publication_panel(fig, style.outputDir, panelName);
end


function fig = publication_figure(figureSizeIn)
    fig = figure('Color', 'w', 'Units', 'inches', ...
        'Position', [1 1 figureSizeIn], 'PaperPositionMode', 'auto');
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 figureSizeIn];
    fig.PaperSize = figureSizeIn;
end


function export_publication_panel(fig, outputDir, panelName)
    drawnow;
    print(fig, fullfile(outputDir, [panelName '.png']), '-dpng', '-r300');
    print(fig, fullfile(outputDir, [panelName '.pdf']), '-dpdf', '-vector');
end


function overlay_class_hatching(ax, xValues, yValues, classMap, ...
        classValue, direction, stride, lineColor, lineWidth)

    mask = (classMap == classValue);

    switch lower(direction)

        case 'vertical'
            for jj = 1:stride:numel(xValues)

                runs = contiguous_true_runs(mask(:,jj));

                for kk = 1:size(runs,1)
                    plot(ax, ...
                        [xValues(jj), xValues(jj)], ...
                        [yValues(runs(kk,1)), yValues(runs(kk,2))], ...
                        '-', ...
                        'Color', lineColor, ...
                        'LineWidth', lineWidth, ...
                        'HandleVisibility', 'off');
                end
            end

        case 'horizontal'
            for ii = 1:stride:numel(yValues)

                runs = contiguous_true_runs(mask(ii,:).');

                for kk = 1:size(runs,1)
                    plot(ax, ...
                        [xValues(runs(kk,1)), xValues(runs(kk,2))], ...
                        [yValues(ii), yValues(ii)], ...
                        '-', ...
                        'Color', lineColor, ...
                        'LineWidth', lineWidth, ...
                        'HandleVisibility', 'off');
                end
            end

        otherwise
            error('Unknown hatching direction: %s', direction);
    end
end


function runs = contiguous_true_runs(maskVector)

    transitions = diff([false; maskVector(:); false]);

    runStarts = find(transitions == 1);
    runEnds   = find(transitions == -1) - 1;

    runs = [runStarts, runEnds];
end
