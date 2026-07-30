%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project : How noise shapes decision-making dynamics in fish schools: 
% Insights from a stochastic model with state-dependent diffusion
% Author  : Deze Liu, Daniel Burbano Lombana
% Lab     : The Swarm Intelligence Lab
% Date    : 08/01/2026
% Description :
% This script generates the numerical results for Fig. 3-6.
% It analyzes noise-induced stationary behavior in a reflected stochastic
% decision model with nonlinear drift and either state-dependent or constant
% diffusion. The script includes four switchable modules:
%   Fig. 3a,b : stationary densities and comparison of density peaks with
%               deterministic stable equilibria;
%   Fig. 4c,d : large-noise limits for constant and scaled state-dependent
%               diffusion;
%   Fig. 5a-d : peak-configuration maps in the (alpha/beta,d) parameter
%               plane for state-dependent and constant diffusion;
%               Monte Carlo stationary-density heat maps as alpha/beta varies;
%   Fig. 6a-d : stationary-density extrema versus noise scale kappa for four
%               drift-bias values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;
clear;
close all;

% Initial settings
RUN.a_b = true;
RUN.c_d = true;
RUN.e_f = true;
RUN.g_h = true;
RUN.i_l = true;

% Figure 3a,b: stationary densities.
AB.xmax = 100;
AB.beta = (0.1/125/2) * 0.3;
AB.b2 = -1.5e-2;
AB.b1 = 1.50;
AB.b0 = 5.0;
AB.useCenteredWell = true;
AB.xlFixed = 0.01;
AB.dt = 1e-3;
AB.simulationTime = 600;
AB.burninTime = AB.simulationTime/3;
AB.nTrajectories = 200;
AB.initialState = 50;
AB.nTheoryGrid = 2000;
AB.nBins = 100;
AB.thinStep = 20;
AB.randomSeed = 1;
AB.d = 40;
AB.rRatios = [0, 3];  % r/(d^3/(3*sqrt(3))): bimodal, unimodal

% Figure 4a,b: large-noise limits.
CD.xmax = 100;
CD.xl = 20;
CD.xh = 80;
CD.beta = 0.1/125/2;
CD.b2 = -1e-4;
CD.b1 = -CD.b2*CD.xmax;
CD.b0 = 10;
CD.nX = 2000;
CD.alphaOverBeta = 0;
CD.runSpecifiedKappa = true;
CD.kappaValues = [10, 15, 30];
CD.curveColors = [0.85 0.1 0.1; 0 0 0; 0.1 0.35 0.9];

% Figure 5a,b: peak-configuration maps.
EF.xmax = 100;
EF.beta = (0.1/125/2) * 0.3;
EF.b2 = -1.5e-2;
EF.b1 = 1.50;
EF.b0 = 5.0;
EF.useCenteredWell = true;
EF.xlFixed = 0.01;
EF.dMin = 1e-3;
EF.dMaxPadding = 0.1;
EF.nD = 220;
EF.nR = 300;
EF.rLimitFactor = 2.6;
EF.nDiffusionCheckPoints = 4000;
EF.nBoundaryPoints = 2000;

% Figure 5c,d: stationary-density heat maps.
GH.xmax = 100;
GH.beta = (0.1/125/2) * 0.3;
GH.b2 = -1.5e-2;
GH.b1 = 1.50;
GH.b0 = 5.0;
GH.useCenteredWell = true;
GH.xlFixed = 0.01;
GH.d = 40;
GH.nR = 101;
GH.rRangeInSN = [-1.5, 1.5];
GH.dt = 1e-3;
GH.simulationTime = 2400;
GH.burninTime = 800;
GH.nTrajectories = 200;
GH.thinStep = 20;
GH.initialState = 50;
GH.nBins = 300;
GH.smoothBins = 31;
GH.nPeakGrid = 1500;
GH.randomSeed = 1;

% Figure 6a-d: extrema branches
IL.xmax = 100;
IL.xl = 20;
IL.xh = 80;
IL.beta = 0.1/125/2;
IL.b2 = -1e-4;
IL.b1 = -IL.b2*IL.xmax;
IL.b0 = 10;
IL.kappaMin = 0.05;
IL.kappaMax = 30;
IL.nKappa = 1200;
IL.rRatios = [0, 0.75, 1.50, 2.25];
IL.kappaReferenceLines = [10, 15, 30];
IL.nPositivityCheckPoints = 2001;
IL.maximumColor = [0.1 0.35 0.9];
IL.minimumColor = [0.85 0.1 0.1];
IL.maximumMarkerSize = 10;
IL.minimumMarkerSize = 5;
IL.foldMarkerSize = 10;

% Figure style
STYLE.figureSizeIn = [7.0 5.25];
STYLE.tickFontSize = 28;
STYLE.labelFontSize = 32;
STYLE.annotationFontSize = 22;
STYLE.axisLineWidth = 2.6;
STYLE.plotLineWidth = 3.2;


%% Fig. 3a,b: Stationary densities
if RUN.a_b
    FS = STYLE.tickFontSize;
    LW = STYLE.plotLineWidth;

    % Parameters
    xmax = AB.xmax;
    beta = AB.beta;

    % Base diffusion coefficients
    b2_sd = AB.b2;
    b1_sd = AB.b1;
    b0_sd = AB.b0;

    use_centered_well = AB.useCenteredWell;
    xl_fixed = AB.xlFixed;

    % Simulation parameters
    dt        = AB.dt;
    T         = AB.simulationTime;
    T_burn    = AB.burninTime;
    Ntraj     = AB.nTrajectories;
    x0        = AB.initialState;
    Nx_grid   = AB.nTheoryGrid;
    Nbins     = AB.nBins;
    thin_step = AB.thinStep;

    rng(AB.randomSeed);

    % Colors
    c_hist    = [143 181 216]/255; % light blue
    c_statepk = [0.85 0.1 0.1];    % red
    c_constpk = [0.1 0.35 0.9]/1;  % blue
    c_eq      = [0 0 0];           % equilibrium circles: black

    d0 = AB.d;

    % deterministic saddle-node boundary: |r| = d^3 / (3*sqrt(3))
    rb = d0^3 / (3*sqrt(3));

    r_bi   = AB.rRatios(1) * rb;
    r_mono = AB.rRatios(2) * rb;

    regimes = cell(1,2);
    regimes{1} = struct('name','Bimodal',   'r',r_bi);
    regimes{2} = struct('name','Unimodal', 'r',r_mono);

    % Quantitative goodness-of-fit results for Fig. 3(a,b)
    nRegimes = numel(regimes);
    
    regimeNames_GOF    = strings(nRegimes,1);
    rValues_GOF        = zeros(nRegimes,1);
    L1_distance_GOF    = zeros(nRegimes,1);
    retainedSamples_GOF = zeros(nRegimes,1);

    % Plot
    for ireg = 1:numel(regimes)
        R = regimes{ireg};

        % same drift, two different diffusion choices
        C_state = peak_equilibrium_comparison__make_case([R.name ' + state diffusion'], ...
            R.r, d0, beta, xmax, use_centered_well, xl_fixed, ...
            'state', b0_sd, b1_sd, b2_sd);

        C_const = peak_equilibrium_comparison__make_case([R.name ' + constant diffusion'], ...
            R.r, d0, beta, xmax, use_centered_well, xl_fixed, ...
            'constant', b0_sd, b1_sd, b2_sd);

        % x-grid for theoretical stationary density
        xg = linspace(0, xmax, Nx_grid);

        % theoretical stationary densities
        rho_state = peak_equilibrium_comparison__stationary_density_general(xg, C_state.f, C_state.g);
        rho_const = peak_equilibrium_comparison__stationary_density_general(xg, C_const.f, C_const.g);

        % find theoretical peaks
        [xpk_state, ~] = peak_equilibrium_comparison__theoretical_peaks(xg, rho_state);
        [xpk_const, ~] = peak_equilibrium_comparison__theoretical_peaks(xg, rho_const);

        % Monte Carlo samples from reflected SDE, state-dependent diffusion only
        samples = peak_equilibrium_comparison__simulate_reflected_rsde( ...
            C_state.f, C_state.g, x0, xmax, ...
            dt, T, T_burn, Ntraj, thin_step);
        
        % L1 goodness-of-fit distance
        edges = linspace(0, xmax, Nbins+1);
        
        % Empirical probability contained in each histogram bin
        empiricalCounts = histcounts(samples, edges);
        empiricalProbability = empiricalCounts(:) / sum(empiricalCounts);
        
        % Theoretical CDF obtained from the stationary density
        cdfState = cumtrapz(xg, rho_state);
        cdfState = cdfState / cdfState(end);
        
        % Enforce exact CDF endpoint values
        cdfState(1)   = 0;
        cdfState(end) = 1;
        
        % Theoretical CDF evaluated at the histogram-bin edges
        theoreticalCdfAtEdges = interp1( ...
            xg, cdfState, edges, 'linear', 'extrap');
        
        theoreticalCdfAtEdges = min( ...
            max(theoreticalCdfAtEdges,0),1);
        
        % Theoretical probability contained in each bin
        theoreticalProbability = diff(theoreticalCdfAtEdges(:));
        
        theoreticalProbability = max(theoreticalProbability,0);
        theoreticalProbability = theoreticalProbability / ...
            sum(theoreticalProbability);
        
        % L1 distance between empirical and theoretical bin probabilities
        L1_distance_GOF(ireg) = sum(abs( ...
            empiricalProbability - theoreticalProbability));
        
        % Save settings and identifiers
        regimeNames_GOF(ireg)     = string(R.name);
        rValues_GOF(ireg)         = R.r;
        retainedSamples_GOF(ireg) = numel(samples);
        
        fprintf('%-9s: L1 = %.6g, retained observations = %d\n', ...
            R.name, L1_distance_GOF(ireg), ...
            retainedSamples_GOF(ireg));
        
        % Histogram used in the figure
        [counts, edges] = histcounts( ...
            samples, edges, 'Normalization', 'pdf');
        
        centers = 0.5*(edges(1:end-1) + edges(2:end));

        % deterministic equilibria
        xeq_all = peak_equilibrium_comparison__deterministic_equilibria(C_state.alpha, C_state.beta, ...
        C_state.xl, C_state.xh, xmax);

        xeq_det = peak_equilibrium_comparison__keep_stable_equilibria(xeq_all, C_state.f, xmax);

        % Plot
        fig = publication_figure(STYLE.figureSizeIn);

        % histogram: blue
        % empirical histogram: blue
        h_hist = bar(centers, counts, 1.0, ...
            'FaceColor', c_hist, ...
            'FaceAlpha', 0.28, ...
            'EdgeColor', 'none');
        hold on;
        
        % theoretical stationary probability density: black
        h_theory = plot(xg, rho_state, 'k-', ...
            'LineWidth', LW);

        yl = ylim;

        % state-dependent diffusion peaks: red dashed lines + red triangles
        for j = 1:numel(xpk_state)
            xline(xpk_state(j), '--', ...
                'Color', c_statepk, ...
                'LineWidth', LW);

            y_state_here = interp1(xg, rho_state, xpk_state(j), 'linear');
            plot(xpk_state(j), y_state_here, '^', ...
                'Color', c_statepk, ...
                'MarkerFaceColor', c_statepk, ...
                'MarkerSize', 10, ...
                'LineWidth', LW);
        end

        % constant-diffusion peaks: gray dashed lines + gray triangles
        for j = 1:numel(xpk_const)
            xline(xpk_const(j), '--', ...
                'Color', c_constpk, ...
                'LineWidth', LW);

            y_const_here = interp1(xg, rho_state, xpk_const(j), 'linear');
            plot(xpk_const(j), y_const_here, '^', ...
                'Color', c_constpk, ...
                'MarkerFaceColor', c_constpk, ...
                'MarkerSize', 10, ...
                'LineWidth', LW);
        end

        % deterministic equilibria: white circles
        if ~isempty(xeq_det)
            yeq = interp1(xg, rho_state, xeq_det, 'linear');
            plot(xeq_det, yeq, 'o', ...
                'Color', c_eq, ...
                'MarkerFaceColor', 'w', ...
                'MarkerSize', 9, ...
                'LineWidth', LW);
        end

        ylim(yl);

        xlabel('$x$','Interpreter','latex');
        ylabel('Probability density');
        xlim([0 xmax]);
        grid on;
        box on;
        style_publication_axes(gca, STYLE);

        % legend handles
        legend([h_theory, h_hist], ...
            {'$\rho^*$', ...
             '$\rho_e$'}, ...
            'Interpreter', 'latex', ...
            'Location', 'best', ...
            'FontSize', STYLE.annotationFontSize);
    end

        % Summary table for Fig. 3(a,b)
        Figure3ab_L1_results = table( ...
            regimeNames_GOF, ...
            rValues_GOF, ...
            repmat(T,nRegimes,1), ...
            repmat(Ntraj,nRegimes,1), ...
            repmat(T_burn,nRegimes,1), ...
            repmat(dt,nRegimes,1), ...
            repmat(thin_step,nRegimes,1), ...
            retainedSamples_GOF, ...
            L1_distance_GOF, ...
            'VariableNames',{ ...
                'Regime', ...
                'r', ...
                'SimulationLength', ...
                'NumberOfTrajectories', ...
                'DiscardedTransient', ...
                'TimeStep', ...
                'ThinningStep', ...
                'RetainedObservations', ...
                'L1_Distance'});
    
        fprintf('\nFig. 3(a,b) quantitative goodness-of-fit results:\n');
        disp(Figure3ab_L1_results);

end


%% Fig. 4a,b: Large-noise limits
if RUN.c_d
    FS = STYLE.tickFontSize;
    LW = STYLE.plotLineWidth;

    % Parameters
    x_max = CD.xmax;

    % Base xl,xh
    xl0 = CD.xl;
    xh0 = CD.xh;

    beta = CD.beta;

    % Base diffusion coefficients
    b2_0 = CD.b2;
    b1_0 = CD.b1;
    b0_0 = CD.b0;

    % x grid
    gridN = CD.nX;
    xgrid = linspace(0, x_max, gridN).';

    alpha_over_beta_fixed = CD.alphaOverBeta;
    alpha_fixed = beta * alpha_over_beta_fixed;

    % Figure
    if CD.runSpecifiedKappa
        alpha_kappa = alpha_fixed;

        % Specify kappas 
        kappa_list = CD.kappaValues;
        colors3 = CD.curveColors;

        % Define quadratic base shape g0(x)=b2 x^2 + b1 x + b0 (>0)
        b0_g0 = b0_0;

        g0_grid = large_noise_limit__g_quad_positive(xgrid, b2_0, b1_0, b0_g0);
        rho_limit_quad = 1./(g0_grid.^2);
        rho_limit_quad = rho_limit_quad / trapz(xgrid, rho_limit_quad);

        rho_limit_const = ones(size(xgrid)) / x_max;

        % Compute densities for specified kappas
        nk = numel(kappa_list);
        rho_const_list = zeros(nk, numel(xgrid));
        rho_quad_list  = zeros(nk, numel(xgrid));

        for i = 1:nk
            kappa = kappa_list(i);

            % Case 1: constant diffusion g(x)=kappa
            rho_const_list(i,:) = large_noise_limit__stationary_density_numeric_constG( ...
                xgrid, alpha_kappa, beta, xl0, xh0, kappa, b0_0).';

            % Case 2: quadratic diffusion g(x)=kappa*g0(x)
            rho_quad_list(i,:) = large_noise_limit__stationary_density_numeric_scaledG0( ...
                xgrid, alpha_kappa, beta, xl0, xh0, b2_0, b1_0, b0_g0, kappa ).';
        end

        fig = publication_figure(STYLE.figureSizeIn);

        % Constant diffusion
        hold on; box on;
        for i = 1:nk
            plot(xgrid, rho_const_list(i,:), 'LineWidth',LW, 'Color', colors3(i,:));
        end
        plot(xgrid, rho_limit_const, 'k--', 'LineWidth',LW);
        xlabel('$x$','Interpreter','latex'); ylabel('$\rho^*(x)$','Interpreter','latex');

        leg = cell(1,nk+1);
        for i = 1:nk
            leg{i} = sprintf('$\\kappa=%.3g$', kappa_list(i));
        end
        leg{nk+1} = '$\kappa \to +\infty$';
        legend(leg,'Interpreter','latex','Location','best');
        style_publication_axes(gca, STYLE); xlim([0 x_max]);
        set(legend, 'FontSize', STYLE.annotationFontSize);

        % State-dependent diffusion
        fig = publication_figure(STYLE.figureSizeIn);
        hold on; box on;
        for i = 1:nk
            plot(xgrid, rho_quad_list(i,:), 'LineWidth',LW, 'Color', colors3(i,:));
        end
        plot(xgrid, rho_limit_quad, 'k--', 'LineWidth',LW);
        xlabel('$x$','Interpreter','latex'); ylabel('$\rho^*(x)$','Interpreter','latex');

        legend(leg,'Interpreter','latex','Location','best');
        style_publication_axes(gca, STYLE); xlim([0 x_max]);
        set(legend, 'FontSize', STYLE.annotationFontSize);
    end
end


%% Fig. 5a,b: Peak-configuration regions
if RUN.e_f
    FS = STYLE.tickFontSize;
    LW = STYLE.plotLineWidth;
    LW_bdy = STYLE.plotLineWidth;

    % Parameters
    xmax = EF.xmax;

    beta = EF.beta;

    % State-dependent diffusion: g(x)=b2*x^2 + b1*x + b0
    b2_sd = EF.b2;
    b1_sd = EF.b1;
    b0_sd = EF.b0;

    % Geometry switch
    use_centered_well = EF.useCenteredWell;
    xl_fixed = EF.xlFixed;

    % d = x_h - x_l
    d_min = EF.dMin;
    d_max = xmax - EF.dMaxPadding;

    Nd = EF.nD;
    Nr = EF.nR;

    d_vals = linspace(d_min, d_max, Nd);

    % alpha/beta range
    r_SN_max = (d_max^3)/(3*sqrt(3));
    r_lim = EF.rLimitFactor * r_SN_max;

    r_vals = linspace(-r_lim, r_lim, Nr);

    % Constant diffusion
    x_tmp = linspace(0, xmax, EF.nDiffusionCheckPoints);
    g_state_tmp = b2_sd*x_tmp.^2 + b1_sd*x_tmp + b0_sd;

    if any(g_state_tmp <= 0)
        error('Nonlinear diffusion g(x) is not strictly positive on [0,xmax].');
    end

    sigma_const = mean(g_state_tmp);

    % Storage
    Class_nl    = zeros(Nd, Nr);   % State-dependent diffusion display classes
    Class_const = zeros(Nd, Nr);   % Constant diffusion display classes

    Nint_nl     = zeros(Nd, Nr);   % Number of interior peaks only
    Nint_const  = zeros(Nd, Nr);

    Bdy0_nl     = false(Nd, Nr);   % x=0 boundary peak
    BdyX_nl     = false(Nd, Nr);   % x=xmax boundary peak
    Bdy0_const  = false(Nd, Nr);
    BdyX_const  = false(Nd, Nr);

    Disc_nl     = nan(Nd, Nr);     % Discriminant of psi(x), state-dependent diffusion
    Disc_const  = nan(Nd, Nr);     % Discriminant of psi(x), constant diffusion

    Psi0_nl     = nan(Nd, Nr);     % psi(0), state-dependent diffusion
    PsiX_nl     = nan(Nd, Nr);     % psi(xmax), state-dependent diffusion
    Psi0_const  = nan(Nd, Nr);     % psi(0), constant diffusion
    PsiX_const  = nan(Nd, Nr);     % psi(xmax), constant diffusion

    % Main sweep
    for ii = 1:Nd
        d = d_vals(ii);

        % Geometry
        if use_centered_well
            xl = (xmax - d)/2;
            xh = xl + d;
        else
            xl = xl_fixed;
            xh = xl + d;
        end

        for jj = 1:Nr
            r = r_vals(jj);
            alpha = beta * r;

            % Drift term
            % f(x) = -( alpha + 2*beta*(x-xl)*(x-xh)*(2*x-xl-xh) )
            f_poly = peak_number_map__drift_poly(alpha, beta, xl, xh);

            % 1) State-dependent diffusion term
            g_poly  = [b2_sd, b1_sd, b0_sd];   % quadratic
            gp_poly = [2*b2_sd, b1_sd];        % linear

            ggprime_poly = conv(g_poly, gp_poly);  % cubic polynomial g(x)g'(x)
            psi_nl_poly  = f_poly - ggprime_poly;  % cubic polynomial psi(x)

            Disc_nl(ii,jj) = peak_number_map__cubic_discriminant(psi_nl_poly);
            Psi0_nl(ii,jj) = polyval(psi_nl_poly, 0);
            PsiX_nl(ii,jj) = polyval(psi_nl_poly, xmax);

            [Class_nl(ii,jj), Nint_nl(ii,jj), Bdy0_nl(ii,jj), BdyX_nl(ii,jj)] = ...
                peak_number_map__classify_stationary_peaks(psi_nl_poly, xmax);

            % 2) Constant diffusion term
            psi_const_poly = f_poly;

            Disc_const(ii,jj) = peak_number_map__cubic_discriminant(psi_const_poly);
            Psi0_const(ii,jj) = polyval(psi_const_poly, 0);
            PsiX_const(ii,jj) = polyval(psi_const_poly, xmax);

            [Class_const(ii,jj), Nint_const(ii,jj), Bdy0_const(ii,jj), BdyX_const(ii,jj)] = ...
                peak_number_map__classify_stationary_peaks(psi_const_poly, xmax);
        end
    end

    % class_map values:
    %   0 = other
    %   1 = one interior peak only
    %   2 = two interior peaks only
    %   3 = only one peak at x = 0
    %   4 = only one peak at x = xmax
    %   5 = one peak at x = 0 plus one interior peak
    %   6 = one peak at x = xmax plus one interior peak
    %   7 = two boundary peaks only
    %   8 = two boundary peaks plus one interior peak
    %   9 = two boundary peaks plus two interior peaks

    c_other     = [1.0000 1.0000 1.0000];   % white, other
    c_gray      = [202 200 201]/255;        % one interior
    c_orange    = [238 184 100]/255;        % two interiors
    
    c_x0        = [143 181 216]/255;        % only x=0
    c_xmax      = [201 116 104]/255;        % only xmax
    c_x0_int    = [0, 119, 188]/255;        % x=0 + interior
    c_xmax_int  = [0.55 0.10 0.10];         % xmax + interior
    
    c_two_bdy   = [0.00 0.00 0.00];         % black: two boundary peaks only
    c_two_bdy_1 = [0.50 0.00 0.50];         % purple: two boundaries + one interior
    c_two_bdy_2 = [0.00 0.60 0.20];         % green: two boundaries + two interiors
    
    myColors = [
        c_other;       % 0 other
        c_gray;        % 1 one interior
        c_orange;      % 2 two interiors
        c_x0;          % 3 only x=0
        c_xmax;        % 4 only xmax
        c_x0_int;      % 5 x=0 + interior
        c_xmax_int;    % 6 xmax + interior
        c_two_bdy;     % 7 two boundary peaks only
        c_two_bdy_1;   % 8 two boundaries + one interior
        c_two_bdy_2    % 9 two boundaries + two interiors
    ];
    
    class_ticks  = 0:9;
    class_labels = {
        'other', ...
        '1 int', ...
        '2 int', ...
        'x=0', ...
        'x_{max}', ...
        'x=0+int', ...
        'x_{max}+int', ...
        '2 bdy', ...
        '2 bdy+1 int', ...
        '2 bdy+2 int'
    };

    % Plot: state-dependent diffusion
    fig = publication_figure(STYLE.figureSizeIn);

    imagesc(r_vals, d_vals, Class_nl);
    set(gca,'YDir','normal');
    hold on;
    
    % Add hatching to selected regions
    % Class 5: deep-blue region, vertical hatching
    overlay_class_hatching(gca, r_vals, d_vals, Class_nl, ...
        5, 'vertical', 8, [0 0 0], (2/3)*LW_bdy);
    
    % Class 6: deep-red region, horizontal hatching
    overlay_class_hatching(gca, r_vals, d_vals, Class_nl, ...
        6, 'horizontal', 6, [0 0 0], (2/3)*LW_bdy);
    
    % Interior peak bifurcation boundary
    contour(r_vals, d_vals, Disc_nl, [0 0], 'k--', 'LineWidth', LW);

    % Boundary peak theoretical boundaries
    contour(r_vals, d_vals, Psi0_nl, [0 0], 'w--', 'LineWidth', LW_bdy);
    contour(r_vals, d_vals, PsiX_nl, [0 0], 'w--', 'LineWidth', LW_bdy);

    hold off;
    axis tight;
    box on;
    grid on;

    xlabel('$\alpha/\beta$','Interpreter','latex');
    ylabel('$d=x_h-x_l$','Interpreter','latex');

    colormap(gca, myColors);
    caxis([-0.5 9.5]);

    style_publication_axes(gca, STYLE);

    % Plot: constant diffusion
    fig = publication_figure(STYLE.figureSizeIn);

    imagesc(r_vals, d_vals, Class_const);
    set(gca,'YDir','normal');
    hold on;
    
    % Add hatching to selected regions
    % Class 5: deep-blue region, vertical hatching
    overlay_class_hatching(gca, r_vals, d_vals, Class_const, ...
        5, 'vertical', 8, [0 0 0], (2/3)*LW_bdy);
    
    % Class 6: deep-red region, horizontal hatching
    overlay_class_hatching(gca, r_vals, d_vals, Class_const, ...
        6, 'horizontal', 6, [0 0 0], (2/3)*LW_bdy);
    
    % Interior peak bifurcation boundary
    d_curve = linspace(d_min, d_max, EF.nBoundaryPoints);
    
    r_SN_pos =  d_curve.^3 ./ (3*sqrt(3));
    r_SN_neg = -d_curve.^3 ./ (3*sqrt(3));
    
    plot(r_SN_pos, d_curve, 'k--', 'LineWidth', LW);
    plot(r_SN_neg, d_curve, 'k--', 'LineWidth', LW);

    % Boundary peak theoretical boundaries
    contour(r_vals, d_vals, Psi0_const, [0 0], 'w--', 'LineWidth', LW_bdy);
    contour(r_vals, d_vals, PsiX_const, [0 0], 'w--', 'LineWidth', LW_bdy);

    hold off;
    axis tight;
    box on;
    grid on;

    xlabel('$\alpha/\beta$','Interpreter','latex');
    ylabel('$d=x_h-x_l$','Interpreter','latex');

    colormap(gca, myColors);
    caxis([-0.5 9.5]);

    style_publication_axes(gca, STYLE);
end


%% Fig. 5c,d: Stationary-density heat maps
if RUN.g_h
    FS = STYLE.tickFontSize;
    xmax = GH.xmax;

    beta = GH.beta;

    b2_sd = GH.b2;
    b1_sd = GH.b1;
    b0_sd = GH.b0;

    use_centered_well = GH.useCenteredWell;
    xl_fixed = GH.xlFixed;
    d0 = GH.d;

    % Parameter sweep
    rb = d0^3 / (3*sqrt(3));

    Nr = GH.nR;
    r_list = linspace(GH.rRangeInSN(1)*rb, GH.rRangeInSN(2)*rb, Nr);

    % Simulation parameters
    dt        = GH.dt;
    T         = GH.simulationTime;
    T_burn    = GH.burninTime;
    Ntraj     = GH.nTrajectories;
    thin_step = GH.thinStep;
    x0        = GH.initialState;

    Nbins     = GH.nBins;
    edges = linspace(0, xmax, Nbins+1);
    x_centers = 0.5*(edges(1:end-1) + edges(2:end));

    smooth_bins = GH.smoothBins;

    % Theoretical peaks
    Nx_peak_grid = GH.nPeakGrid;
    xg_peak = linspace(0, xmax, Nx_peak_grid);

    rng(GH.randomSeed);

    RHO_sim_state = zeros(Nr, Nbins);
    RHO_sim_const = zeros(Nr, Nbins);

    xpk_state_all = cell(Nr,1);
    xpk_const_all = cell(Nr,1);


    %% Monte Carlo simulations
    for ir = 1:Nr
        r = r_list(ir);

        fprintf('Simulating %d / %d, alpha/beta = %.2f\n', ir, Nr, r);

        C_state = simulation_density_map__make_case('state diffusion', ...
            r, d0, beta, xmax, use_centered_well, xl_fixed, ...
            'state', b0_sd, b1_sd, b2_sd);

        C_const = simulation_density_map__make_case('constant diffusion', ...
            r, d0, beta, xmax, use_centered_well, xl_fixed, ...
            'constant', b0_sd, b1_sd, b2_sd);

        rho_state = simulation_density_map__simulate_density_reflected_rsde(C_state.f, C_state.g, ...
            x0, xmax, dt, T, T_burn, Ntraj, thin_step, edges);

        rho_const = simulation_density_map__simulate_density_reflected_rsde(C_const.f, C_const.g, ...
            x0, xmax, dt, T, T_burn, Ntraj, thin_step, edges);

        rho_state = simulation_density_map__smooth_density_1d(x_centers, rho_state, smooth_bins);
        rho_const = simulation_density_map__smooth_density_1d(x_centers, rho_const, smooth_bins);

        RHO_sim_state(ir,:) = rho_state;
        RHO_sim_const(ir,:) = rho_const;

        rho_state_theory = simulation_density_map__stationary_density_general(xg_peak, C_state.f, C_state.g);
        rho_const_theory = simulation_density_map__stationary_density_general(xg_peak, C_const.f, C_const.g);

        [xpk_state, ~] = simulation_density_map__theoretical_peaks(xg_peak, rho_state_theory);
        [xpk_const, ~] = simulation_density_map__theoretical_peaks(xg_peak, rho_const_theory);

        xpk_state_all{ir} = xpk_state;
        xpk_const_all{ir} = xpk_const;
    end

    % Difference map
    D_sim = RHO_sim_state - RHO_sim_const;

    rho_color_max = max([RHO_sim_state(:); RHO_sim_const(:)]);
    diff_abs_max  = max(abs(D_sim(:)));

    idx_peak_plot = 1:Nr;

    % colormaps
    c_none = [1.0000 1.0000 1.0000];       % white
    c_x0   = [143 181 216]/255;            % light blue
    c_deepblue = [0, 119, 188]/255;        % deep blue
    c_xmax = [201 116 104]/255;            % light red
    c_deepred = [0.55 0.10 0.10];          % deep red

    cmap_rho = simulation_density_map__interp_cmap_code1([ ...
        c_none;
        c_x0;
        c_deepblue;
        c_xmax;
        c_deepred], 256);

    cmap_diff = simulation_density_map__interp_cmap_code1([ ...
        c_deepblue;
        c_x0;
        c_none;
        c_xmax;
        c_deepred], 256);

    % Plot
    % Figure: constant diffusion numerical density
    fig = publication_figure(STYLE.figureSizeIn);

    imagesc(x_centers, r_list, RHO_sim_const);
    set(gca,'YDir','normal');
    hold on;

    % peak locations
    for ii = 1:numel(idx_peak_plot)
        ir = idx_peak_plot(ii);
        xpk = xpk_const_all{ir};
        if ~isempty(xpk)
            plot(xpk, r_list(ir)*ones(size(xpk)), 'w.', 'MarkerSize', 18);
        end
    end

    xlabel('$x$','Interpreter','latex');
    ylabel('$\alpha/\beta$','Interpreter','latex');

    style_publication_axes(gca, STYLE);
    xlim([0 xmax]);

    colormap(gca, cmap_rho);
    caxis([0 rho_color_max]);

    cb1 = colorbar;
    cb1.Title.String = '$\rho$';
    cb1.Title.Interpreter = 'latex';
    cb1.Title.FontSize = STYLE.annotationFontSize;
    cb1.FontSize = STYLE.tickFontSize;
    cb1.FontWeight = 'bold';
    cb1.LineWidth = STYLE.axisLineWidth;

    grid on;
    box on;

    % Figure: state-dependent diffusion numerical density
    fig = publication_figure(STYLE.figureSizeIn);

    imagesc(x_centers, r_list, RHO_sim_state);
    set(gca,'YDir','normal');
    hold on;

    % peak locations
    for ii = 1:numel(idx_peak_plot)
        ir = idx_peak_plot(ii);
        xpk = xpk_state_all{ir};
        if ~isempty(xpk)
            plot(xpk, r_list(ir)*ones(size(xpk)), 'w.', 'MarkerSize', 18);
        end
    end

    xlabel('$x$','Interpreter','latex');
    ylabel('$\alpha/\beta$','Interpreter','latex');

    style_publication_axes(gca, STYLE);
    xlim([0 xmax]);

    colormap(gca, cmap_rho);
    caxis([0 rho_color_max]);

    cb2 = colorbar;
    cb2.Title.String = '$\rho$';
    cb2.Title.Interpreter = 'latex';
    cb2.Title.FontSize = STYLE.annotationFontSize;
    cb2.FontSize = STYLE.tickFontSize;
    cb2.FontWeight = 'bold';
    cb2.LineWidth = STYLE.axisLineWidth;

    grid on;
    box on;
end


%% Fig. 6a-d: Stationary-density extrema versus noise scale
if RUN.i_l
    x_l = IL.xl;
    x_h = IL.xh;
    x_max = IL.xmax;
    beta = IL.beta;
    b2 = IL.b2;
    b1 = IL.b1;
    b0 = IL.b0;

    x_check = linspace(0, x_max, IL.nPositivityCheckPoints);
    g0_check = b2*x_check.^2 + b1*x_check + b0;
    assert(min(g0_check) > 0, ...
        'The Figure 3i-l diffusion amplitude must be positive on [0,xmax].');

    kappa_values = linspace(IL.kappaMin, IL.kappaMax, IL.nKappa);
    d = x_h - x_l;
    r_SN = d^3/(3*sqrt(3));
    r_values = r_SN * IL.rRatios;
    alpha_values = beta * r_values;

    fprintf('Figure 3i-l extrema sweep\n');
    fprintf('  beta = %.8g, r_SN = %.8g\n', beta, r_SN);

    for caseIndex = 1:numel(alpha_values)
        alpha = alpha_values(caseIndex);

        kappa_maxima = [];
        x_maxima = [];
        kappa_minima = [];
        x_minima = [];
        kappa_deg = [];
        x_deg = [];
        delta_pk_values = zeros(size(kappa_values));

        for i = 1:numel(kappa_values)
            kappa = kappa_values(i);
            delta_pk_values(i) = figure3il__peak_discriminant( ...
                kappa, alpha, beta, x_l, x_h, b0, b1, b2);

            [roots_inside, root_type] = figure3il__stationary_extrema( ...
                kappa, alpha, beta, x_l, x_h, b0, b1, b2, x_max);

            for j = 1:numel(roots_inside)
                if root_type(j) == 1
                    kappa_maxima(end+1) = kappa;
                    x_maxima(end+1) = roots_inside(j);
                elseif root_type(j) == -1
                    kappa_minima(end+1) = kappa;
                    x_minima(end+1) = roots_inside(j);
                else
                    kappa_deg(end+1) = kappa;
                    x_deg(end+1) = roots_inside(j);
                end
            end
        end

        crossing_indices = find( ...
            sign(delta_pk_values(1:end-1)) .* ...
            sign(delta_pk_values(2:end)) < 0);

        for crossingIndex = crossing_indices
            kappa_bracket = kappa_values([crossingIndex, crossingIndex+1]);
            kappa_fold = fzero( ...
                @(kappa) figure3il__peak_discriminant( ...
                    kappa, alpha, beta, x_l, x_h, b0, b1, b2), ...
                kappa_bracket);
            x_fold = figure3il__repeated_interior_root( ...
                kappa_fold, alpha, beta, x_l, x_h, ...
                b0, b1, b2, x_max);

            if isfinite(x_fold)
                kappa_deg(end+1) = kappa_fold;
                x_deg(end+1) = x_fold;
            end
        end

        fig = publication_figure(STYLE.figureSizeIn);
        ax = axes(fig);
        hold(ax, 'on');
        box(ax, 'on');
        grid(ax, 'on');

        hMax = plot(ax, nan, nan, '.', ...
            'Color', IL.maximumColor, ...
            'MarkerSize', IL.maximumMarkerSize, ...
            'DisplayName', 'Local maxima');
        hMin = plot(ax, nan, nan, 'o', ...
            'Color', IL.minimumColor, ...
            'MarkerSize', IL.minimumMarkerSize, ...
            'LineWidth', STYLE.plotLineWidth, ...
            'DisplayName', 'Local minima');
        hDeg = plot(ax, nan, nan, 'kx', ...
            'MarkerSize', IL.foldMarkerSize, ...
            'LineWidth', STYLE.plotLineWidth, ...
            'DisplayName', 'Fold point ($\Delta_{pk}=0$)');

        plot(ax, kappa_maxima, x_maxima, '.', ...
            'Color', IL.maximumColor, ...
            'MarkerSize', IL.maximumMarkerSize, ...
            'HandleVisibility', 'off');
        plot(ax, kappa_minima, x_minima, 'o', ...
            'Color', IL.minimumColor, ...
            'MarkerSize', IL.minimumMarkerSize, ...
            'LineWidth', STYLE.plotLineWidth, ...
            'HandleVisibility', 'off');
        plot(ax, kappa_deg, x_deg, 'kx', ...
            'MarkerSize', IL.foldMarkerSize, ...
            'LineWidth', STYLE.plotLineWidth, ...
            'HandleVisibility', 'off');

        for kappaReference = IL.kappaReferenceLines
            xline(ax, kappaReference, ':', ...
                'Color', [0.45 0.45 0.45], ...
                'HandleVisibility', 'off');
        end

        xlabel(ax, '$\kappa$', 'Interpreter', 'latex');
        ylabel(ax, '$x_{pk}$', 'Interpreter', 'latex');
        xlim(ax, [IL.kappaMin, IL.kappaMax]);
        ylim(ax, [0, x_max]);
        style_publication_axes(ax, STYLE);

        if caseIndex == numel(alpha_values)
            legend(ax, [hMax, hMin, hDeg], ...
                'Interpreter', 'latex', ...
                'Location', 'best', ...
                'FontSize', STYLE.annotationFontSize);
        end
    end
end





function C = peak_equilibrium_comparison__make_case(name, r, d, beta, xmax, use_centered_well, xl_fixed, ...
    diff_type, b0_sd, b1_sd, b2_sd)

    if use_centered_well
        xl = (xmax - d)/2;
        xh = xl + d;
    else
        xl = xl_fixed;
        xh = xl + d;
        % xl = 30;
        % xh = 70;
    end

    alpha = beta * r;

    % drift: f(x) = - dV/dx
    f = @(x) -( alpha + 2*beta*(x-xl).*(x-xh).*(2*x-xl-xh) );

    % state-dependent diffusion
    g_state = @(x) b2_sd*x.^2 + b1_sd*x + b0_sd;

    switch lower(diff_type)
        case 'state'
            g = @(x) g_state(x);

        case 'constant'
            x_tmp = linspace(0, xmax, 2000);
            sigma = mean(g_state(x_tmp));
            g = @(x) sigma + 0*x;

        otherwise
            error('Unknown diffusion type. Use ''state'' or ''constant''.');
    end

    % positivity check
    x_chk = linspace(0, xmax, 2000);
    gx_chk = g(x_chk);
    if any(gx_chk <= 0)
        error('Diffusion g(x) must stay strictly positive on [0,xmax].');
    end

    C = struct();
    C.name  = name;
    C.r     = r;
    C.d     = d;
    C.beta  = beta;
    C.alpha = alpha;
    C.xl    = xl;
    C.xh    = xh;
    C.f     = f;
    C.g     = g;
end


function rho = peak_equilibrium_comparison__stationary_density_general(x, f, g)
    gx = g(x);
    q  = 2*f(x) ./ (gx.^2);

    I = cumtrapz(x, q);

    log_rho = I - 2*log(gx);

    log_rho = log_rho - max(log_rho);

    rho = exp(log_rho);

    rho(~isfinite(rho)) = 0;
    rho(rho < 0) = 0;

    Z = trapz(x, rho);
    if Z <= 0
        error('Failed to normalize theoretical rho*(x).');
    end

    rho = rho / Z;
end


function samples = peak_equilibrium_comparison__simulate_reflected_rsde(f, g, x0, xmax, ...
    dt, T, T_burn, Ntraj, thin_step)

    Nt      = round(T/dt);
    Nt_burn = round(T_burn/dt);

    samples = [];

    for n = 1:Ntraj
        x = x0;

        for k = 1:Nt
            xi = randn;
            x_new = x + f(x)*dt + g(x)*sqrt(dt)*xi;

            % reflected step into [0, xmax]
            x = peak_equilibrium_comparison__reflect_to_interval(x_new, 0, xmax);

            if k > Nt_burn && mod(k, thin_step) == 0
                samples(end+1,1) = x;
            end
        end
    end
end


function x = peak_equilibrium_comparison__reflect_to_interval(x, a, b)
    while x < a || x > b
        if x < a
            x = 2*a - x;
        elseif x > b
            x = 2*b - x;
        end
    end

    x = min(max(x, a), b);
end


function xeq = peak_equilibrium_comparison__deterministic_equilibria(alpha, beta, xl, xh, xmax)
    % deterministic equilibria solve f(x)=0
    a3 = -4*beta;
    a2 =  6*beta*(xl + xh);
    a1 = -2*beta*(xl^2 + xh^2 + 4*xl*xh);
    a0 = -alpha + 2*beta*(xl^2*xh + xh^2*xl);

    rr = roots([a3 a2 a1 a0]);
    rr = rr(abs(imag(rr)) < 1e-10);
    rr = real(rr);

    rr = sort(rr);
    rr = rr(rr >= 0 & rr <= xmax);

    if isempty(rr)
        xeq = [];
        return;
    end

    tol = 1e-7;
    keep = true(size(rr));
    for k = 2:numel(rr)
        if abs(rr(k)-rr(k-1)) < tol
            keep(k) = false;
        end
    end
    xeq = rr(keep);
end


function [xpk, ypk] = peak_equilibrium_comparison__theoretical_peaks(x, rho)
    % local maxima finder
    xpk = [];
    ypk = [];

    if numel(x) < 3
        return;
    end

    idx = find(rho(2:end-1) >= rho(1:end-2) & ...
               rho(2:end-1) >= rho(3:end)) + 1;

    if isempty(idx)
        return;
    end

    thr = 0.00 * max(rho);
    idx = idx(rho(idx) >= thr);

    xpk = x(idx);
    ypk = rho(idx);
end


function xeq_stable = peak_equilibrium_comparison__keep_stable_equilibria(xeq, f, xmax)

    if isempty(xeq)
        xeq_stable = [];
        return;
    end

    h = 1e-5 * xmax;
    is_stable = false(size(xeq));

    for k = 1:numel(xeq)
        x0 = xeq(k);

        xp = min(x0 + h, xmax);
        xm = max(x0 - h, 0);

        df = (f(xp) - f(xm)) / (xp - xm);

        is_stable(k) = df < 0;
    end

    xeq_stable = xeq(is_stable);
end


function rho = large_noise_limit__stationary_density_numeric_constG(xgrid, alpha, beta, xl, xh, kappa, b0_const)
    % constant diffusion g(x)=kappa*b0
    fx = large_noise_limit__f_from_potential(xgrid, alpha, beta, xl, xh);
    g2 = (kappa^2) * (b0_const^2) * ones(size(xgrid));
    integrand = 2*fx ./ g2;
    I = cumtrapz(xgrid, integrand);
    I = I - max(I);
    unnorm = (1./g2) .* exp(I);
    Z = trapz(xgrid, unnorm);
    rho = unnorm / Z;
end


function rho = large_noise_limit__stationary_density_numeric_scaledG0(xgrid, alpha, beta, xl, xh, b2, b1, b0, kappa)
    % quadratic diffusion g(x)=kappa*g0(x), where g0(x)=b2x^2+b1x+b0 (must be >0)
    fx = large_noise_limit__f_from_potential(xgrid, alpha, beta, xl, xh);
    g0 = large_noise_limit__g_quad_positive(xgrid, b2, b1, b0);
    g2 = (kappa^2) * (g0.^2);
    integrand = 2*fx ./ g2;
    I = cumtrapz(xgrid, integrand);
    I = I - max(I);
    unnorm = (1./g2) .* exp(I);
    Z = trapz(xgrid, unnorm);
    rho = unnorm / Z;
end


function fx = large_noise_limit__f_from_potential(x, alpha, beta, xl, xh)
    A = x - xl;
    B = x - xh;
    d_term = 2 .* A .* B .* (2.*x - xl - xh);
    Vprime = alpha + beta .* d_term;
    fx = -Vprime;
end


function gx = large_noise_limit__g_quad_positive(x, b2, b1, b0)
    gx = b2.*x.^2 + b1.*x + b0;
    if any(gx <= 0)
        error('g(x) <= 0 detected. Fix (b2,b1,b0) construction instead of clamping.');
    end
end


function f_poly = peak_number_map__drift_poly(alpha, beta, xl, xh)
    p1 = conv([1, -xl], [1, -xh]);      % (x-xl)(x-xh)
    p2 = [2, -(xl + xh)];               % (2x-xl-xh)
    p3 = conv(p1, p2);                  % cubic

    f_poly = -2*beta * p3;
    f_poly(end) = f_poly(end) - alpha;
end


function [class_id, nInt, isLeftBdyPeak, isRightBdyPeak] = peak_number_map__classify_stationary_peaks(psi_poly, xmax)

    tol_bdy = 1e-10;

    nInt = peak_number_map__count_internal_peaks_from_psi(psi_poly, xmax);

    psi0 = polyval(psi_poly, 0);
    psiX = polyval(psi_poly, xmax);

    isLeftBdyPeak  = (psi0 <=  tol_bdy);
    isRightBdyPeak = (psiX >= -tol_bdy);

    % Default: other
    class_id = 0;

    if ~isLeftBdyPeak && ~isRightBdyPeak
        if nInt == 1
            class_id = 1;
        elseif nInt == 2
            class_id = 2;
        end

    elseif isLeftBdyPeak && ~isRightBdyPeak
        if nInt == 0
            class_id = 3;
        elseif nInt == 1
            class_id = 5;
        end

    elseif ~isLeftBdyPeak && isRightBdyPeak
        if nInt == 0
            class_id = 4;
        elseif nInt == 1
            class_id = 6;
        end

    elseif isLeftBdyPeak && isRightBdyPeak
        if nInt == 0
            class_id = 7;   % two boundary peaks only
        elseif nInt == 1
            class_id = 8;   % two boundaries + one interior peak
        elseif nInt == 2
            class_id = 9;   % two boundaries + two interior peaks
        end
    end
end


function npeak = peak_number_map__count_internal_peaks_from_psi(psi_poly, xmax)
    tol_imag = 1e-9;
    tol_x    = 1e-7;
    tol_der  = 1e-8;

    % roots of psi
    rr = roots(psi_poly);

    % keep real roots only
    rr = rr(abs(imag(rr)) < tol_imag);
    rr = real(rr);

    % keep only internal roots
    rr = rr(rr > tol_x & rr < xmax - tol_x);

    if isempty(rr)
        npeak = 0;
        return;
    end

    % cluster nearly repeated roots
    rr = peak_number_map__unique_tol(sort(rr), 1e-6);

    % derivative of psi
    dpsi_poly = polyder(psi_poly);

    % count only maxima: psi'(x*) < 0
    isPeak = false(size(rr));
    for k = 1:length(rr)
        val = polyval(dpsi_poly, rr(k));
        if val < -tol_der
            isPeak(k) = true;
        end
    end

    npeak = sum(isPeak);
end


function x_unique = peak_number_map__unique_tol(x, tol)
    if isempty(x)
        x_unique = x;
        return;
    end

    x = sort(x(:));
    x_unique = x(1);

    for k = 2:length(x)
        if abs(x(k) - x_unique(end)) > tol
            x_unique(end+1,1) = x(k);
        else
            x_unique(end) = 0.5*(x_unique(end) + x(k));
        end
    end
end


function Delta = peak_number_map__cubic_discriminant(p)
% Discriminant of cubic polynomial
    a = p(1);
    b = p(2);
    c = p(3);
    d = p(4);

    Delta = 18*a*b*c*d ...
          - 4*b^3*d ...
          + b^2*c^2 ...
          - 4*a*c^3 ...
          - 27*a^2*d^2;
end


function C = simulation_density_map__make_case(name, r, d, beta, xmax, use_centered_well, xl_fixed, ...
    diff_type, b0_sd, b1_sd, b2_sd)

    if use_centered_well
        xl = (xmax - d)/2;
        xh = xl + d;
    else
        xl = xl_fixed;
        xh = xl + d;
    end

    alpha = beta * r;

    % drift: f(x) = -dV/dx
    f = @(x) -( alpha + 2*beta*(x-xl).*(x-xh).*(2*x-xl-xh) );

    % nonlinear state-dependent diffusion
    g_state = @(x) b2_sd*x.^2 + b1_sd*x + b0_sd;

    switch lower(diff_type)
        case 'state'
            g = @(x) g_state(x);

        case 'constant'
            x_tmp = linspace(0, xmax, 2000);
            sigma = mean(g_state(x_tmp));
            g = @(x) sigma + 0*x;

        otherwise
            error('Unknown diffusion type. Use ''state'' or ''constant''.');
    end

    % positivity check
    x_chk = linspace(0, xmax, 2000);
    gx_chk = g(x_chk);

    if any(gx_chk <= 0)
        error('Diffusion g(x) must stay strictly positive on [0,xmax].');
    end

    C = struct();
    C.name  = name;
    C.r     = r;
    C.d     = d;
    C.beta  = beta;
    C.alpha = alpha;
    C.xl    = xl;
    C.xh    = xh;
    C.f     = f;
    C.g     = g;
end


function rho = simulation_density_map__simulate_density_reflected_rsde(f, g, x0, xmax, ...
    dt, T, T_burn, Ntraj, thin_step, edges)

    Nt      = round(T/dt);
    Nt_burn = round(T_burn/dt);

    x = x0 * ones(Ntraj,1);

    counts = zeros(1, numel(edges)-1);
    n_saved = 0;

    sqrt_dt = sqrt(dt);

    for k = 1:Nt
        xi = randn(Ntraj,1);

        x_new = x + f(x)*dt + g(x)*sqrt_dt.*xi;

        x = simulation_density_map__reflect_to_interval_vectorized(x_new, 0, xmax);

        if k > Nt_burn && mod(k, thin_step) == 0
            counts = counts + histcounts(x, edges);
            n_saved = n_saved + Ntraj;
        end
    end

    bin_width = edges(2) - edges(1);
    rho = counts / (n_saved * bin_width);

    centers = 0.5*(edges(1:end-1) + edges(2:end));
    Z = trapz(centers, rho);

    if Z > 0
        rho = rho / Z;
    end
end


function x_ref = simulation_density_map__reflect_to_interval_vectorized(x, a, b)
    L = b - a;
    y = mod(x - a, 2*L);
    x_ref = a + min(y, 2*L - y);

    x_ref = min(max(x_ref, a), b);
end


function rho_s = simulation_density_map__smooth_density_1d(x, rho, smooth_bins)
    if smooth_bins <= 1
        rho_s = rho;
        return;
    end

    kernel = ones(1, smooth_bins) / smooth_bins;
    rho_s = conv(rho, kernel, 'same');

    Z = trapz(x, rho_s);

    if Z > 0
        rho_s = rho_s / Z;
    end
end


function rho = simulation_density_map__stationary_density_general(x, f, g)
    gx = g(x);
    q  = 2*f(x) ./ (gx.^2);

    I = cumtrapz(x, q);

    log_rho = I - 2*log(gx);
    log_rho = log_rho - max(log_rho);

    rho = exp(log_rho);

    rho(~isfinite(rho)) = 0;
    rho(rho < 0) = 0;

    Z = trapz(x, rho);

    if Z <= 0
        error('Failed to normalize theoretical rho*(x).');
    end

    rho = rho / Z;
end


function [xpk, ypk] = simulation_density_map__theoretical_peaks(x, rho)
    xpk = [];
    ypk = [];

    if numel(x) < 3
        return;
    end

    idx = find(rho(2:end-1) >= rho(1:end-2) & ...
               rho(2:end-1) >= rho(3:end)) + 1;

    if isempty(idx)
        return;
    end

    thr = 0.00 * max(rho);
    idx = idx(rho(idx) >= thr);

    xpk = x(idx);
    ypk = rho(idx);
end


function cmap = simulation_density_map__interp_cmap_code1(anchorColors, m)
    if nargin < 2
        m = 256;
    end

    n = size(anchorColors, 1);

    x  = linspace(1, m, n);
    xi = 1:m;

    cmap = [
        interp1(x, anchorColors(:,1), xi)', ...
        interp1(x, anchorColors(:,2), xi)', ...
        interp1(x, anchorColors(:,3), xi)'
    ];

    cmap = min(max(cmap, 0), 1);
end


function fig = publication_figure(figureSizeIn)
    fig = figure('Color', 'w', 'Units', 'inches', ...
        'Position', [1 1 figureSizeIn], 'PaperPositionMode', 'auto');
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 figureSizeIn];
    fig.PaperSize = figureSizeIn;
end


function style_publication_axes(ax, style)
    ax.FontSize = style.tickFontSize;
    ax.FontWeight = 'bold';
    ax.LineWidth = style.axisLineWidth;
    ax.TickDir = 'out';
    ax.XLabel.FontSize = style.labelFontSize;
    ax.YLabel.FontSize = style.labelFontSize;
    ax.XLabel.FontWeight = 'bold';
    ax.YLabel.FontWeight = 'bold';
    ax.Title.FontSize = style.annotationFontSize;
end


function [rootsInside, rootType] = figure3il__stationary_extrema( ...
    kappa, alpha, beta, xl, xh, b0, b1, b2, xmax)

    coefficients = figure3il__peak_polynomial_coefficients( ...
        kappa, alpha, beta, xl, xh, b0, b1, b2);
    candidateRoots = roots(coefficients);
    realTolerance = 1e-8;
    candidateRoots = real(candidateRoots( ...
        abs(imag(candidateRoots)) < realTolerance));
    rootsInside = sort(candidateRoots( ...
        candidateRoots > 0 & candidateRoots < xmax));
    rootType = zeros(size(rootsInside));

    for j = 1:numel(rootsInside)
        xpk = rootsInside(j);
        fPrime = -4*beta*( ...
            3*(xpk-(xl+xh)/2)^2-(xh-xl)^2/4);
        gValue = kappa*(b2*xpk^2 + b1*xpk + b0);
        gPrime = kappa*(2*b2*xpk + b1);
        gSecond = 2*kappa*b2;
        psiPrime = fPrime - gPrime^2 - gValue*gSecond;

        if psiPrime < -1e-7
            rootType(j) = 1;
        elseif psiPrime > 1e-7
            rootType(j) = -1;
        end
    end
end


function coefficients = figure3il__peak_polynomial_coefficients( ...
    kappa, alpha, beta, xl, xh, b0, b1, b2)

    B0 = kappa*b0;
    B1 = kappa*b1;
    B2 = kappa*b2;

    eta3 = -2*(2*beta + B2^2);
    eta2 = 6*beta*(xl+xh) - 3*B1*B2;
    eta1 = -2*beta*(xl^2 + xh^2 + 4*xl*xh) ...
        - (B1^2 + 2*B0*B2);
    eta0 = -alpha + 2*beta*(xl^2*xh + xh^2*xl) - B0*B1;
    coefficients = [eta3, eta2, eta1, eta0];
end


function deltaPk = figure3il__peak_discriminant( ...
    kappa, alpha, beta, xl, xh, b0, b1, b2)

    coefficients = figure3il__peak_polynomial_coefficients( ...
        kappa, alpha, beta, xl, xh, b0, b1, b2);
    a = coefficients(1);
    b = coefficients(2);
    c = coefficients(3);
    d = coefficients(4);
    deltaPk = 18*a*b*c*d - 4*b^3*d + b^2*c^2 ...
        - 4*a*c^3 - 27*a^2*d^2;
end


function xFold = figure3il__repeated_interior_root( ...
    kappa, alpha, beta, xl, xh, b0, b1, b2, xmax)

    coefficients = figure3il__peak_polynomial_coefficients( ...
        kappa, alpha, beta, xl, xh, b0, b1, b2);
    derivativeCoefficients = [ ...
        3*coefficients(1), 2*coefficients(2), coefficients(3)];
    candidates = roots(derivativeCoefficients);
    candidates = real(candidates(abs(imag(candidates)) < 1e-7));
    candidates = candidates(candidates > 0 & candidates < xmax);

    if isempty(candidates)
        xFold = nan;
        return;
    end

    residuals = abs(polyval(coefficients, candidates));
    [~, bestIndex] = min(residuals);
    xFold = candidates(bestIndex);
end


function overlay_class_hatching(ax, xValues, yValues, classMap, ...
        classValue, direction, stride, lineColor, lineWidth)

    mask = (classMap == classValue);

    switch direction
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
    % Find every continuous true interval
    transitions = diff([false; maskVector(:); false]);

    runStarts = find(transitions == 1);
    runEnds   = find(transitions == -1) - 1;

    runs = [runStarts, runEnds];
end
