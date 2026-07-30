%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project : How noise shapes decision-making dynamics in fish schools: 
% Insights from a stochastic model with state-dependent diffusion
% Author  : Deze Liu, Daniel Burbano Lombana
% Lab     : The Swarm Intelligence Lab
% Date    : 08/01/2026
% Description :
% This script generates the results for Fig. 9.
%           nominal local discriminant sensitivity and Monte
%           Carlo peak-change probability for the Bright and Dark nominal sets.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;
clear;
close all;

% Initial settings
RUN.bright = true;
RUN.dark   = true;

% Bright empirical condition
CONDITIONS(1).name = "Bright";
CONDITIONS(1).a3Source = -5.244444e-5;
CONDITIONS(1).a2Source =  0.004555556;
CONDITIONS(1).a1Source = -0.109000000;
CONDITIONS(1).a0Source =  0.545666700;
CONDITIONS(1).xmax = 60;
CONDITIONS(1).b2 = -0.0012595889;
CONDITIONS(1).b1 =  0.0733333;
CONDITIONS(1).b0 =  0.7057778;

% Dark empirical condition
CONDITIONS(2).name = "Dark";
CONDITIONS(2).a3Source = -2.000000e-5;
CONDITIONS(2).a2Source =  0.002000000;
CONDITIONS(2).a1Source = -0.062000000;
CONDITIONS(2).a0Source =  0.561000000;
CONDITIONS(2).xmax = 60;
CONDITIONS(2).b2 = -0.0009;
CONDITIONS(2).b1 =  0.038;
CONDITIONS(2).b0 =  1.287;

% Convert the supplied cubic drift coefficients
% [a3,a2,a1,a0] into [alpha,beta,xl,xh].
for conditionIndex = 1:numel(CONDITIONS)
    derived = derivePotentialParameters(CONDITIONS(conditionIndex));

    CONDITIONS(conditionIndex).alpha = derived.alpha;
    CONDITIONS(conditionIndex).beta  = derived.beta;
    CONDITIONS(conditionIndex).xl    = derived.xl;
    CONDITIONS(conditionIndex).xh    = derived.xh;
end

% Unscaled SDE: dX = f(X)dt + g(X)dW
ANALYSIS.driftScale = 1;
ANALYSIS.diffusionScale = 1;

ANALYSIS.parameterNames = ["alpha","beta","b0","b1","b2"];
ANALYSIS.parameterLabels = { ...
    '$\alpha$', '$\beta$', '$b_0$', '$b_1$', '$b_2$'};

% Relative finite-difference step for nominal local sensitivity
ANALYSIS.relativeSensitivityStep = 0.010;

% Figure style
STYLE.figureSizeIn = [7.0 5.25];
STYLE.tickFontSize = 28;
STYLE.labelFontSize = 32;
STYLE.annotationFontSize = 28;
STYLE.axisLineWidth = 2.6;

% Plot colors
COLORS.bright = [228 137 87]/255;
COLORS.dark   = [66 54 153]/255;

% Run normalized local sensitivity analysis
for conditionIndex = 1:numel(CONDITIONS)
    condition = CONDITIONS(conditionIndex);

    shouldRun = ...
        (condition.name == "Bright" && RUN.bright) || ...
        (condition.name == "Dark"   && RUN.dark);

    if shouldRun
        runNormalizedSensitivityAnalysis( ...
            condition,ANALYSIS,STYLE,COLORS);
    end
end





function runNormalizedSensitivityAnalysis( ...
        condition,analysis,style,colors)

    empiricalCondition = condition.name;
    x_l = condition.xl;
    x_h = condition.xh;
    x_max = condition.xmax;

    % theta = [alpha, beta, b0, b1, b2]
    theta0 = [ ...
        condition.alpha, ...
        condition.beta, ...
        condition.b0, ...
        condition.b1, ...
        condition.b2];

    driftScale = analysis.driftScale;
    diffusionScale = analysis.diffusionScale;
    parameterNames = analysis.parameterNames;
    parameterLabels = analysis.parameterLabels;
    relativeStep = analysis.relativeSensitivityStep;

    if driftScale <= 0
        error('analysis.driftScale must be strictly positive.');
    end

    if diffusionScale <= 0
        error('analysis.diffusionScale must be strictly positive.');
    end

    if relativeStep <= 0 || relativeStep >= 1
        error('analysis.relativeSensitivityStep must lie in (0,1).');
    end

    if ~isDiffusionPositive(theta0,x_max)
        error('%s nominal diffusion is not positive on [0,xmax].', ...
            empiricalCondition);
    end

    fprintf('\n%s empirical condition:\n',empiricalCondition);
    fprintf('SDE: dX = %.1f f(X)dt + %.1f g(X)dW\n', ...
        driftScale,diffusionScale);
    fprintf('alpha = %.15g\n',condition.alpha);
    fprintf('beta  = %.15g\n',condition.beta);
    fprintf('xl    = %.15g\n',condition.xl);
    fprintf('xh    = %.15g\n',condition.xh);

    nominalDiscriminant = peakDiscriminant( ...
        theta0,x_l,x_h,driftScale,diffusionScale);

    numberParameters = numel(theta0);
    scaledDiscriminantSensitivity = nan(1,numberParameters);

    for parameterIndex = 1:numberParameters
        thetaPlus = theta0;
        thetaMinus = theta0;

        thetaPlus(parameterIndex) = ...
            theta0(parameterIndex)*(1+relativeStep);
        thetaMinus(parameterIndex) = ...
            theta0(parameterIndex)*(1-relativeStep);

        % The perturbed diffusion must remain positive on the full state
        % interval when b0, b1, or b2 is varied.
        if parameterIndex >= 3
            if ~isDiffusionPositive(thetaPlus,x_max) || ...
               ~isDiffusionPositive(thetaMinus,x_max)
                error([ ...
                    '%s: the local perturbation of %s violates ', ...
                    'diffusion positivity.'], ...
                    empiricalCondition,parameterNames(parameterIndex));
            end
        end

        deltaPlus = peakDiscriminant( ...
            thetaPlus,x_l,x_h,driftScale,diffusionScale);
        deltaMinus = peakDiscriminant( ...
            thetaMinus,x_l,x_h,driftScale,diffusionScale);

        % Central relative finite difference:
        % |theta_i * d(Delta_pk)/d(theta_i)|
        scaledDiscriminantSensitivity(parameterIndex) = ...
            abs(deltaPlus-deltaMinus)/(2*relativeStep);
    end

    maximumSensitivity = max( ...
        scaledDiscriminantSensitivity,[],'omitnan');

    if isempty(maximumSensitivity) || ...
            ~isfinite(maximumSensitivity) || ...
            maximumSensitivity <= 0
        error('%s normalized sensitivity could not be computed.', ...
            empiricalCondition);
    end

    normalizedDiscriminantSensitivity = ...
        scaledDiscriminantSensitivity/maximumSensitivity;

    resultsTable = table( ...
        parameterNames', ...
        theta0', ...
        scaledDiscriminantSensitivity', ...
        normalizedDiscriminantSensitivity', ...
        'VariableNames',{ ...
            'Parameter', ...
            'NominalValue', ...
            'ScaledSensitivity', ...
            'NormalizedSensitivity'});

    fprintf('Nominal Delta_pk = %.12g\n',nominalDiscriminant);
    disp(resultsTable);

    if empiricalCondition == "Bright"
        conditionColor = colors.bright;
    else
        conditionColor = colors.dark;
    end

    fig = publicationFigure( ...
        style.figureSizeIn, ...
        sprintf('%s normalized local discriminant sensitivity', ...
        empiricalCondition));

    ax = axes('Parent',fig);

    bar(ax,normalizedDiscriminantSensitivity, ...
        'FaceColor',conditionColor, ...
        'EdgeColor',conditionColor, ...
        'LineWidth',style.axisLineWidth);

    grid(ax,'on');
    box(ax,'on');

    set(ax, ...
        'XTick',1:numberParameters, ...
        'XTickLabel',parameterLabels, ...
        'TickLabelInterpreter','latex');

    xlim(ax,[0.5,numberParameters+0.5]);
    ylim(ax,[0,1.08]);

    ylabel(ax,'$\widetilde{S}^{\Delta}_i$', ...
        'Interpreter','latex');

    applyPublicationAxesStyle(ax,style);
end


function coefficients = peakPolynomialCoefficients( ...
        theta,x_l,x_h,driftScale,diffusionScale)

    alpha = theta(1);
    beta  = theta(2);
    b0    = theta(3);
    b1    = theta(4);
    b2    = theta(5);

    effectiveDiffusionScale = ...
        diffusionScale/sqrt(driftScale);

    B0 = effectiveDiffusionScale*b0;
    B1 = effectiveDiffusionScale*b1;
    B2 = effectiveDiffusionScale*b2;

    eta3 = -2*(2*beta+B2^2);
    eta2 = 6*beta*(x_l+x_h)-3*B1*B2;
    eta1 = -2*beta*( ...
        x_l^2+x_h^2+4*x_l*x_h) ...
        -(B1^2+2*B0*B2);
    eta0 = -alpha ...
        +2*beta*(x_l^2*x_h+x_h^2*x_l) ...
        -B0*B1;

    coefficients = [eta3,eta2,eta1,eta0];
end


function deltaPk = peakDiscriminant( ...
        theta,x_l,x_h,driftScale,diffusionScale)

    coefficients = peakPolynomialCoefficients( ...
        theta,x_l,x_h,driftScale,diffusionScale);

    eta3 = coefficients(1);
    eta2 = coefficients(2);
    eta1 = coefficients(3);
    eta0 = coefficients(4);

    coefficientScale = max([ ...
        1,abs(eta3),abs(eta2),abs(eta1),abs(eta0)]);

    if abs(eta3) <= 100*eps(coefficientScale)
        error('The peak polynomial is not numerically cubic.');
    end

    pPk = ...
        (3*eta3*eta1-eta2^2)/(3*eta3^2);

    qPk = ...
        (2*eta2^3 ...
        -9*eta3*eta2*eta1 ...
        +27*eta3^2*eta0)/(27*eta3^3);

    % Normalized discriminant of the depressed cubic.
    deltaPk = -(4*pPk^3+27*qPk^2);
end


function positive = isDiffusionPositive(theta,xMax)

    b0 = theta(3);
    b1 = theta(4);
    b2 = theta(5);

    if any(~isfinite([b0,b1,b2]))
        positive = false;
        return
    end

    % A quadratic reaches its minimum on a compact interval at an endpoint
    % or, when applicable, at its interior vertex.
    candidateX = [0,xMax];

    coefficientScale = max([1,abs(b0),abs(b1),abs(b2)]);

    if abs(b2) > 100*eps(coefficientScale)
        vertexX = -b1/(2*b2);

        if vertexX > 0 && vertexX < xMax
            candidateX(end+1) = vertexX;
        end
    end

    gValues = b2*candidateX.^2+b1*candidateX+b0;
    positive = all(gValues > 0);
end


function fig = publicationFigure(figureSizeIn,figureName)

    fig = figure( ...
        'Color','w', ...
        'Units','inches', ...
        'Position',[1 1 figureSizeIn], ...
        'PaperPositionMode','auto', ...
        'Name',figureName, ...
        'NumberTitle','off');

    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 figureSizeIn];
    fig.PaperSize = figureSizeIn;
end


function applyPublicationAxesStyle(ax,style)

    ax.FontSize = style.tickFontSize;
    ax.FontWeight = 'bold';
    ax.LineWidth = style.axisLineWidth;
    ax.TickDir = 'out';

    ax.XLabel.FontSize = style.labelFontSize;
    ax.YLabel.FontSize = style.labelFontSize;
    ax.XLabel.FontWeight = 'bold';
    ax.YLabel.FontWeight = 'bold';
end


function condition = derivePotentialParameters(condition)
% Convert
%
% f(x) = a3*x^3 + a2*x^2 + a1*x + a0
%
% into
%
% V(x) = alpha*x + beta*(x-xl)^2*(x-xh)^2,
%
% with f(x) = -dV/dx.

    a3 = condition.a3Source;
    a2 = condition.a2Source;
    a1 = condition.a1Source;
    a0 = condition.a0Source;

    if abs(a3) < 1e-14
        error('%s: a3Source is too close to zero.',condition.name);
    end

    beta = -a3/4;

    % a2 = 6*beta*(xl+xh)
    rootSum = a2/(6*beta);

    % a1 = -2*beta*((xl+xh)^2 + 2*xl*xh)
    rootProduct = ( ...
        -a1/(2*beta)-rootSum^2)/2;

    rootDiscriminant = rootSum^2-4*rootProduct;

    tolerance = 1e-12*max([ ...
        1,abs(rootSum^2),abs(4*rootProduct)]);

    if rootDiscriminant < -tolerance
        error([ ...
            '%s: the cubic drift does not produce ', ...
            'two real xl and xh values.'],condition.name);
    end

    rootDiscriminant = max(rootDiscriminant,0);

    potentialRoots = sort([ ...
        (rootSum-sqrt(rootDiscriminant))/2, ...
        (rootSum+sqrt(rootDiscriminant))/2]);

    xl = potentialRoots(1);
    xh = potentialRoots(2);

    % a0 = -alpha + 2*beta*(xl+xh)*xl*xh
    alpha = -a0+2*beta*rootSum*rootProduct;

    condition.alpha = alpha;
    condition.beta = beta;
    condition.xl = xl;
    condition.xh = xh;

    recoveredCoefficients = [ ...
        -4*beta, ...
        6*beta*(xl+xh), ...
        -2*beta*(xl^2+xh^2+4*xl*xh), ...
        -alpha+2*beta*(xl^2*xh+xh^2*xl)];

    originalCoefficients = [a3,a2,a1,a0];
    maximumReconstructionError = max(abs( ...
        originalCoefficients-recoveredCoefficients));

    fprintf('\n%s derived potential parameters:\n',condition.name);
    fprintf('  alpha = %.15g\n',condition.alpha);
    fprintf('  beta  = %.15g\n',condition.beta);
    fprintf('  xl    = %.15g\n',condition.xl);
    fprintf('  xh    = %.15g\n',condition.xh);
    fprintf('  maximum coefficient error = %.3e\n', ...
        maximumReconstructionError);
end
