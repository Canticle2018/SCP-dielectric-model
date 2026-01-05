function fit_cc_quadratic_c()
% Single-relaxation Cole–Cole model
% DE(c), tau(c), and beta(c) are parameterized as quadratic functions of salinity c (mol/L).
% eps_inf is treated as a global parameter and optimized jointly
% (einf_init is used only for initialization and is not fixed).
% Frequency band is defined by FMIN_GHZ / FMAX_GHZ.
% Objective: mean squared relative error + mild endpoint/interval-boundary
% penalties + a very small monotonic trend regularization.

    clc; clear;

    % ====== Frequency band settings (modify here globally) ======
    FMIN_GHZ = 1;      % Lower frequency bound (GHz)
    FMAX_GHZ = 6;      % Upper frequency bound (GHz)
    band.min_hz = FMIN_GHZ * 1e9;
    band.max_hz = FMAX_GHZ * 1e9;

    % ====== Concentration settings (one-to-one with data files) ======
    % c is used for fitting and modeling (mol/L) and serves as the
    % independent variable in the quadratic parameterization.
    c = [0.11086, 0.19193, 0.40412, 0.59740, ...
         0.79496, 0.97854, 1.4581, 1.9381, ...
         2.3937, 2.8316, 3.7118, 4.5212]';  % mol/L

    % concs is used only for parsing file names (g/L).
    concs = [6.52, 11.34, 24.18, 36.15, 48.65, 60.11, 92.70, 126.63, ...
             160.49, 194.59, 268.01, 341.59]; % g/L
    concs = sanitize_numeric_column(concs);

    % Consistency check
    M = numel(concs);
    assert(M == numel(c), ...
        'The lengths of c and concs must match and correspond one-to-one.');


    cmin = min(c);  cmax = max(c);   % --> Use the range of c

    % ====== Data loading ======
    ds.freq = cell(M,1); ds.epsr = cell(M,1); ds.epsi = cell(M,1);
    for i = 1:M
        CgL = concs(i);  % g/L values encoded in the file names
        fname = fullfile('D:', 'Dielectric_Permittivity', 'solution', 'no_sigma', ...
                         sprintf('%.2f_LC_25℃_avg.csv', CgL));
        data = csvread(fname, 12, 0);       
        ds.freq{i} = data(:,1) * 1e9;       
        ds.epsr{i} = data(:,2);
        ds.epsi{i} = data(:,3);
    end

    % ====== range ======
    rg.DE   = [5, 90];             % Δε
    rg.tau  = [6.0e-12, 9.0e-12];  % tau range
    rg.beta = [0.00, 0.35];        % beta
    rg.einf = [3.0, 8.0];

    % ====== Tabulated priors for initialization: quadratic polyfit over c (mol/L) ======
    eps_s = [ ...
        (76.48+76.79)/2, (75.45+75.71)/2, (72.59+73.22)/2, (70.29+70.80)/2, ...
        (68.22+68.75)/2, (66.46+66.82)/2, (62.18+62.50)/2, (58.21+58.39)/2, ...
        (55.24+55.08)/2, (52.43+52.75)/2, (47.37+47.99)/2, (43.48+44.63)/2 ]';
    tau_ps = [ ...
        (8.02+8.07)/2, (8.01+7.99)/2, (7.81+7.92)/2, (7.69+7.72)/2, ...
        (7.57+7.69)/2, (7.41+7.52)/2, (7.35+7.32)/2, (7.08+7.04)/2, ...
        (6.86+6.96)/2, (6.74+7.04)/2, (6.63+7.14)/2, (6.29+6.70)/2 ]';
    beta0 = [ ...
        (0.000+0.016)/2, (0.024+0.010)/2, (0.009+0.029)/2, (0.021+0.030)/2, ...
        (0.038+0.025)/2, (0.045+0.036)/2, (0.032+0.044)/2, (0.055+0.059)/2, ...
        (0.088+0.080)/2, (0.112+0.085)/2, (0.150+0.114)/2, (0.199+0.196)/2 ]';

    % eps_inf used only as an initial guess (will be optimized, not fixed)

    einf_init = 5.65;
    DE_tab    = max(eps_s - einf_init, rg.DE(1));
    tau_tab   = tau_ps * 1e-12;
    beta_tab  = beta0;

    % ——> Use c as the independent variable
    pDE   = polyfit(c, DE_tab,  2);  % return [a2 a1 a0]
    pTAU  = polyfit(c, tau_tab, 2);
    pBE   = polyfit(c, beta_tab,2);

    % ====== Parameter vector: theta = [a0 a1 a2 b0 b1 b2 g0 g1 g2 z_einf] ======
    invsig = @(y,lo,hi) log((y-lo)./(hi-y));
    theta0 = [ pDE(3) pDE(2) pDE(1), ...
               pTAU(3) pTAU(2) pTAU(1), ...
               pBE(3)  pBE(2)  pBE(1), ...
               invsig(einf_init, rg.einf(1), rg.einf(2)) ]';

    % ====== Optimization ======
    opts = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',5e4, ...
                        'MaxIterations', 1000,'SpecifyObjectiveGradient',false);
    problem.objective = @(th) obj_quad_c(th, ds, c, [cmin cmax], rg, band); % ——> 传入 c
    problem.x0        = theta0;
    problem.solver    = 'fmincon';
    problem.options   = opts;
    L = [-Inf*ones(9,1); -20];
    U = [ Inf*ones(9,1);  20];
    problem.lb = L; problem.ub = U;

    [theta_opt, ~] = fmincon(problem);

    % ====== Unpack, report, and save results (including eps_inf saved separately) ======
    params = unpack_quad(theta_opt, rg);

    fprintf('\n== frequency range：%.2g–%.2g GHz ==\n', FMIN_GHZ, FMAX_GHZ);
    fprintf('== global eps_inf: %.6f ==\n', params.einf);
    fprintf('DE(c)   = %.6f  + (%.6f)*c  + (%.6e)*c^2\n', params.a0, params.a1, params.a2);
    fprintf('tau(c)  = %.6e + (%.6e)*c + (%.6e)*c^2   [unit: s]\n', params.b0, params.b1, params.b2);
    fprintf('beta(c) = %.6f  + (%.6f)*c  + (%.6e)*c^2\n', params.g0, params.g1, params.g2);

    % ——> Compute parameters for each concentration using c (tau reported in ps)
    DE_seq  = params.a0 + params.a1*c + params.a2*(c.^2);
    tau_seq = params.b0 + params.b1*c + params.b2*(c.^2);
    be_seq  = params.g0 + params.g1*c + params.g2*(c.^2);

    % Save the main table: first column DE; second column Tau (ps); third column Beta
    T = [DE_seq, tau_seq*1e12, be_seq];
    dlmwrite('best_quad_cc_seq_c.txt', T, 'delimiter','\t','precision', 8);

    % Save eps_inf to a separate file for convenient access in plotting scripts
    dlmwrite('best_quad_epsinf.txt', params.einf, 'delimiter','\t','precision', 12);

    % Save the fitted equations
    fid = fopen('best_quad_cc_equations.txt','w');
    fprintf(fid, 'band = %.2g–%.2g GHz\n', FMIN_GHZ, FMAX_GHZ);
    fprintf(fid, 'eps_inf = %.12f\n', params.einf);
    fprintf(fid, 'DE(c)   = %.12f + (%.12f)*c + (%.12e)*c^2\n', params.a0, params.a1, params.a2);
    fprintf(fid, 'tau(c)  = %.12e + (%.12e)*c + (%.12e)*c^2   (unit: seconds; table saved in ps)\n', ...
                 params.b0, params.b1, params.b2);
    fprintf(fid, 'beta(c) = %.12f + (%.12f)*c + (%.12e)*c^2\n', params.g0, params.g1, params.g2);
    fclose(fid);

    % Optionally print a copy to the console for easy verification
    disp('output [DE, tau(ps), beta]：');
    disp(T);
end

% ================= Utility / objective functions =================
function v = sanitize_numeric_column(x)
    if iscell(x),      v = str2double(string(x));
    elseif isstring(x),v = str2double(x);
    elseif iscategorical(x), v = str2double(string(x));
    else,              v = x;
    end
    v = double(v(:));
    assert(all(isfinite(v)) && ~isempty(v), 'concs 需要是有限的数值向量');
end

function y = sigmoid(x)
    y = 1./(1+exp(-x));
end

function h = hinge(x)
    h = max(0,x).^2;
end

function params = unpack_quad(theta, rg)
% theta = [a0 a1 a2 b0 b1 b2 g0 g1 g2 z_einf]
    params.a0 = theta(1); params.a1 = theta(2); params.a2 = theta(3);
    params.b0 = theta(4); params.b1 = theta(5); params.b2 = theta(6);
    params.g0 = theta(7); params.g1 = theta(8); params.g2 = theta(9);
    z_einf    = theta(10);
    params.einf = rg.einf(1) + (rg.einf(2)-rg.einf(1)) * sigmoid(z_einf);
end

function [J, grad] = obj_quad_c(theta, ds, cvec, crange, rg, band)
% Objective: mean squared relative error + mild interval-boundary penalties
% + a very small monotonic trend regularization + a small high-frequency anchor.
% Note: the independent variable here is c (mol/L).
    params = unpack_quad(theta, rg);
    cmin = crange(1); cmax = crange(2);

    % Primary error term (restricted to the specified frequency band)
    M = numel(ds.freq);
    sse = zeros(M,1);

    for k = 1:M
        f_all = ds.freq{k};
        mask = (f_all >= band.min_hz) & (f_all <= band.max_hz);
        f = f_all(mask); w = 2*pi*f;

        eps_meas = ds.epsr{k}(mask) - 1i*ds.epsi{k}(mask);

        ck  = cvec(k);                       % ——> 用 c(k)
        DE = params.a0 + params.a1*ck + params.a2*ck^2;
        t  = params.b0 + params.b1*ck + params.b2*ck^2;
        be = params.g0 + params.g1*ck + params.g2*ck^2;

        eps_model = params.einf + DE ./ (1 + (1i*w*t).^(1 - be));

        res_re = real(eps_model) - real(eps_meas);
        res_im = imag(eps_model) - imag(eps_meas);

        denom_re = max(abs(real(eps_meas)), 1e-6);
        denom_im = max(abs(imag(eps_meas)), 1e-6);

        rel2 = sum( (abs(res_re)./denom_re).^2 + (abs(res_im)./denom_im).^2 );
        dof  = max(2*numel(f) - 4, 1);
        sse(k) = rel2 / dof;

        % Small high-frequency anchor (last 10% of in-band frequencies) to suppress
        % compensatory drift in eps_inf
        Nhf = max(3, round(0.1*numel(f)));
        idx = (numel(f)-Nhf+1):numel(f);
        hf_bias = mean( (real(eps_model(idx)) - real(eps_meas(idx))) ./ ...
                        max(abs(real(eps_meas(idx))),1e-6) );
        sse(k) = sse(k) + 0.01*hf_bias^2;
    end

    J = mean(sse);

    % Interval boundary/range penalty: sample on a grid over c ∈ [cmin, cmax]
    % to prevent the quadratic terms from exceeding plausible bounds.
    cgrid = linspace(cmin, cmax, 21);
    DEg = params.a0 + params.a1*cgrid + params.a2*(cgrid.^2);
    tg  = params.b0 + params.b1*cgrid + params.b2*(cgrid.^2);
    beg = params.g0 + params.g1*cgrid + params.g2*(cgrid.^2);

    w_range = 1e-3;
    J = J + w_range * ( sum(hinge(DEg - rg.DE(2)) + hinge(rg.DE(1) - DEg)) ...
                      + sum(hinge(tg  - rg.tau(2)) + hinge(rg.tau(1) - tg )) ...
                      + sum(hinge(beg - rg.beta(2))+ hinge(rg.beta(1) - beg)) );

    % % Very small monotonic trend regularization (encourages DE↓, tau↓, beta↑);
    % set w_mono = 0 to disable if undesired.
    dDE = params.a1 + 2*params.a2*cgrid;   % d(DE)/dc
    dt  = params.b1 + 2*params.b2*cgrid;
    dbe = params.g1 + 2*params.g2*cgrid;

    w_mono = 1e-3;
    J = J + w_mono * ( sum(hinge(dDE)) + sum(hinge(dt)) + sum(hinge(-dbe)) );

    grad = [];  % Use numerical gradients
end
