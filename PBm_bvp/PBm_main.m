%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Batch EDL solver (objective: err_record = normalized Δσ only)
%
% Workflow:
%  1) First pass: use parfor to run all samples once with the initial x_max;
%  2) Refinement: select samples with err >= threshold, and reprocess them
%     using a standard for-loop by scanning
%           x_max ∈ [4.0, 1.5] nm with a step of 0.1 nm,
%     retaining the x_max that yields the minimum err.
%
% Outputs:
%  - result/ss_sm.txt:
%       out_data for each sample
%  - all_eps_r0.txt:
%       [ss, sm, eps_r0]
%  - all_dthick_nm.txt:
%       [ss, sm, d_thick_nm]
%  - all_objs.txt:
%       each row [ss, sm, err_record]
%  - all_delta_sigma_xmax.txt:
%       [ss, sm, err_record, x_max_nm]
%  - Command window:
%       prints err for all samples and reports the minimum and maximum values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;

% ======================= Basic parameters (DY sample as an example) =======================
T    = 298.15;      % Temperature (K)
Sand = 8.8;         % Sand content (%)
Clay = 8.3;         % Clay content (%)

ERR_THRESH       = 0.005;    % Acceptance threshold: err_record < 1%
XMAX_SCAN_START  = 4.0e-9;  % Scan start point (m)
XMAX_SCAN_STOP   = 1.5e-9;  % Scan end point (m)
XMAX_SCAN_STEP   = -0.1e-9; % Scan step size (m)

% ======================= Parallel configuration (used for the first pass only) =======================
ENABLE_PARALLEL = true;
PAR_WORKERS     = 8;

useParallel = false;
pool = gcp('nocreate');
if ENABLE_PARALLEL
    try
        if isempty(pool) || pool.NumWorkers ~= PAR_WORKERS
            if ~isempty(pool), delete(pool); end
            parpool('local', PAR_WORKERS);
        end
        useParallel = true;
        fprintf('[INFO] Parallel pool is active with %d workers.\n', gcp('nocreate').NumWorkers);

    catch ME
        warning('[WARN] Parallel pool initialization failed: %s\n[WARN] Reverting to serial execution.', ME.message);
        if ~isempty(gcp('nocreate')), delete(gcp('nocreate')); end
        useParallel = false;
    end
else
    if ~isempty(pool)
        delete(pool);
        fprintf('[INFO] Parallel pool shut down successfully.\n');
    end
    useParallel = false;
    fprintf('[INFO] Parallel pool is closed; falling back to serial for-loop execution.\n');
end

% ======================= 数据读取 =======================
% all_mol.txt: [ss(‰), sm(%), c0(mol/L)]
data   = load('all_mol.txt');
n      = size(data, 1);
c0_vec = data(:,3);

% ======================= Output containers =======================
output_data_all = cell(n, 1);
err_vec         = nan(n, 1);   % = err_record（scalar; normalized Δσ）
eps_r0_vec      = nan(n, 1);
dthick_nm_vec   = nan(n, 1);
meta_eps        = nan(n, 3);   % [ss sm eps_r0]
meta_d          = nan(n, 3);   % [ss sm d_thick_nm]
xmax_used       = zeros(n, 1); % final x_max（m）

% ======================= output =======================
current_folder = pwd;
output_folder  = fullfile(current_folder, 'result/');
if ~exist(output_folder, 'dir'); mkdir(output_folder); end

% ======================= Initialize x_max (segmented by c0; used for the initial run only) =======================
xmax_init = zeros(n,1);
for i = 1:n
    c0 = c0_vec(i);
    if c0 < 1
        xmax_init(i) = 8e-9;
    elseif c0 < 3
        xmax_init(i) = 3e-9;
    else
        xmax_init(i) = 2.0e-9;
    end
end

% ======================= Pass-0: Initial pass (parfor) =======================
fprintf('==== Pass-0: Full initial computation (n=%d) ====\n', n);
t0 = tic;
[output_data_all, eps_r0_vec, dthick_nm_vec, err_vec, ...
 meta_eps, meta_d, xmax_used] = ...
    run_batch_once(1:n, useParallel, data, Sand, Clay, T, xmax_init, ...
                   output_data_all, eps_r0_vec, dthick_nm_vec, err_vec, ...
                   meta_eps, meta_d, xmax_used);
fprintf('Pass-0 completed in %.2f s.\n', toc(t0));

% ======================= 修正阶段：串行扫描（for 循环逐个修正） =======================
bad_idx = find(err_vec >= ERR_THRESH | ~isfinite(err_vec));
fprintf('Number of samples requiring scan-based refinement: %d / %d (threshold %.2f%%)\n', ...
        numel(bad_idx), n, 100*ERR_THRESH);
fprintf('  Scan range: %.1f nm to %.1f nm with step %.1f nm\n', ...
        XMAX_SCAN_START*1e9, XMAX_SCAN_STOP*1e9, XMAX_SCAN_STEP*1e9);


for jj = 1:numel(bad_idx)
    m = bad_idx(jj);
    [out_m, eps0_m, dnm_m, err_m, meta_eps_m, meta_d_m, x_used_m] = ...
        scan_best_by_err(m, data, Sand, Clay, T, ...
                         XMAX_SCAN_START, XMAX_SCAN_STOP, XMAX_SCAN_STEP);

    % If the scan fails (err_m is non-finite), retain the initial-pass result and issue a warning
    if ~isfinite(err_m)
        warning('Sample m=%d: no feasible solution found during scan; retaining initial-pass result.', m);
        continue;
    end

    % Write back the optimal result
    output_data_all{m} = out_m;
    eps_r0_vec(m)      = eps0_m;
    dthick_nm_vec(m)   = dnm_m;
    err_vec(m)         = err_m;
    meta_eps(m,:)      = meta_eps_m;
    meta_d(m,:)        = meta_d_m;
    xmax_used(m)       = x_used_m;
end

% ======================= Summary statistics and extrema =======================
fprintf('\n===== Summary of delta_sigma (err_record) for all samples =====\n');
for m = 1:n
    ss = data(m,1); sm = data(m,2);
    fprintf('m=%d (ss=%g, sm=%g) -> err=%.9g\n', m, ss, sm, err_vec(m));
end
mask = isfinite(err_vec);
if any(mask)
    idx_valid = find(mask);
    [vmin, kmin] = min(err_vec(mask)); imin = idx_valid(kmin);
    [vmax, kmax] = max(err_vec(mask)); imax = idx_valid(kmax);
    fprintf('---- min err ----\n');
    fprintf('m=%d (ss=%g, sm=%g) -> err=%.9g, x_max=%.6f nm\n', ...
        imin, data(imin,1), data(imin,2), vmin, xmax_used(imin)*1e9);
    fprintf('---- max err ----\n');
    fprintf('m=%d (ss=%g, sm=%g) -> err=%.9g, x_max=%.6f nm\n', ...
        imax, data(imax,1), data(imax,2), vmax, xmax_used(imax)*1e9);
end
fprintf('===== all err values =====\n\n');

% ======================= Unified saving =======================
% 1) Spatial profiles (out_data) for each sample
for m = 1:n
    ss = data(m,1); sm = data(m,2);
    out_data = output_data_all{m};
    if isnumeric(out_data) && ~isempty(out_data)
        filename = sprintf('%s%d_%d.txt', output_folder, ss, sm);
        save(filename, 'out_data', '-ascii');
    else
        warning('Skipping save for m=%d: out_data null。', m);
    end
end

% 2) Export tabulated values of eps_r0 and d_thick
output_eps_r0    = meta_eps;   % [ss  sm  eps_r0]
output_dthick_nm = meta_d;     % [ss  sm  d_thick_nm]
save('all_eps_r0.txt',    'output_eps_r0',    '-ascii');
save('all_dthick_nm.txt', 'output_dthick_nm', '-ascii');

% 3) every row：ss sm err_record
fid = fopen('all_objs.txt','w');
for m = 1:n
    ss = data(m,1); sm = data(m,2);
    fprintf(fid, '%g %g %.9g\n', ss, sm, err_vec(m));
end
fclose(fid);

% 4)  err_record 与 x_max（nm）corresponding to (ss, sm) 
fid2 = fopen('all_delta_sigma_xmax.txt','w');
for m = 1:n
    ss = data(m,1); sm = data(m,2);
    fprintf(fid2, '%.6g %.6g %.9g %.6f\n', ss, sm, err_vec(m), xmax_used(m)*1e9);
end
fclose(fid2);

fprintf('all saved.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subfunctions below %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Single call: given x_max, return out_data / epsr0 / d / err, etc.
function [out_data, eps_r0_final, dthick_nm, err_scalar, ...
          meta_row_eps, meta_row_d, x_used] = ...
    compute_once_with_xmax(m, data, Sand, Clay, T, x_in)

    ss_mg = data(m,1); sm_pct = data(m,2); c0 = data(m,3);
    ss = ss_mg / 1000; sm = sm_pct / 100;
    [sigma0, eps_r] = get_para(Sand, Clay, c0);

    % Call the solver function (err_record = normalized Δσ; may be a scalar or a vector)
    [out_data, eps_r0_final, err_record, d_thick_m] = ...
        new_solveEDL_sigma_iter_iterative( ...
            x_in, ss*1000, sm*100, sigma0, c0, eps_r, T);

    % Use only the last err value
    if numel(err_record) >= 1
        err_scalar = err_record(end);
    else
        err_scalar = NaN;
    end

    dthick_nm    = d_thick_m * 1e9;
    meta_row_eps = [ss_mg, sm_pct, eps_r0_final];
    meta_row_d   = [ss_mg, sm_pct, dthick_nm];
    x_used       = x_in;
end

% Initial batch pass: first evaluation for idx_list using their respective x_max values (parfor inside)
function [output_data_all, eps_r0_vec, dthick_nm_vec, err_vec, ...
          meta_eps, meta_d, xmax_used] = ...
    run_batch_once(idx_list, useParallel, data, Sand, Clay, T, xmax_init, ...
                   output_data_all, eps_r0_vec, dthick_nm_vec, err_vec, ...
                   meta_eps, meta_d, xmax_used)

    if isempty(idx_list), return; end
    K = numel(idx_list);

    tmp_out   = cell(K,1);
    tmp_eps0  = nan(K,1);
    tmp_dnm   = nan(K,1);
    tmp_err   = nan(K,1);
    tmp_metae = nan(K,3);
    tmp_metad = nan(K,3);
    tmp_x     = nan(K,1);

    if useParallel
        parfor kk = 1:K
            m = idx_list(kk);
            [tmp_out{kk}, tmp_eps0(kk), tmp_dnm(kk), tmp_err(kk), ...
             tmp_metae(kk,:), tmp_metad(kk,:), tmp_x(kk)] = ...
                 compute_once_with_xmax(m, data, Sand, Clay, T, xmax_init(m));
        end
    else
        for kk = 1:K
            m = idx_list(kk);
            [tmp_out{kk}, tmp_eps0(kk), tmp_dnm(kk), tmp_err(kk), ...
             tmp_metae(kk,:), tmp_metad(kk,:), tmp_x(kk)] = ...
                 compute_once_with_xmax(m, data, Sand, Clay, T, xmax_init(m));
        end
    end

    for kk = 1:K
        m = idx_list(kk);
        output_data_all{m} = tmp_out{kk};
        eps_r0_vec(m)      = tmp_eps0(kk);
        dthick_nm_vec(m)   = tmp_dnm(kk);
        err_vec(m)         = tmp_err(kk);
        meta_eps(m,:)      = tmp_metae(kk,:);
        meta_d(m,:)        = tmp_metad(kk,:);
        xmax_used(m)       = tmp_x(kk);
    end
end

% Serial scan: x_max ∈ [x_start, x_stop] (step x_step < 0), select the minimum err
function [out_best, eps0_best, dnm_best, err_best, meta_eps_best, meta_d_best, x_best] = ...
    scan_best_by_err(m, data, Sand, Clay, T, x_start, x_stop, x_step)

    xs = x_start:x_step:x_stop;   
    if xs(end) ~= x_stop
        xs = [xs, x_stop];        % Ensure inclusion of the endpoint
    end

    out_best      = [];
    eps0_best     = NaN;
    dnm_best      = NaN;
    err_best      = Inf;
    meta_eps_best = [data(m,1), data(m,2), NaN];
    meta_d_best   = [data(m,1), data(m,2), NaN];
    x_best        = NaN;

    for j = 1:numel(xs)
        x_try = xs(j);
        try
            [out_j, eps0_j, dnm_j, err_j, meta_eps_j, meta_d_j, ~] = ...
                compute_once_with_xmax(m, data, Sand, Clay, T, x_try);
        catch
            err_j = Inf;  
        end

        if isfinite(err_j) && (err_j < err_best)
            err_best      = err_j;
            out_best      = out_j;
            eps0_best     = eps0_j;
            dnm_best      = dnm_j;
            meta_eps_best = meta_eps_j;
            meta_d_best   = meta_d_j;
            x_best        = x_try;
        end
    end

    % If all attempts fail, return NaN/empty; the caller will retain the initial-pass result and issue a warning
    if ~isfinite(err_best)
        out_best = [];
        eps0_best = NaN; dnm_best = NaN; x_best = NaN;
        meta_eps_best(3) = NaN; meta_d_best(3) = NaN;
    end
end

% Empirical relations for obtaining sigma0 and bulk-phase eps_r (for use by the new solver)
function [sigma0, eps_r] = get_para(Sand, Clay, c0)
    % SSA
    SSA = 0.0360 * Clay^2.08 + 49.73;   % m^2/g
    As_m2kg = SSA * 1000;               % m^2/kg
    % CEC
    CEC_cmolkg = 0.0067 * Clay^2.00 + 12.25; % cmol/kg
    CEC_molkg  = CEC_cmolkg * 1e-2;          % mol/kg
    % sigma0 (C/m^2)
    F = 96485;
    sigma0 = CEC_molkg * F / As_m2kg;
    % Gavish εr
    eps_b0 = 78.54; eps_ms = 30.08; alpha  = 11.5;
    v = 3 * alpha * c0 / (eps_b0 - eps_ms);
    L = coth(v) - 1 / v;
    eps_r = eps_b0 - (eps_b0 - eps_ms) * L;
end
