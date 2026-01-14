clc; clear;

%% =========================================================
%  SCP LUT generator (L/C band)
%  Output: LUT/DY.txt, C03.txt, ...
% =========================================================

% ---------- Frequency bands ----------
Bands = [1.4, 5.4];   % GHz (L / C)

% ---------- Fixed model parameters ----------
x_param = struct( ...
    'alpha',   0.586, ...
    'a_rho',   0.00759, ...
    'a_sigma', 2.253, ...
    'p_sigma', 4, ...
    'k',       0.5, ...
    'd',       0 );

disp('Using fixed SCP parameters:');
disp(x_param);

% ---------- Samples ----------
sampleNames = {'DY','C03','C10','C24','C26','C28'};
samples = build_samples(sampleNames);

% ---------- LUT output directory ----------
lutDir = fullfile(pwd,'LUT');
if ~exist(lutDir,'dir')
    mkdir(lutDir);
end

%% =========================================================
%  Main loop: sample → (SS,SM,c0) → eps*
%% =========================================================
for si = 1:numel(samples)
    S = samples{si};

    % Load (SS, SM, c0)
    DATA = load_sample_mol(S.name);   % [SS, SM, c0]
    nRow = size(DATA,1);

    % Output matrix
    % [SS SM ReL ImL ReC ImC]
    outMat = zeros(nRow, 6);

    for i = 1:nRow
        ss = DATA(i,1);
        sm = DATA(i,2);
        c0 = DATA(i,3);

        eps_L = eps_calcu_sample(S, ss, sm, c0, Bands(1), x_param);
        eps_C = eps_calcu_sample(S, ss, sm, c0, Bands(2), x_param);

        outMat(i,:) = [
            ss, sm, ...
            real(eps_L), abs(imag(eps_L)), ...
            real(eps_C), abs(imag(eps_C))
        ];
    end

    % ---------- Write file ----------
    outFile = fullfile(lutDir, [S.name,'.txt']);
    fid = fopen(outFile,'w');

    % Header (FIRST LINE)
    fprintf(fid, ...
        'SS\tSM\teps_real_L\teps_imag_L\teps_real_C\teps_imag_C\n');

    % Data
    % for k = 1:nRow
    %     fprintf(fid, ...
    %         '%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n', outMat(k,:));
    % end

    for k = 1:nRow
        fprintf(fid, ...
            '%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
            round(outMat(k,1)), round(outMat(k,2)), outMat(k,3:6));
    end


    fclose(fid);
    fprintf('✔ LUT generated: %s (%d rows)\n', outFile, nRow);
end

disp('All LUT files generated successfully.');

%% =========================================================
%  ======== Functions (UNMODIFIED CORE) ========
%% =========================================================

function eps_soil = eps_calcu_sample(S, ss, sm, c0, Band, x)
    Sand = S.Sand; Clay = S.Clay;
    rho_b_eff = S.rho_b - x.a_rho * S.TOC;

    [~,~,SSA] = get_para(Sand,Clay,c0);

    fname = sprintf('%g_%g.txt', ss, sm);
    filePath = fullfile(pwd,S.name,'test',fname);
    A = load(filePath);

    xgrid   = A(:,1);
    c_plus  = A(:,3);
    c_minus = A(:,4);
    E       = abs(A(:,5));
    Nor_E   = E / max(E);

    [r1_f,r2_f,s_f] = DCC(c0,Band);
    eps_f = 5.65 + r1_f + r2_f + s_f;

    idx2 = find(Nor_E < 0.005,1);
    if isempty(idx2), idx2 = numel(xgrid); end

    d_max = sm*0.01 / (SSA*rho_b_eff*1e6);
    [~,idx1] = min(abs(xgrid-d_max));
    stop_idx = max(2, min(idx1, idx2));

    eps_b = 0; d0 = 0; tiny = 1e-12;

    for m = 2:stop_idx
        c_eff = c_minus(m) + x.k*(c_plus(m)-c_minus(m)) + x.d;
        c_eff = max(c_eff, tiny);

        [r1,r2,s] = DCC(c_eff,Band);
        g_sigma = exp(-x.a_sigma * Nor_E(m).^x.p_sigma);
        eps_loc = 5.65 + r1 + r2 + g_sigma*s;

        dx = xgrid(m) - xgrid(m-1);
        eps_b = eps_b + eps_loc * dx;
        d0 = d0 + dx;
    end

    if d0 == 0
        eps_soil = eps_f;
        return;
    end

    eps_b = eps_b / d0;
    W_b = SSA*d0*rho_b_eff*1e6;
    W_f = max(0, sm/100 - W_b);

    eps_soil = mix_eps_alpha(eps_b, eps_f, W_b, W_f, rho_b_eff, x.alpha);
end

function [relax1, relax2, s] = DCC(c, Band)
    epsilon0 = 8.854187817e-12;
    omega = 2*pi*Band*1e9;

    a0=71.848030960126; a1=-11.710697339478; a2=9.705806851115e-01;
    b0=8.081607866047e-12; b1=-6.202183992817e-13; b2=6.549275688871e-14;
    g0=0.013517149127; g1=0.012928376385; g2=5.867094912158e-03;

    DE   = a0 + a1*c + a2*c.^2;
    tau  = b0 + b1*c + b2*c.^2;
    beta = g0 + g1*c + g2*c.^2;

    relax1 = DE ./ (1 + (1i*omega*tau).^(1 - beta));
    relax2 = 0;

    sigma = 9.65*c - 0.89*c.^1.94;
    s = sigma ./ (1i*omega*epsilon0);
end

function ep_soil = mix_eps_alpha(ep_b, ep_f, W_b, W_f, rho_b, alpha_mix)
    V_b = W_b;
    V_ss = rho_b / 2.65;
    V_f = min(W_f, 1 - V_ss - V_b);
    V_air = 1 - V_ss - V_f - V_b;

    ep_ss = 4.7 - 1i*0.2;
    ep_a  = 1;

    s_e = ep_a^alpha_mix*V_air + ep_ss^alpha_mix*V_ss + ...
          ep_b^alpha_mix*V_b  + ep_f^alpha_mix*V_f;

    ep_soil = s_e.^(1/alpha_mix);
end

function [sigma0, eps_r, SSA] = get_para(Sand, Clay, c0)
    SSA = 0.0360 * Clay^2.08 + 49.73;
    As_m2kg = SSA*1000;

    CEC = (0.0067 * Clay^2.00 + 12.25)*1e-2;
    F = 96485;
    sigma0 = CEC*F/As_m2kg;

    eps_b0 = 78.54; eps_ms = 30.08; alpha = 11.5;
    v = 3*alpha*c0/(eps_b0-eps_ms);
    L = coth(v) - 1/v;
    eps_r = eps_b0 - (eps_b0-eps_ms)*L;
end

function DATA = load_sample_mol(sampleName)
    persistent CACHE
    if isempty(CACHE)
        CACHE = containers.Map('KeyType','char','ValueType','any');
    end
    if isKey(CACHE, sampleName)
        DATA = CACHE(sampleName);
        return;
    end

    mol_path = fullfile(pwd, sampleName, 'all_mol.txt');
    if ~exist(mol_path,'file')
        error('Missing all_mol.txt for sample %s', sampleName);
    end
    DATA = load(mol_path);
    CACHE(sampleName) = DATA;
end

function samples = build_samples(names)
    db = struct( ...
        'DY',  struct('name','DY', 'Sand',8.8, 'Clay',8.3, 'rho_b',1.40,'TOC',0.05), ...
        'C03', struct('name','C03','Sand',67.9,'Clay',8.9, 'rho_b',1.47,'TOC',0.80), ...
        'C10', struct('name','C10','Sand',29.4,'Clay',24.5,'rho_b',1.59,'TOC',2.70), ...
        'C24', struct('name','C24','Sand',17.4,'Clay',19.2,'rho_b',1.57,'TOC',1.74), ...
        'C26', struct('name','C26','Sand',18.5,'Clay',15.9,'rho_b',1.46,'TOC',1.13), ...
        'C28', struct('name','C28','Sand',10.5,'Clay',31.3,'rho_b',1.39,'TOC',1.30) );

    samples = cell(numel(names),1);
    for i = 1:numel(names)
        samples{i} = db.(names{i});
    end
end
