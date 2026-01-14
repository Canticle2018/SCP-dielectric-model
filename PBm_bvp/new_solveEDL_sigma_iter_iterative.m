%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Core function for solving BVP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out_data, eps_r0_final, err_record, d_thick] = new_solveEDL_sigma_iter_iterative( ...
    xmax, ss, sm, sigma0, c0_molL, eps_r_unused, T)

    % ===== constant =====
    e    = 1.602e-19;          % C
    kB   = 1.381e-23;          % J/K
    NA   = 6.022e23;           % 1/mol
    eps0 = 8.854e-12;          % F/m

    sigma_s = -abs(sigma0);    % Free surface charge density with sign (C/m^2)

    % mol/L -> mol/m^3
    c0 = c0_molL * 1000;

    % Ionic volume (Bikerman/Fermi; sphere volume divided by the close-packing limit of 0.74)
    a_Na = 0.716e-9; v_Na = (pi/6) * a_Na^3 / 0.74; 
    a_Cl = 0.664e-9; v_Cl = (pi/6) * a_Cl^3 / 0.74;

    % ===== solve BVP =====
    xmesh   = linspace(0, xmax, 3000);
    initFun = @(x) [ -0.05*(1 - x/xmax),  sigma_s*(1 - x/xmax) ];
    solinit = bvpinit(xmesh, initFun);
    opts    = bvpset('RelTol',1e-5,'AbsTol',1e-7,'NMax',20000);
    sol     = bvp4c(@edlODE, @edlBC, solinit, opts);

    % ===== Post-processing =====
    x_fine = linspace(0, xmax, 3000);
    Y      = deval(sol, x_fine);
    phi    = Y(1, :);           % V
    D      = Y(2, :);           % C/m^2

    N = numel(x_fine);
    E        = zeros(1,N);      % V/m
    rho_full = zeros(1,N);      % C/m^3
    c_plusL  = zeros(1,N);      % mol/L
    c_minusL = zeros(1,N);      % mol/L

    for m = 1:N
        u = e*phi(m)/(kB*T);
        [~, rho_m3, eps_loc, c_plus_m3, c_minus_m3] = c_rho_eps_from_u(u);
        E(m)        = D(m) / max(eps_loc,1e-18);   
        rho_full(m) = rho_m3;
        c_plusL(m)  = c_plus_m3  / 1000;           % mol/L
        c_minusL(m) = c_minus_m3 / 1000;
    end

    % ===== bound water thickness =====
    eta = 0.005;                   
    E0  = abs(E(1));
    idx = find(abs(E) >= eta*E0, 1, 'last');
    if isempty(idx), idx = 1; end
    d_thick = x_fine(idx);

    % ===== Charge conservation (strict finite-domain formulation) =====
    sigma_int_full = trapz(x_fine, rho_full);     
    D_L            = D(end);
    delta_sigma    = sigma_int_full - (D_L - sigma_s);   
    delta_sigma_abs= abs(delta_sigma);                    

    % Normalized error: returned as err_record (scalar)
    err_record = delta_sigma_abs / max(abs(sigma_s), 1e-30);

    % ===== eps_r(0) (single algorithm retained) =====
    u0 = e*phi(1)/(kB*T);
    [c0_loc_m3, ~, ~, ~, ~] = c_rho_eps_from_u(u0);  % c(0) in mol/m^3
    eps_r0_final = eps_r_local_func(c0_loc_m3/1000); % convert to mol/L 

    % ===== Output matrix (consistent with the original interface) =====
    out_data = [x_fine(:), phi(:), c_plusL(:), c_minusL(:), E(:)];

    % ===== Logging (minimal; essential information only) =====
    fprintf('ss = %.2f, sm = %.2f, c0 = %.4g mol/L ****************\n', ss, sm, c0_molL);
    fprintf('eps_r(0) = %.3f\n', eps_r0_final);
    fprintf('delta_sigma_abs = %.3e C/m^2  (err = %.3e)\n', delta_sigma_abs, err_record);

    % ====================== ODE & BC (first-order mixed 耳每D formulation) ======================
    function dydx = edlODE(~,y)
        phi_loc = y(1);          
        D_loc   = y(2);          
        u = softclip(e*phi_loc/(kB*T), 80);
        [~, rho_m3, eps_loc, ~, ~] = c_rho_eps_from_u(u);
        dydx = [ - D_loc / max(eps_loc,1e-18) ;   
                  + rho_m3 ];                     
    end

    function res = edlBC(ya, yb)
        % bound: D(0) = 考sㄛ耳(L) = 0
        res = [ ya(2) - sigma_s ;
                yb(1) ];
    end

    % ====================== Utility functions ======================
    function us = softclip(u, umax)
        us = umax * tanh(u/umax);
    end

    % Given u, return: c, 老, 汍, and c㊣ 
    % (c and c㊣ in mol/m^3; 老 in C/m^3; 汍 in F/m)
    function [c_m3, rho_m3, eps_loc, c_plus_m3, c_minus_m3] = c_rho_eps_from_u(uq)
        A  = v_Na*NA*c0;                    % = v_+ N_A c0
        B  = v_Cl*NA*c0;                    % = v_- N_A c0
        eu = exp(+uq);                      % e^{+u}
        em = exp(-uq);                      % e^{-u}

        % Bikerman/FermiㄗNa+↙emㄛCl-↙euㄘ
        den = 1 + A*(em - 1) + B*(eu - 1);
        den = max(den, 1e-18);

        c_plus_m3  = c0 * em / den;         % c_+ = c0 e^{-u}/den
        c_minus_m3 = c0 * eu / den;         % c_- = c0 e^{+u}/den
        c_m3       = 0.5*(c_plus_m3 + c_minus_m3);

        rho_m3  = e*NA*(c_plus_m3 - c_minus_m3);
        eps_loc = eps0 * eps_r_local_func(c_m3/1000);   % c (mol/m^3) ↙ mol/L
    end
end

% ===== Gavish relative permittivity (c: mol/L) =====
function eps_r = eps_r_local_func(c)
    eps_w = 78.54; eps_ms = 30.08; alpha = 11.5;
    beta = eps_w - eps_ms;
    v = 3 * alpha * c / beta;
    L = coth(v) - 1./v;
    eps_r = eps_w - beta .* L;
end
