function design = Observer_LTR(params,rho_values,frequency)
%OBSERVER_LTR Design observers along a loop-transfer-recovery sweep.
%   Q_o(rho) = Q_o0 + (1/rho)*B*B'
% Then L=P*C'/R. Smaller rho drives the input return ratio toward the
% full-state design, at the cost of a higher observer gain.

    if nargin < 2 || isempty(rho_values)
        rho_values = [1e2,1,1e-2,1e-4];
    end
    if nargin < 3 || isempty(frequency)
        frequency = logspace(-2,2.5,260);
    end

    rho_values = rho_values(:).';
    frequency = frequency(:);
    if any(~isfinite(rho_values)) || any(rho_values <= 0)
        error('Observer_LTR:InvalidRecoveryParameter', ...
            'Every recovery parameter must be positive and finite.');
    end
    if any(~isfinite(frequency)) || any(frequency <= 0)
        error('Observer_LTR:InvalidFrequencyGrid', ...
            'The frequency grid must be positive and finite.');
    end

    A = params.A;
    B = params.B;
    C = params.Cm;
    n = size(A,1);
    p = size(C,1);
    m = size(B,2);
    if size(C,2) ~= n || p >= n
        error('Observer_LTR:OutputDimension', ...
            'Task 5 requires C in R^(p-by-n) with p<n.');
    end

    % Nominal full-state feedback is u=-Ksf*x.
    Ksf = params.G\[params.Kp,params.Kd];
    observer_noise_covariance = eye(p);
    process_noise_floor = 1e-6*eye(n);

    n_rho = numel(rho_values);
    n_frequency = numel(frequency);
    L_history = zeros(n,p,n_rho);
    observer_spectral_abscissa = zeros(n_rho,1);
    observer_gain_norm = zeros(n_rho,1);
    riccati_residual = zeros(n_rho,1);
    full_state_loop_sv = zeros(n_frequency,m);
    observer_loop_sv = zeros(n_frequency,m,n_rho);
    recovery_log_sv_rmse = zeros(n_rho,1);

    identity_state = eye(n);
    for frequency_index = 1:n_frequency
        s = 1i*frequency(frequency_index);
        full_state_loop = Ksf*((s*identity_state-A)\B);
        full_state_loop_sv(frequency_index,:) = ...
            svd(full_state_loop).';
    end

    for rho_index = 1:n_rho
        rho = rho_values(rho_index);
        process_noise_covariance = process_noise_floor ...
            +(1/rho)*(B*B.');
        P_observer = care(A.',C.',process_noise_covariance, ...
            observer_noise_covariance);
        L = P_observer*C.'/observer_noise_covariance;
        L_history(:,:,rho_index) = L;

        observer_matrix = A-L*C;
        observer_spectral_abscissa(rho_index) = ...
            max(real(eig(observer_matrix)));
        observer_gain_norm(rho_index) = norm(L,2);
        riccati_residual(rho_index) = norm( ...
            A*P_observer+P_observer*A.' ...
            -P_observer*C.' ...
                *(observer_noise_covariance\(C*P_observer)) ...
            +process_noise_covariance,'fro');

        controller_matrix = A-B*Ksf-L*C;
        for frequency_index = 1:n_frequency
            s = 1i*frequency(frequency_index);
            plant_transfer = C*((s*identity_state-A)\B);
            observer_controller = Ksf ...
                *((s*identity_state-controller_matrix)\L);
            observer_loop = observer_controller*plant_transfer;
            observer_loop_sv(frequency_index,:,rho_index) = ...
                svd(observer_loop).';
        end

        log_difference = log10(max( ...
            observer_loop_sv(:,:,rho_index),realmin)) ...
            -log10(max(full_state_loop_sv,realmin));
        recovery_log_sv_rmse(rho_index) = sqrt( ...
            mean(log_difference.^2,'all'));
    end

    design = struct();
    design.C = C;
    design.Ksf = Ksf;
    design.rho = rho_values;
    design.frequency = frequency;
    design.L = L_history;
    design.observer_spectral_abscissa = ...
        observer_spectral_abscissa;
    design.observer_gain_norm = observer_gain_norm;
    design.riccati_residual = riccati_residual;
    design.full_state_loop_sv = full_state_loop_sv;
    design.observer_loop_sv = observer_loop_sv;
    design.recovery_log_sv_rmse = recovery_log_sv_rmse;
end
