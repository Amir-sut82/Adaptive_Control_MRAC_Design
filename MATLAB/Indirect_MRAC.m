function result = Indirect_MRAC(params,command_type)
%INDIRECT_MRAC Simulate online identification and certainty equivalence.
% The series-parallel predictor estimates W from
%   W'*omega = Lambda*u + Theta_plant_T*Phi(x),
%   omega     = [u; Phi(x)].
%   x_predictor_dot = A*x + B*W_hat'*omega + L_id*epsilon
%   W_hat_dot       = Gamma_W*omega*epsilon'*P_id*B.
% Controller gains use a regularized inverse of Lambda_hat at each step.

    if nargin < 2
        error('Indirect_MRAC:MissingCommandType', ...
            'Specify command_type as ''step'' or ''multisine''.');
    end

    command_type = lower(char(command_type));
    if ~ismember(command_type,{'step','multisine'})
        error('Indirect_MRAC:UnknownCommandType', ...
            'command_type must be ''step'' or ''multisine''.');
    end
    if ~isfield(params,'indirect')
        error('Indirect_MRAC:MissingConfiguration', ...
            'The params.indirect estimator configuration is required.');
    end

    Kx_matching_T = params.G\[-params.Kp,-params.Kd];
    Kr_matching_T = params.G\params.Kp;
    W_hat0 = [params.indirect.Lambda_hat0.'; ...
              params.indirect.Theta_plant_hat_T0.'];

    z0 = [params.x0;
          params.xm0;
          params.indirect.x_predictor0;
          W_hat0(:)];

    [t,z] = ode45(@closed_loop_rhs,[0 params.t_final],z0, ...
        params.ode_opts);

    if t(end) < params.t_final*(1-10*eps) || any(~isfinite(z),'all')
        error('Indirect_MRAC:IntegrationFailure', ...
            'The ODE solver did not produce a finite trajectory over the requested interval.');
    end

    n_time = numel(t);
    x = z(:,1:6);
    xm = z(:,7:12);
    x_predictor = z(:,13:18);
    e = x-xm;
    prediction_error = zeros(n_time,6);

    r_history = zeros(n_time,3);
    u_history = zeros(n_time,3);
    Phi_history = zeros(n_time,6);
    omega_history = zeros(n_time,9);
    minimum_used_singular_value = zeros(n_time,1);

    Lambda_hat_history = zeros(n_time,3,3);
    Lambda_used_history = zeros(n_time,3,3);
    Theta_plant_hat_history = zeros(n_time,6,3);
    Theta_hat_history = zeros(n_time,6,3);
    Kx_hat_history = zeros(n_time,6,3);
    Kr_hat_history = zeros(n_time,3,3);

    Lambda_error_norm = zeros(n_time,1);
    Theta_plant_error_norm = zeros(n_time,1);
    Theta_error_norm = zeros(n_time,1);
    Kx_error_norm = zeros(n_time,1);
    Kr_error_norm = zeros(n_time,1);
    W_error_norm = zeros(n_time,1);

    for k = 1:n_time
        W_hat = reshape(z(k,19:end),9,3);
        [Lambda_hat,Theta_plant_hat_T] = unpack_identified_matrix(W_hat);
        [Lambda_inverse,Lambda_used,min_singular_value] = ...
            regularized_inverse(Lambda_hat);

        x_true = x(k,:).';
        xk = measured_state(t(k),x_true);
        rk = command_signal(t(k));
        Phi = mrac_regressor(xk);

        Kx_hat_T = Lambda_inverse*Kx_matching_T;
        Kr_hat_T = Lambda_inverse*Kr_matching_T;
        Theta_hat_T = Lambda_inverse*Theta_plant_hat_T;
        uk = Kx_hat_T*xk+Kr_hat_T*rk-Theta_hat_T*Phi;
        [~,~,~,Lambda_true,Theta_plant_T] = ...
            Plant_Dynamics(t(k),x_true,uk,params);
        [Kx_star,Kr_star,Theta_star] = ...
            ideal_parameters(Lambda_true,Theta_plant_T);
        W_star = [Lambda_true.';Theta_plant_T.'];

        r_history(k,:) = rk.';
        u_history(k,:) = uk.';
        Phi_history(k,:) = Phi.';
        omega_history(k,:) = [uk;Phi].';
        minimum_used_singular_value(k) = min_singular_value;

        Lambda_hat_history(k,:,:) = Lambda_hat;
        Lambda_used_history(k,:,:) = Lambda_used;
        Theta_plant_hat_history(k,:,:) = Theta_plant_hat_T.';
        Theta_hat_history(k,:,:) = Theta_hat_T.';
        Kx_hat_history(k,:,:) = Kx_hat_T.';
        Kr_hat_history(k,:,:) = Kr_hat_T.';
        prediction_error(k,:) = (xk-x_predictor(k,:).').';

        Lambda_error_norm(k) = norm( ...
            Lambda_hat-Lambda_true,'fro');
        Theta_plant_error_norm(k) = norm( ...
            Theta_plant_hat_T-Theta_plant_T,'fro');
        Theta_error_norm(k) = norm( ...
            Theta_hat_T.'-Theta_star,'fro');
        Kx_error_norm(k) = norm(Kx_hat_T.'-Kx_star,'fro');
        Kr_error_norm(k) = norm(Kr_hat_T.'-Kr_star,'fro');
        W_error_norm(k) = norm(W_hat-W_star,'fro');
    end

    [pe_min_eigenvalue_history,information_matrix] = ...
        excitation_diagnostics(t,omega_history,params.indirect.pe_window);
    information_singular_values = svd(information_matrix);
    if information_singular_values(1) == 0
        information_effective_rank = 0;
        information_condition_number = Inf;
    else
        information_effective_rank = sum(information_singular_values > ...
            1e-6*information_singular_values(1));
        if information_singular_values(end) <= eps( ...
                information_singular_values(1))
            information_condition_number = Inf;
        else
            information_condition_number = ...
                information_singular_values(1)/ ...
                information_singular_values(end);
        end
    end

    result = struct();
    result.label = command_label();
    result.command_type = command_type;
    result.controller = 'Indirect MRAC';
    result.t = t;
    result.x = x;
    result.xm = xm;
    result.x_predictor = x_predictor;
    result.e = e;
    result.prediction_error = prediction_error;
    result.r = r_history;
    result.u = u_history;
    result.Phi = Phi_history;
    result.omega = omega_history;

    result.Lambda_hat = Lambda_hat_history;
    result.Lambda_used = Lambda_used_history;
    result.Theta_plant_hat = Theta_plant_hat_history;
    result.Theta_hat = Theta_hat_history;
    result.Kx_hat = Kx_hat_history;
    result.Kr_hat = Kr_hat_history;

    result.Lambda_star = params.Lambda_true;
    result.Theta_plant_star = params.Theta_plant_T.';
    result.Theta_star = params.Theta_star;
    result.Kx_star = params.Kx_star;
    result.Kr_star = params.Kr_star;

    result.Lambda_error_norm = Lambda_error_norm;
    result.Theta_plant_error_norm = Theta_plant_error_norm;
    result.Theta_error_norm = Theta_error_norm;
    result.Kx_error_norm = Kx_error_norm;
    result.Kr_error_norm = Kr_error_norm;
    result.W_error_norm = W_error_norm;
    result.prediction_error_norm = sqrt(sum(prediction_error.^2,2));
    result.minimum_used_singular_value = minimum_used_singular_value;

    result.pe_window = params.indirect.pe_window;
    result.pe_min_eigenvalue = pe_min_eigenvalue_history;
    result.information_matrix = information_matrix;
    result.information_singular_values = information_singular_values;
    result.information_effective_rank = information_effective_rank;
    result.information_condition_number = information_condition_number;

    result.rms_tracking_error = sqrt( ...
        trapz(t,sum(e.^2,2))/(t(end)-t(1)));
    result.maximum_control_effort = max(sqrt(sum(u_history.^2,2)));
    result.final_tracking_error = norm(e(end,:));
    result.final_prediction_error = norm(prediction_error(end,:));
    result.final_Lambda_estimation_error = Lambda_error_norm(end);
    result.final_Theta_plant_estimation_error = ...
        Theta_plant_error_norm(end);
    result.final_Theta_estimation_error = Theta_error_norm(end);
    result.final_Kx_estimation_error = Kx_error_norm(end);
    result.final_Kr_estimation_error = Kr_error_norm(end);
    result.final_parameter_error = sqrt(Kx_error_norm(end)^2 ...
        +Kr_error_norm(end)^2+Theta_error_norm(end)^2);

    function zdot = closed_loop_rhs(t_local,z_local)
        x_local = z_local(1:6);
        xm_local = z_local(7:12);
        x_predictor_local = z_local(13:18);
        W_hat_local = reshape(z_local(19:end),9,3);

        [Lambda_hat_local,Theta_plant_hat_T_local] = ...
            unpack_identified_matrix(W_hat_local);
        Lambda_inverse_local = regularized_inverse(Lambda_hat_local);

        r_local = command_signal(t_local);
        x_measured = measured_state(t_local,x_local);
        Phi_local = mrac_regressor(x_measured);

        Kx_hat_T_local = Lambda_inverse_local*Kx_matching_T;
        Kr_hat_T_local = Lambda_inverse_local*Kr_matching_T;
        Theta_hat_T_local = ...
            Lambda_inverse_local*Theta_plant_hat_T_local;

        u_local = Kx_hat_T_local*x_measured ...
            + Kr_hat_T_local*r_local ...
            - Theta_hat_T_local*Phi_local;

        xdot_local = Plant_Dynamics(t_local,x_local,u_local,params);
        xmdot_local = params.Am*xm_local+params.Bm*r_local;

        omega_local = [u_local;Phi_local];
        prediction_error_local = x_measured-x_predictor_local;
        x_predictor_dot = params.A*x_measured ...
            + params.B*(W_hat_local.'*omega_local) ...
            + params.indirect.L_id*prediction_error_local;

        W_hat_dot = params.indirect.Gamma_W*omega_local ...
            *(prediction_error_local.'*params.indirect.P_id*params.B);

        zdot = [xdot_local;
                xmdot_local;
                x_predictor_dot;
                W_hat_dot(:)];
    end

    function [Lambda_hat,Theta_plant_hat_T] = ...
            unpack_identified_matrix(W_hat)
        Lambda_hat = W_hat(1:3,:).';
        Theta_plant_hat_T = W_hat(4:9,:).';
    end

    function [Lambda_inverse,Lambda_used,min_singular_value] = ...
            regularized_inverse(Lambda_hat)
        [U,S,V] = svd(Lambda_hat);
        singular_values = diag(S);
        regularized_values = max(singular_values, ...
            params.indirect.lambda_singular_value_floor);
        Lambda_inverse = V*diag(1./regularized_values)*U.';
        Lambda_used = U*diag(regularized_values)*V.';
        min_singular_value = min(regularized_values);
    end

    function Phi = mrac_regressor(x_local)
        p = x_local(4);
        q = x_local(5);
        r_body = x_local(6);
        Phi = [q*r_body;
               p*r_body;
               p*q;
               p;
               q;
               r_body];
    end

    function x_measured = measured_state(t_local,x_true)
        x_measured = x_true;
        if isfield(params,'task6') ...
                && isfield(params.task6,'state_noise')
            x_measured = x_true+params.task6.state_noise(t_local);
        end
    end

    function [Kx_star,Kr_star,Theta_star] = ...
            ideal_parameters(Lambda_true,Theta_plant_T)
        Kx_star = (Lambda_true\Kx_matching_T).';
        Kr_star = (Lambda_true\Kr_matching_T).';
        Theta_star = (Lambda_true\Theta_plant_T).';
    end

    function r_command = command_signal(t_local)
        if t_local < params.command.start_time
            r_command = zeros(3,1);
            return;
        end

        switch command_type
            case 'step'
                r_command = params.command.step_value;
            case 'multisine'
                tau = t_local-params.command.start_time;
                tone_phase = ...
                    params.command.multisine_frequency*tau ...
                    + params.command.multisine_phase;
                r_command = sum( ...
                    params.command.multisine_amplitude.*sin(tone_phase),2);
        end
    end

    function label = command_label()
        switch command_type
            case 'step'
                label = 'Step command';
            case 'multisine'
                label = 'Multisine command';
        end
    end
end

function [minimum_eigenvalue_history,information_matrix] = ...
        excitation_diagnostics(t,omega_history,window_length)
% Compute a sliding-window finite-time excitation indicator.

    n_time = numel(t);
    n_regressor = size(omega_history,2);
    cumulative = zeros(n_regressor,n_regressor,n_time);

    for k = 2:n_time
        previous_outer = omega_history(k-1,:).'*omega_history(k-1,:);
        current_outer = omega_history(k,:).'*omega_history(k,:);
        cumulative(:,:,k) = cumulative(:,:,k-1) ...
            + 0.5*(t(k)-t(k-1))*(previous_outer+current_outer);
    end

    minimum_eigenvalue_history = zeros(n_time,1);
    window_start = 1;
    for k = 2:n_time
        while window_start < k ...
                && t(k)-t(window_start) > window_length
            window_start = window_start+1;
        end

        duration = t(k)-t(window_start);
        if duration <= 0
            continue;
        end

        window_gramian = (cumulative(:,:,k)- ...
            cumulative(:,:,window_start))/duration;
        window_gramian = (window_gramian+window_gramian.')/2;
        minimum_eigenvalue_history(k) = max(0,min(eig(window_gramian)));
    end

    total_duration = t(end)-t(1);
    information_matrix = cumulative(:,:,end)/total_duration;
    information_matrix = (information_matrix+information_matrix.')/2;
end
