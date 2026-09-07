function result = Direct_MRAC(params,command_type,robust_method)
%DIRECT_MRAC Simulate the direct controller and its robust variants.
%   u = Kx_hat'*x + Kr_hat'*r - Theta_hat'*Phi(x)
%   Kx_hat_dot    = -Gamma_x*x*e'*P*B
%   Kr_hat_dot    = -Gamma_r*r*e'*P*B
%   Theta_hat_dot =  Gamma_theta*Phi(x)*e'*P*B
% robust_method selects baseline, deadzone, sigma, emod, or projection.

    if nargin < 2
        error('Direct_MRAC:MissingCommandType', ...
            'Specify command_type as ''step'' or ''multisine''.');
    end
    if nargin < 3
        robust_method = 'nominal';
    end

    command_type = lower(char(command_type));
    robust_method = lower(char(robust_method));
    if ~ismember(command_type,{'step','multisine'})
        error('Direct_MRAC:UnknownCommandType', ...
            'command_type must be ''step'' or ''multisine''.');
    end
    valid_methods = {'nominal','baseline','deadzone', ...
        'sigma','emod','projection'};
    if ~ismember(robust_method,valid_methods)
        error('Direct_MRAC:UnknownRobustMethod', ...
            'Unknown robust adaptive-law selection: %s.',robust_method);
    end

    if strcmp(robust_method,'nominal') ...
            && ~strcmpi(params.disturbance.type,'none')
        error('Direct_MRAC:NonNominalDisturbance', ...
            'Task 1 requires eta(x,t)=0. Set disturbance.type to ''none''.');
    end
    if ~strcmp(robust_method,'nominal') && ~isfield(params,'task4')
        error('Direct_MRAC:MissingTask4Configuration', ...
            'Robust simulations require params.task4.');
    end

    z0 = [params.x0;
          params.xm0;
          params.Kx_hat0(:);
          params.Kr_hat0(:);
          params.Theta_hat0(:)];

    [t,z] = ode45(@closed_loop_rhs,[0 params.t_final],z0, ...
        params.ode_opts);

    n_time = numel(t);
    x  = z(:,1:6);
    xm = z(:,7:12);
    e  = x-xm;

    r_history = zeros(n_time,3);
    u_history = zeros(n_time,3);
    Phi_history = zeros(n_time,6);
    eta_history = zeros(n_time,6);
    V_history = zeros(n_time,1);
    Vdot_theory = zeros(n_time,1);
    Kx_error_norm = zeros(n_time,1);
    Kr_error_norm = zeros(n_time,1);
    Theta_error_norm = zeros(n_time,1);
    parameter_estimate_norm = zeros(n_time,1);
    parameter_error_norm = zeros(n_time,1);
    adaptation_active = true(n_time,1);
    projection_active_count = zeros(n_time,1);

    Kx_hat_history = zeros(n_time,6,3);
    Kr_hat_history = zeros(n_time,3,3);
    Theta_hat_history = zeros(n_time,6,3);

    for k = 1:n_time
        [Kx_hat,Kr_hat,Theta_hat] = unpack_estimates(z(k,:).');

        x_true = x(k,:).';
        xk = measured_state(t(k),x_true);
        ek = xk-xm(k,:).';
        rk = command_signal(t(k));
        Phi = mrac_regressor(xk);

        uk = Kx_hat.'*xk + Kr_hat.'*rk - Theta_hat.'*Phi;
        [~,~,eta,Lambda_true,Theta_plant_T] = ...
            Plant_Dynamics(t(k),x_true,uk,params);
        [Kx_star,Kr_star,Theta_star] = ...
            ideal_parameters(Lambda_true,Theta_plant_T);

        Kx_tilde = Kx_hat-Kx_star;
        Kr_tilde = Kr_hat-Kr_star;
        Theta_tilde = Theta_hat-Theta_star;

        r_history(k,:) = rk.';
        u_history(k,:) = uk.';
        Phi_history(k,:) = Phi.';
        eta_history(k,:) = eta.';

        Kx_hat_history(k,:,:) = Kx_hat;
        Kr_hat_history(k,:,:) = Kr_hat;
        Theta_hat_history(k,:,:) = Theta_hat;

        V_history(k) = ek.'*params.P*ek ...
            + trace(Kx_tilde.'*(params.Gamma_x\Kx_tilde) ...
                *Lambda_true) ...
            + trace(Kr_tilde.'*(params.Gamma_r\Kr_tilde) ...
                *Lambda_true) ...
            + trace(Theta_tilde.'*(params.Gamma_theta\Theta_tilde) ...
                *Lambda_true);

        Vdot_theory(k) = -ek.'*params.Q*ek;
        Kx_error_norm(k) = norm(Kx_tilde,'fro');
        Kr_error_norm(k) = norm(Kr_tilde,'fro');
        Theta_error_norm(k) = norm(Theta_tilde,'fro');
        parameter_estimate_norm(k) = sqrt( ...
            norm(Kx_hat,'fro')^2+norm(Kr_hat,'fro')^2 ...
            +norm(Theta_hat,'fro')^2);
        parameter_error_norm(k) = sqrt( ...
            Kx_error_norm(k)^2+Kr_error_norm(k)^2 ...
            +Theta_error_norm(k)^2);

        if strcmp(robust_method,'deadzone')
            adaptation_active(k) = norm(ek,2) ...
                > params.task4.deadzone_epsilon;
        elseif strcmp(robust_method,'projection')
            projection_active_count(k) = ...
                count_projected_entries(Kx_hat,Kr_hat,Theta_hat, ...
                    xk,rk,Phi,ek);
        end
    end

    if t(end) < params.t_final*(1-10*eps) || any(~isfinite(z),'all')
        error('Direct_MRAC:IntegrationFailure', ...
            'The ODE solver did not produce a finite trajectory over the requested interval.');
    end

    result = struct();
    result.label = command_label();
    result.command_type = command_type;
    result.controller = controller_label();
    result.robust_method = robust_method;
    result.t = t;
    result.x = x;
    result.xm = xm;
    result.e = e;
    result.r = r_history;
    result.u = u_history;
    result.Phi = Phi_history;
    result.eta = eta_history;
    result.V = V_history;
    result.Vdot_theory = Vdot_theory;
    result.Kx_error_norm = Kx_error_norm;
    result.Kr_error_norm = Kr_error_norm;
    result.Theta_error_norm = Theta_error_norm;
    result.parameter_estimate_norm = parameter_estimate_norm;
    result.parameter_error_norm = parameter_error_norm;
    result.adaptation_active = adaptation_active;
    result.projection_active_count = projection_active_count;

    result.Kx_hat = Kx_hat_history;
    result.Kr_hat = Kr_hat_history;
    result.Theta_hat = Theta_hat_history;

    result.Kx_star = params.Kx_star;
    result.Kr_star = params.Kr_star;
    result.Theta_star = params.Theta_star;

    result.rms_tracking_error = sqrt( ...
        trapz(t,sum(e.^2,2))/(t(end)-t(1)));
    result.maximum_control_effort = max(sqrt(sum(u_history.^2,2)));
    result.final_tracking_error = norm(e(end,:));
    result.final_Kx_estimation_error = Kx_error_norm(end);
    result.final_Kr_estimation_error = Kr_error_norm(end);
    result.final_Theta_estimation_error = Theta_error_norm(end);
    if isfield(params,'task4') && isfield(params.task4,'tail_window')
        tail_window = params.task4.tail_window;
    else
        tail_window = min(10,0.25*(t(end)-t(1)));
    end
    tail_mask = t >= t(end)-tail_window;
    tail_time = t(tail_mask);
    tail_error = e(tail_mask,:);
    result.tail_rms_tracking_error = sqrt(trapz(tail_time, ...
        sum(tail_error.^2,2))/(tail_time(end)-tail_time(1)));
    result.peak_tracking_error = max(vecnorm(e,2,2));
    result.final_parameter_error = parameter_error_norm(end);
    result.maximum_parameter_error = max(parameter_error_norm);
    result.maximum_parameter_estimate_norm = ...
        max(parameter_estimate_norm);
    initial_parameter_vector = [params.Kx_hat0(:); ...
        params.Kr_hat0(:);params.Theta_hat0(:)];
    final_parameter_vector = [reshape(Kx_hat_history(end,:,:),[],1); ...
        reshape(Kr_hat_history(end,:,:),[],1); ...
        reshape(Theta_hat_history(end,:,:),[],1)];
    result.parameter_excursion_from_initial = norm( ...
        final_parameter_vector-initial_parameter_vector,2);
    result.adaptation_duty_fraction = mean(adaptation_active);
    result.projection_active_fraction = mean( ...
        projection_active_count > 0);

    function zdot = closed_loop_rhs(t_local,z_local)
        x_local = z_local(1:6);
        xm_local = z_local(7:12);

        [Kx_hat_local,Kr_hat_local,Theta_hat_local] = ...
            unpack_estimates(z_local);

        r_local = command_signal(t_local);
        x_measured = measured_state(t_local,x_local);
        e_local = x_measured-xm_local;
        Phi_local = mrac_regressor(x_measured);

        u_local = Kx_hat_local.'*x_measured ...
            + Kr_hat_local.'*r_local ...
            - Theta_hat_local.'*Phi_local;

        xdot_local = Plant_Dynamics(t_local,x_local,u_local,params);
        xmdot_local = params.Am*xm_local+params.Bm*r_local;

        adaptation_signal = e_local.'*params.P*params.B;

        raw_Kx_dot = -params.Gamma_x*x_measured*adaptation_signal;
        raw_Kr_dot = -params.Gamma_r*r_local*adaptation_signal;
        raw_Theta_dot = params.Gamma_theta ...
            *Phi_local*adaptation_signal;

        switch robust_method
            case {'nominal','baseline'}
                Kx_hat_dot = raw_Kx_dot;
                Kr_hat_dot = raw_Kr_dot;
                Theta_hat_dot = raw_Theta_dot;

            case 'deadzone'
                [Kx_hat_dot,Kr_hat_dot,Theta_hat_dot] = ...
                    Deadzone_MRAC(raw_Kx_dot,raw_Kr_dot, ...
                        raw_Theta_dot,e_local, ...
                        params.task4.deadzone_epsilon);

            case 'sigma'
                [Kx_hat_dot,Kr_hat_dot,Theta_hat_dot] = ...
                    Sigma_MRAC(raw_Kx_dot,raw_Kr_dot, ...
                        raw_Theta_dot,Kx_hat_local, ...
                        Kr_hat_local,Theta_hat_local,params);

            case 'emod'
                [Kx_hat_dot,Kr_hat_dot,Theta_hat_dot] = ...
                    eMod_MRAC(raw_Kx_dot,raw_Kr_dot, ...
                        raw_Theta_dot,Kx_hat_local, ...
                        Kr_hat_local,Theta_hat_local,e_local,params);

            case 'projection'
                Kx_hat_dot = Projection_Operator(Kx_hat_local, ...
                    raw_Kx_dot,params.task4.projection.Kx_lower, ...
                    params.task4.projection.Kx_upper);
                Kr_hat_dot = Projection_Operator(Kr_hat_local, ...
                    raw_Kr_dot,params.task4.projection.Kr_lower, ...
                    params.task4.projection.Kr_upper);
                Theta_hat_dot = Projection_Operator(Theta_hat_local, ...
                    raw_Theta_dot, ...
                    params.task4.projection.Theta_lower, ...
                    params.task4.projection.Theta_upper);
        end

        zdot = [xdot_local;
                xmdot_local;
                Kx_hat_dot(:);
                Kr_hat_dot(:);
                Theta_hat_dot(:)];
    end

    function [Kx_hat,Kr_hat,Theta_hat] = unpack_estimates(z_local)
        index = 13;

        Kx_hat = reshape(z_local(index:index+17),6,3);
        index = index+18;

        Kr_hat = reshape(z_local(index:index+8),3,3);
        index = index+9;

        Theta_hat = reshape(z_local(index:index+17),6,3);
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
        Kx_star = (Lambda_true\ ...
            (params.G\[-params.Kp,-params.Kd])).';
        Kr_star = (Lambda_true\(params.G\params.Kp)).';
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

    function label = controller_label()
        switch robust_method
            case 'nominal'
                label = 'Direct MRAC';
            case 'baseline'
                label = 'Baseline MRAC';
            case 'deadzone'
                label = 'Dead-zone MRAC';
            case 'sigma'
                label = 'Sigma-modification';
            case 'emod'
                label = 'e-modification';
            case 'projection'
                label = 'Projection MRAC';
        end
    end

    function active_count = count_projected_entries( ...
            Kx_hat,Kr_hat,Theta_hat,x_local,r_local,Phi_local,e_local)
        adaptation_signal_local = e_local.'*params.P*params.B;
        raw_Kx_dot_local = -params.Gamma_x ...
            *x_local*adaptation_signal_local;
        raw_Kr_dot_local = -params.Gamma_r ...
            *r_local*adaptation_signal_local;
        raw_Theta_dot_local = params.Gamma_theta ...
            *Phi_local*adaptation_signal_local;
        [~,Kx_count] = Projection_Operator(Kx_hat, ...
            raw_Kx_dot_local,params.task4.projection.Kx_lower, ...
            params.task4.projection.Kx_upper);
        [~,Kr_count] = Projection_Operator(Kr_hat, ...
            raw_Kr_dot_local,params.task4.projection.Kr_lower, ...
            params.task4.projection.Kr_upper);
        [~,Theta_count] = Projection_Operator(Theta_hat, ...
            raw_Theta_dot_local, ...
            params.task4.projection.Theta_lower, ...
            params.task4.projection.Theta_upper);
        active_count = Kx_count+Kr_count+Theta_count;
    end
end
