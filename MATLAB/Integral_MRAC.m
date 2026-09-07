function result = Integral_MRAC(params,command_type,controller_mode)
%INTEGRAL_MRAC Compare direct MRAC with and without integral action.
%   eI_dot = C*x-Cm*xm
%   xbar   = [x;eI]
%   u      = Kx_hat'*x + Ki_hat'*eI + Kr_hat'*r ...
%            - Theta_hat'*Phi(x).
%   Kbar_hat_dot = -Gamma_bar*xbar*ebar'*Pbar*Bbar,
% where Kbar_hat=[Kx_hat;Ki_hat] and ebar=[x-xm;eI].
% Baseline mode fixes Ki=0 but otherwise uses the same experiment.

    if nargin < 3
        error('Integral_MRAC:MissingMode', ...
            'Specify controller_mode as ''baseline'' or ''integral''.');
    end

    command_type = lower(char(command_type));
    controller_mode = lower(char(controller_mode));
    if ~ismember(command_type,{'step','multisine'})
        error('Integral_MRAC:UnknownCommandType', ...
            'command_type must be ''step'' or ''multisine''.');
    end
    if ~ismember(controller_mode,{'baseline','integral'})
        error('Integral_MRAC:UnknownMode', ...
            'controller_mode must be ''baseline'' or ''integral''.');
    end
    if ~isfield(params,'task3')
        error('Integral_MRAC:MissingConfiguration', ...
            'The params.task3 augmented-model configuration is required.');
    end

    z0 = [params.x0;
          params.xm0;
          params.task3.eI0;
          params.Kx_hat0(:);
          params.task3.Ki_hat0(:);
          params.Kr_hat0(:);
          params.Theta_hat0(:)];

    [t,z] = ode45(@closed_loop_rhs,[0 params.t_final],z0, ...
        params.ode_opts);

    if t(end) < params.t_final*(1-10*eps) || any(~isfinite(z),'all')
        error('Integral_MRAC:IntegrationFailure', ...
            'The ODE solver did not produce a finite trajectory over the requested interval.');
    end

    n_time = numel(t);
    x = z(:,1:6);
    xm = z(:,7:12);
    eI = z(:,13:15);
    e = x-xm;
    output_error = e*params.Cm.';

    r_history = zeros(n_time,3);
    u_history = zeros(n_time,3);
    Phi_history = zeros(n_time,6);
    eta_history = zeros(n_time,6);
    Kx_hat_history = zeros(n_time,6,3);
    Ki_hat_history = zeros(n_time,3,3);
    Kr_hat_history = zeros(n_time,3,3);
    Theta_hat_history = zeros(n_time,6,3);
    Kx_error_norm = zeros(n_time,1);
    Ki_error_norm = zeros(n_time,1);
    Kr_error_norm = zeros(n_time,1);
    Theta_error_norm = zeros(n_time,1);

    for k = 1:n_time
        [Kx_hat,Ki_hat,Kr_hat,Theta_hat] = ...
            unpack_estimates(z(k,:).');
        x_true = x(k,:).';
        xk = measured_state(t(k),x_true);
        eIk = eI(k,:).';
        rk = command_signal(t(k));
        Phi = mrac_regressor(xk);

        if strcmp(controller_mode,'integral')
            uk = Kx_hat.'*xk+Ki_hat.'*eIk ...
                + Kr_hat.'*rk-Theta_hat.'*Phi;
        else
            uk = Kx_hat.'*xk+Kr_hat.'*rk-Theta_hat.'*Phi;
        end

        [~,~,eta,Lambda_true,Theta_plant_T] = ...
            Plant_Dynamics(t(k),x_true,uk,params);
        [Kx_star,Ki_star,Kr_star,Theta_star] = ...
            ideal_parameters(Lambda_true,Theta_plant_T);

        r_history(k,:) = rk.';
        u_history(k,:) = uk.';
        Phi_history(k,:) = Phi.';
        eta_history(k,:) = eta.';
        Kx_hat_history(k,:,:) = Kx_hat;
        Ki_hat_history(k,:,:) = Ki_hat;
        Kr_hat_history(k,:,:) = Kr_hat;
        Theta_hat_history(k,:,:) = Theta_hat;
        Kx_error_norm(k) = norm(Kx_hat-Kx_star,'fro');
        Kr_error_norm(k) = norm(Kr_hat-Kr_star,'fro');
        Theta_error_norm(k) = norm(Theta_hat-Theta_star,'fro');
        if strcmp(controller_mode,'integral')
            Ki_error_norm(k) = norm(Ki_hat-Ki_star,'fro');
        else
            Ki_error_norm(k) = NaN;
        end
    end

    disturbance_mask = t >= params.disturbance.start_time;
    steady_start = max(params.disturbance.start_time, ...
        t(end)-params.task3.steady_state_window);
    steady_mask = t >= steady_start;

    result = struct();
    result.controller_mode = controller_mode;
    result.command_type = command_type;
    if strcmp(controller_mode,'integral')
        result.controller = 'Integral-augmented MRAC';
        result.label = 'Integral-augmented MRAC';
    else
        result.controller = 'Baseline direct MRAC';
        result.label = 'Baseline direct MRAC';
    end

    result.t = t;
    result.x = x;
    result.xm = xm;
    result.e = e;
    result.eI = eI;
    result.output_error = output_error;
    result.r = r_history;
    result.u = u_history;
    result.Phi = Phi_history;
    result.eta = eta_history;
    result.Kx_hat = Kx_hat_history;
    result.Ki_hat = Ki_hat_history;
    result.Kr_hat = Kr_hat_history;
    result.Theta_hat = Theta_hat_history;
    result.Kx_error_norm = Kx_error_norm;
    result.Ki_error_norm = Ki_error_norm;
    result.Kr_error_norm = Kr_error_norm;
    result.Theta_error_norm = Theta_error_norm;
    result.Kx_star = params.Kx_star;
    result.Ki_star = params.task3.Ki_star;
    result.Kr_star = params.Kr_star;
    result.Theta_star = params.Theta_star;
    result.disturbance_start_time = params.disturbance.start_time;

    result.rms_tracking_error = sqrt( ...
        trapz(t,sum(e.^2,2))/(t(end)-t(1)));
    result.rms_output_error = sqrt( ...
        trapz(t,sum(output_error.^2,2))/(t(end)-t(1)));
    result.post_disturbance_rms_output_error = rms_over_mask( ...
        t,output_error,disturbance_mask);
    result.steady_state_rms_output_error = rms_over_mask( ...
        t,output_error,steady_mask);
    result.maximum_control_effort = max(sqrt(sum(u_history.^2,2)));
    result.peak_post_disturbance_output_error = max( ...
        sqrt(sum(output_error(disturbance_mask,:).^2,2)));
    result.final_output_error = norm(output_error(end,:));
    result.final_integral_state = norm(eI(end,:));
    result.final_Ki_estimation_error = Ki_error_norm(end);
    result.final_Theta_estimation_error = Theta_error_norm(end);
    result.final_parameter_error = sqrt(Kx_error_norm(end)^2 ...
        +Kr_error_norm(end)^2+Theta_error_norm(end)^2 ...
        +max(0,Ki_error_norm(end))^2);

    function zdot = closed_loop_rhs(t_local,z_local)
        x_local = z_local(1:6);
        xm_local = z_local(7:12);
        eI_local = z_local(13:15);
        [Kx_hat_local,Ki_hat_local,Kr_hat_local,Theta_hat_local] = ...
            unpack_estimates(z_local);

        r_local = command_signal(t_local);
        x_measured = measured_state(t_local,x_local);
        Phi_local = mrac_regressor(x_measured);
        e_local = x_measured-xm_local;

        if strcmp(controller_mode,'integral')
            u_local = Kx_hat_local.'*x_measured ...
                + Ki_hat_local.'*eI_local ...
                + Kr_hat_local.'*r_local ...
                - Theta_hat_local.'*Phi_local;
        else
            u_local = Kx_hat_local.'*x_measured ...
                + Kr_hat_local.'*r_local ...
                - Theta_hat_local.'*Phi_local;
        end

        xdot_local = Plant_Dynamics(t_local,x_local,u_local,params);
        xmdot_local = params.Am*xm_local+params.Bm*r_local;
        eI_dot = params.Cm*x_measured-params.Cm*xm_local;

        if strcmp(controller_mode,'integral')
            ebar = [e_local;eI_local];
            xbar = [x_measured;eI_local];
            adaptation_signal = ...
                ebar.'*params.task3.Pbar*params.task3.Bbar;
            Kbar_hat_dot = -params.task3.Gamma_bar ...
                *xbar*adaptation_signal;
            Kx_hat_dot = Kbar_hat_dot(1:6,:);
            Ki_hat_dot = Kbar_hat_dot(7:9,:);
            Kr_hat_dot = -params.Gamma_r*r_local*adaptation_signal;
            Theta_hat_dot = params.Gamma_theta ...
                *Phi_local*adaptation_signal;
        else
            adaptation_signal = e_local.'*params.P*params.B;
            Kx_hat_dot = -params.Gamma_x ...
                *x_measured*adaptation_signal;
            Ki_hat_dot = zeros(3,3);
            Kr_hat_dot = -params.Gamma_r*r_local*adaptation_signal;
            Theta_hat_dot = params.Gamma_theta ...
                *Phi_local*adaptation_signal;
        end

        zdot = [xdot_local;
                xmdot_local;
                eI_dot;
                Kx_hat_dot(:);
                Ki_hat_dot(:);
                Kr_hat_dot(:);
                Theta_hat_dot(:)];
    end

    function [Kx_hat,Ki_hat,Kr_hat,Theta_hat] = ...
            unpack_estimates(z_local)
        index = 16;
        Kx_hat = reshape(z_local(index:index+17),6,3);
        index = index+18;
        Ki_hat = reshape(z_local(index:index+8),3,3);
        index = index+9;
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

    function [Kx_star,Ki_star,Kr_star,Theta_star] = ...
            ideal_parameters(Lambda_true,Theta_plant_T)
        Kx_star = (Lambda_true\ ...
            (params.G\[-params.Kp,-params.Kd])).';
        Ki_star = (Lambda_true\ ...
            (params.G\(-params.task3.Ki_des))).';
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
end

function value = rms_over_mask(t,signal,mask)
% Continuous-time vector RMS over a logical time mask.

    selected_time = t(mask);
    selected_signal = signal(mask,:);
    duration = selected_time(end)-selected_time(1);
    value = sqrt(trapz(selected_time, ...
        sum(selected_signal.^2,2))/duration);
end
