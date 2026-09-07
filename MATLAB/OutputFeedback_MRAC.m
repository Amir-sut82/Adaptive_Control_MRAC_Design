function result = OutputFeedback_MRAC(params,command_type,configuration)
%OUTPUTFEEDBACK_MRAC Simulate the Task-5 output-feedback controllers.
%   fullstate       - full-state benchmark used only as the LTR target
%   output_baseline - output feedback with the low-recovery observer
%   ltr             - output feedback with the selected LTR observer
%   augmented       - LTR output feedback plus Adaptive_Augmentation

    if nargin < 3
        error('OutputFeedback_MRAC:MissingConfiguration', ...
            'Specify the Task-5 controller configuration.');
    end

    command_type = lower(char(command_type));
    configuration = lower(char(configuration));
    if ~ismember(command_type,{'step','multisine'})
        error('OutputFeedback_MRAC:UnknownCommandType', ...
            'command_type must be ''step'' or ''multisine''.');
    end
    valid_configurations = {'fullstate','output_baseline', ...
        'ltr','augmented'};
    if ~ismember(configuration,valid_configurations)
        error('OutputFeedback_MRAC:UnknownConfiguration', ...
            'Unknown Task-5 configuration: %s.',configuration);
    end
    if ~isfield(params,'task5')
        error('OutputFeedback_MRAC:MissingTask5Configuration', ...
            'params.task5 is required.');
    end

    uses_true_state = strcmp(configuration,'fullstate');
    uses_augmentation = strcmp(configuration,'augmented');
    if strcmp(configuration,'output_baseline')
        observer_gain = params.task5.L_baseline;
    else
        observer_gain = params.task5.L_selected;
    end

    z0 = [params.x0;
          params.xm0;
          params.task5.xhat0;
          params.Kx_hat0(:);
          params.Kr_hat0(:);
          params.Theta_hat0(:);
          params.task5.W_hat0(:)];

    [t,z] = ode45(@closed_loop_rhs,[0 params.t_final],z0, ...
        params.ode_opts);
    if t(end) < params.t_final*(1-10*eps) || any(~isfinite(z),'all')
        error('OutputFeedback_MRAC:IntegrationFailure', ...
            'The Task-5 ODE solution is incomplete or nonfinite.');
    end

    n_time = numel(t);
    x = z(:,1:6);
    xm = z(:,7:12);
    xhat = z(:,13:18);
    e = x-xm;
    ehat_tracking = xhat-xm;
    observer_error = x-xhat;

    r_history = zeros(n_time,3);
    y_history = zeros(n_time,3);
    y_measured_history = zeros(n_time,3);
    measurement_noise_history = zeros(n_time,3);
    innovation_history = zeros(n_time,3);
    u_history = zeros(n_time,3);
    u_augmentation_history = zeros(n_time,3);
    Kx_hat_history = zeros(n_time,6,3);
    Kr_hat_history = zeros(n_time,3,3);
    Theta_hat_history = zeros(n_time,6,3);
    W_hat_history = zeros(n_time,4,3);
    parameter_error_norm = zeros(n_time,1);
    Theta_error_norm = zeros(n_time,1);
    W_error_norm = zeros(n_time,1);

    for time_index = 1:n_time
        [Kx_hat,Kr_hat,Theta_hat,W_hat] = ...
            unpack_estimates(z(time_index,:).');
        xk = x(time_index,:).';
        xhatk = xhat(time_index,:).';
        xmk = xm(time_index,:).';
        rk = command_signal(t(time_index));
        noise = measurement_noise(t(time_index));
        yk = params.Cm*xk;
        y_measured = yk+noise;
        innovation = y_measured-params.Cm*xhatk;

        if uses_true_state
            controller_state = xk;
        else
            controller_state = xhatk;
        end
        controller_error = controller_state-xmk;
        Phi_controller = mrac_regressor(controller_state);

        [u_augmentation,~,~] = Adaptive_Augmentation( ...
            controller_state,controller_error,W_hat,params, ...
            uses_augmentation);
        uk = Kx_hat.'*controller_state ...
            +Kr_hat.'*rk-Theta_hat.'*Phi_controller ...
            +u_augmentation;
        [~,~,~,Lambda_true,Theta_plant_T] = ...
            Plant_Dynamics(t(time_index),xk,uk,params);
        [Kx_star,Kr_star,Theta_star] = ...
            ideal_parameters(Lambda_true,Theta_plant_T);
        W_star = (Lambda_true\params.task5.W_true_T).';

        Kx_error = Kx_hat-Kx_star;
        Kr_error = Kr_hat-Kr_star;
        Theta_error = Theta_hat-Theta_star;
        W_error = W_hat-W_star;

        r_history(time_index,:) = rk.';
        y_history(time_index,:) = yk.';
        y_measured_history(time_index,:) = y_measured.';
        measurement_noise_history(time_index,:) = noise.';
        innovation_history(time_index,:) = innovation.';
        u_history(time_index,:) = uk.';
        u_augmentation_history(time_index,:) = u_augmentation.';
        Kx_hat_history(time_index,:,:) = Kx_hat;
        Kr_hat_history(time_index,:,:) = Kr_hat;
        Theta_hat_history(time_index,:,:) = Theta_hat;
        W_hat_history(time_index,:,:) = W_hat;
        parameter_error_norm(time_index) = sqrt( ...
            norm(Kx_error,'fro')^2+norm(Kr_error,'fro')^2 ...
            +norm(Theta_error,'fro')^2);
        Theta_error_norm(time_index) = norm(Theta_error,'fro');
        W_error_norm(time_index) = norm(W_error,'fro');
    end

    tail_start = max(t(1),t(end)-params.task5.tail_window);
    tail_mask = t >= tail_start;

    result = struct();
    result.configuration = configuration;
    result.controller = controller_label();
    result.command_type = command_type;
    result.label = command_label();
    result.uses_true_state = uses_true_state;
    result.uses_augmentation = uses_augmentation;
    result.observer_gain = observer_gain;
    result.t = t;
    result.x = x;
    result.xm = xm;
    result.xhat = xhat;
    result.e = e;
    result.ehat_tracking = ehat_tracking;
    result.observer_error = observer_error;
    result.r = r_history;
    result.y = y_history;
    result.y_measured = y_measured_history;
    result.measurement_noise = measurement_noise_history;
    result.innovation = innovation_history;
    result.u = u_history;
    result.u_augmentation = u_augmentation_history;
    result.Kx_hat = Kx_hat_history;
    result.Kr_hat = Kr_hat_history;
    result.Theta_hat = Theta_hat_history;
    result.W_hat = W_hat_history;
    result.parameter_error_norm = parameter_error_norm;
    result.Theta_error_norm = Theta_error_norm;
    result.W_error_norm = W_error_norm;
    result.Kx_star = params.Kx_star;
    result.Kr_star = params.Kr_star;
    result.Theta_star = params.Theta_star;
    result.W_star = params.task5.W_star;

    result.rms_tracking_error = vector_rms(t,e);
    result.tail_rms_tracking_error = vector_rms( ...
        t(tail_mask),e(tail_mask,:));
    result.rms_observer_error = vector_rms(t,observer_error);
    result.tail_rms_observer_error = vector_rms( ...
        t(tail_mask),observer_error(tail_mask,:));
    result.rms_innovation = vector_rms(t,innovation_history);
    result.maximum_control_effort = max(vecnorm(u_history,2,2));
    result.maximum_augmentation_effort = max( ...
        vecnorm(u_augmentation_history,2,2));
    result.final_parameter_error = parameter_error_norm(end);
    result.final_Theta_estimation_error = Theta_error_norm(end);
    if uses_augmentation
        result.final_augmentation_error = W_error_norm(end);
    else
        result.final_augmentation_error = NaN;
    end

    function zdot = closed_loop_rhs(t_local,z_local)
        x_local = z_local(1:6);
        xm_local = z_local(7:12);
        xhat_local = z_local(13:18);
        [Kx_hat_local,Kr_hat_local,Theta_hat_local,W_hat_local] = ...
            unpack_estimates(z_local);

        r_local = command_signal(t_local);
        noise_local = measurement_noise(t_local);
        y_measured_local = params.Cm*x_local+noise_local;
        innovation_local = ...
            y_measured_local-params.Cm*xhat_local;

        if uses_true_state
            controller_state_local = x_local;
        else
            controller_state_local = xhat_local;
        end
        controller_error_local = controller_state_local-xm_local;
        Phi_controller_local = ...
            mrac_regressor(controller_state_local);

        [u_augmentation_local,W_hat_dot,~] = ...
            Adaptive_Augmentation(controller_state_local, ...
                controller_error_local,W_hat_local,params, ...
                uses_augmentation);
        u_local = Kx_hat_local.'*controller_state_local ...
            +Kr_hat_local.'*r_local ...
            -Theta_hat_local.'*Phi_controller_local ...
            +u_augmentation_local;

        [xdot_local,~,~] = Plant_Dynamics( ...
            t_local,x_local,u_local,params);
        [~,~,psi_true] = Adaptive_Augmentation( ...
            x_local,zeros(6,1),zeros(4,3),params,false);
        xdot_local = xdot_local+params.B ...
            *(params.task5.W_true_T*psi_true);

        xmdot_local = params.Am*xm_local+params.Bm*r_local;
        xhat_dot = params.A*xhat_local+params.B*u_local ...
            +observer_gain*innovation_local;

        adaptation_signal = controller_error_local.' ...
            *params.P*params.B;
        Kx_hat_dot = -params.Gamma_x ...
            *controller_state_local*adaptation_signal;
        Kr_hat_dot = -params.Gamma_r ...
            *r_local*adaptation_signal;
        Theta_hat_dot = params.Gamma_theta ...
            *Phi_controller_local*adaptation_signal;

        zdot = [xdot_local;
                xmdot_local;
                xhat_dot;
                Kx_hat_dot(:);
                Kr_hat_dot(:);
                Theta_hat_dot(:);
                W_hat_dot(:)];
    end

    function [Kx_hat,Kr_hat,Theta_hat,W_hat] = ...
            unpack_estimates(z_local)
        index = 19;
        Kx_hat = reshape(z_local(index:index+17),6,3);
        index = index+18;
        Kr_hat = reshape(z_local(index:index+8),3,3);
        index = index+9;
        Theta_hat = reshape(z_local(index:index+17),6,3);
        index = index+18;
        W_hat = reshape(z_local(index:index+11),4,3);
    end

    function Phi = mrac_regressor(state)
        p = state(4);
        q = state(5);
        r_body = state(6);
        Phi = [q*r_body;
               p*r_body;
               p*q;
               p;
               q;
               r_body];
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
                tone_phase = params.command.multisine_frequency*tau ...
                    +params.command.multisine_phase;
                r_command = sum( ...
                    params.command.multisine_amplitude ...
                    .*sin(tone_phase),2);
        end
    end

    function noise = measurement_noise(t_local)
        if isfield(params,'task6') ...
                && isfield(params.task6,'output_noise')
            noise = params.task6.output_noise(t_local);
        else
            noise = params.task5.measurement_noise_amplitude ...
                .*sin(params.task5.measurement_noise_frequency*t_local ...
                    +params.task5.measurement_noise_phase);
        end
    end

    function [Kx_star,Kr_star,Theta_star] = ...
            ideal_parameters(Lambda_true,Theta_plant_T)
        Kx_star = (Lambda_true\ ...
            (params.G\[-params.Kp,-params.Kd])).';
        Kr_star = (Lambda_true\(params.G\params.Kp)).';
        Theta_star = (Lambda_true\Theta_plant_T).';
    end

    function label = controller_label()
        switch configuration
            case 'fullstate'
                label = 'Full-state MRAC benchmark';
            case 'output_baseline'
                label = 'Output-feedback MRAC';
            case 'ltr'
                label = 'Output-feedback MRAC + LTR';
            case 'augmented'
                label = 'LTR + adaptive augmentation';
        end
    end

    function label = command_label()
        if strcmp(command_type,'step')
            label = 'Step command';
        else
            label = 'Multisine command';
        end
    end
end

function value = vector_rms(t,signal)
% Continuous-time RMS of a vector signal.

    duration = t(end)-t(1);
    value = sqrt(trapz(t,sum(signal.^2,2))/duration);
end
