function [xdot,Phi,eta,Lambda_true,Theta_plant_T] = ...
        Plant_Dynamics(t,x,u,params)
%PLANT_DYNAMICS Evaluate the nonlinear attitude plant.
% xdot=A*x+B*(Lambda*u+Theta_plant_T*Phi(x))+eta(x,t), where
% x=[phi;theta;psi;p;q;r] and Phi=[q*r;p*r;p*q;p;q;r].

    x = x(:);
    u = u(:);

    if numel(x) ~= 6
        error('Plant_Dynamics:StateDimension', ...
            'x must be a 6-by-1 vector.');
    end
    if numel(u) ~= 3
        error('Plant_Dynamics:InputDimension', ...
            'u must contain the three active inputs [u2;u3;u4].');
    end
    if any(~isfinite([x;u]))
        error('Plant_Dynamics:NonFiniteInput', ...
            'The state and control input must be finite.');
    end

    p = x(4);
    q = x(5);
    r = x(6);

    Phi = [q*r;
           p*r;
           p*q;
           p;
           q;
           r];

    eta = bounded_disturbance(t,params);
    [Lambda_true,Theta_plant_T] = scheduled_parameters(t,params);

    xdot = params.A*x + params.B*( ...
        Lambda_true*u + Theta_plant_T*Phi) + eta;
end

function [Lambda_true,Theta_plant_T] = scheduled_parameters(t,params)
% Apply the Task 6 payload change after its scheduled time.

    Lambda_true = params.Lambda_true;
    Theta_plant_T = params.Theta_plant_T;

    if ~isfield(params,'task6') ...
            || ~isfield(params.task6,'parameter_step')
        return;
    end

    step = params.task6.parameter_step;
    if step.enabled && t >= step.time
        Lambda_true = step.Lambda_post;
        Theta_plant_T = step.Theta_plant_T_post;
    end
end

function eta = bounded_disturbance(t,params)
% Return the disturbance selected for the current experiment.

    if isfield(params.disturbance,'start_time') ...
            && t < params.disturbance.start_time
        eta = zeros(6,1);
        return;
    end

    switch lower(params.disturbance.type)
        case 'none'
            eta = zeros(6,1);

        case 'constant'
            eta = params.disturbance.constant(:);

        case 'sine'
            eta = params.disturbance.amplitude(:).*sin( ...
                params.disturbance.frequency(:)*t);

        case 'multisine'
            amplitude = params.disturbance.amplitude;
            frequency = params.disturbance.frequency;
            if ~isequal(size(amplitude),size(frequency)) ...
                    || size(amplitude,1) ~= 6
                error('Plant_Dynamics:MultisineDimension', ...
                    ['Multisine disturbance amplitude and frequency ' ...
                     'must be equally sized 6-by-N matrices.']);
            end
            if isfield(params.disturbance,'phase')
                phase = params.disturbance.phase;
            else
                phase = zeros(size(amplitude));
            end
            if ~isequal(size(phase),size(amplitude))
                error('Plant_Dynamics:MultisinePhaseDimension', ...
                    'Multisine disturbance phase must match amplitude.');
            end
            eta = sum(amplitude.*sin(frequency*t+phase),2);

        otherwise
            error('Plant_Dynamics:UnknownDisturbance', ...
                'Unknown disturbance type: %s',params.disturbance.type);
    end

    if numel(eta) ~= 6 || any(~isfinite(eta))
        error('Plant_Dynamics:DisturbanceDimension', ...
            'eta(x,t) must be a finite 6-by-1 vector.');
    end
end
