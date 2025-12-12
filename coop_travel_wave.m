function coop11_18()
    num_masses = 250;
    total_mass = 1;
    tension_force = 3;
    string_length = 10;
    damping_coeff = 0.001;
    dx = string_length/(num_masses+1);
    amplitude_Uf = 0.1;
    % omega_Uf = pi/3;

    %generate the struct
    string_params = struct();
    string_params.n = num_masses;
    string_params.M = total_mass;
    string_params.Tf = tension_force;
    string_params.L = string_length;
    string_params.c = damping_coeff;
    string_params.dx = dx;
    string_params.wave = sqrt(tension_force/(total_mass/string_length));

    %modal analysis
    [M_mat,K_mat] = construct_2nd_order_matrices(string_params);
    %Use MATLAB to solve the generalized eigenvalue problem
    [Ur_mat,lambda_mat] = eig(K_mat,M_mat);
    mode_num = 1;
    omega_n = sqrt(lambda_mat(mode_num,mode_num));

    %list of x points (including the two endpoints)
    xlist = linspace(0,string_length,num_masses+2);
    w = 0.25;
    h = 1;
    Uf_func = @(t_in) triangle_pulse(t_in,w,h);
    dUfdt_func = @(t_in) triangle_pulse_derivative(t_in,w,h);
    string_params.Uf_func = Uf_func;
    string_params.dUfdt_func = dUfdt_func;

    %load string_params into rate function
    my_rate_func = @(t_in,V_in) string_rate_func01(t_in,V_in,string_params);
    %initial conditions
    U0 = zeros(string_params.n, 1);
    dUdt0 = zeros(string_params.n, 1);
    V0 = [U0;dUdt0];
    tspan = [0,4*string_length/string_params.wave];
    %run the integration
    DP = make_DP_tableau();
    h_ref = 1e-5;
    err_des = 1e-12;
    [t_list,V_list,~,~,~,~] = rk_variable(my_rate_func,tspan,V0,h_ref,DP,5,err_des);

    %generate an animation of the system
    filename = "testavi.avi";
    
    record_animation_avi(filename, 30, t_list, V_list, string_params,0.25,w)
end

%INPUTS
function dVdt = string_rate_func01(t,V,string_params)
    n = string_params.n; %number of masses
    m = string_params.M; %total mass attached to the string
    Uf_func = string_params.Uf_func; %function describing motion of end point
    dUfdt_func = string_params.dUfdt_func; %time derivative of Uf
    Tf = string_params.Tf; %tension in string
    L = string_params.L; %length of string
    c = string_params.c; %damping coefficient
    dx = string_params.dx; %horizontal spacing between masses
    %unpack state variable
    U = V(1:n);
    dUdt = V((n+1):(2*n));
    Uf = Uf_func(t);
    dUfdt = dUfdt_func(t);
    %compute acceleration
    m_per_ball = m/n;
    
    %matrices:
    % % d2Udt2 = (Tf/dx*K*U+IC)\M;
    % Uminus = circshift(U,[-1,0]);
    % Uminus(end) = 0;
    Uminus = [U(2:end);0];
    Uplus = [0;U(1:end-1)];

    dUdtminus = [dUdt(2:end);0];
    dUdtplus = [0;dUdt(1:end-1)];
    
    % Uplus = circshift(U,[-1,0]);
    % Uplus(1) = 0;
    ForcingA = [zeros(n-1,1);Uf];
    UtotA = -2*U+Uminus+Uplus+ForcingA;

    ForcingB = [zeros(n-1,1);dUfdt];
    UtotB = -2*dUdt+dUdtminus+dUdtplus+ForcingB;
    
    d2Udt2 = (Tf/dx*UtotA+c/dx*UtotB)/m_per_ball;
    %assemble state derivative
    dVdt = [dUdt;d2Udt2];
end

%make M and K
function [M_mat,K_mat] = construct_2nd_order_matrices(string_params)
    n = string_params.n;
    num_per_ball = string_params.M/n;
    M_mat = num_per_ball*eye(n);
    Q = -2*eye(n)+[zeros(n-1,1),eye(n-1);zeros(1,n);]+[zeros(1,n);eye(n-1),zeros(n-1,1);];
    K_mat = -string_params.Tf/string_params.dx*Q;
end

%%%NUMERICAL INTEGRATION
function DP = make_DP_tableau()
    DP.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
    DP.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0; ...
            5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
    DP.A = [0,0,0,0,0,0,0; ...
            1/5,0,0,0,0,0,0; ...
            3/40,9/40,0,0,0,0,0; ...
            44/45,-56/15,32/9,0,0,0,0; ...
            19372/6561,-25360/2187,64448/6561,-212/729,0,0,0; ...
            9017/3168,-355/33,46732/5247,49/176,-5103/18656,0,0; ...
            35/384,0,500/1113,125/192,-2187/6784,11/84,0];
end

% embedded RK step for adaptive DP
function [XB1, XB2, num_evals] = rk_step_embedded(rate,t,XA,h,BT)
    A = BT.A;
    B = BT.B;
    C = BT.C;
    s = numel(C);
    n = numel(XA);
    K = zeros(n,s);
    for i = 1:s
        a = A(i,1:i-1);
        sum_prev = K(:,1:i-1)*a';
        K(:,i) = rate(t + C(i)*h, XA + h*sum_prev);
    end
    XB1 = XA + h*(K*B(1,:)');
    XB2 = XA + h*(K*B(2,:)');
    num_evals = s;
end

function [XB, num_evals, h_next, redo] = rk_step_adaptive(rate,t,XA,h,BT,p,err_des)
    alpha = 1.5;
    redo = false;

    [XB1, XB2, num_evals] = rk_step_embedded(rate,t,XA,h,BT);

    eps_c = norm(XB1 - XB2);
    temp  = (err_des/eps_c)^(1/p);
    h_next = min(0.9*temp, alpha)*h;
    XB = XB1;

    if err_des < eps_c
        redo = true;
    end
end

% variable-step integration wrapper
function [t_list,X_list,h_avg,num_fails,num_evals,h_rec] = rk_variable(rate,tspan,X0,h_ref,BT,p,err_des)

    ti = tspan(1); tf = tspan(2);
    N  = ceil((tf - ti)/h_ref);
    h  = (tf - ti)/N;

    t_list = ti;
    X_list = X0;
    h_rec  = [];
    num_fails = 0;
    num_evals = 0;

    t_now = ti;
    XA = X0;

    while t_now < tf - 1e-14
        redo = true;
        while redo
            h_prev = h;
            h = min(h, tf - t_now);
            [XB, adds, h, redo] = rk_step_adaptive(rate, t_now, XA, h, BT, p, err_des);
            if redo
                num_fails = num_fails + 1;
            end
            num_evals = num_evals + adds;
        end

        h_rec(end+1,1) = h_prev;
        t_now = t_now + h_prev;
        t_list(end+1,1) = t_now;
        X_list(:,end+1) = XB;
        XA = XB;
    end

    h_avg = mean(h_rec);
end

%%%PULSE
%triangle pulse function
function res = triangle_pulse(t,w,h)
t = t*(2/w);
res = 1-min(1*abs(t-1),1);
res = h*res;
end

%triangle pulse function (derivative)
function res = triangle_pulse_derivative(t,w,h)
t = t*(2/w);
res = -sign(t-1).*(abs(t-1)<1);
res = (2*h/w)*res;
end


%%%animation
% AVI recorder: simple box + straight springs animation
function record_animation_avi(filename, frame_rate, t_list, X_list, string_params, rate_factor, w)

    v = VideoWriter(filename, 'Motion JPEG AVI');
    v.FrameRate = frame_rate;
    spf = 1/frame_rate;
    open(v);
    t_end = t_list(end);
    tstart = tic;

    % xs = [string_params.dx*(line-1), string_params.dx*(line)];
    xs = linspace(0,string_params.L,string_params.n+2);

    minX = min(min(X_list(1:string_params.n,:)));
    maxX = max(max(X_list(1:string_params.n,:)));

    Y_max = max(abs(minX),abs(maxX));
    

    fig = figure('Color','white');
    hold on
    axis([min(xs),max(xs),-1.2*Y_max,1.2*Y_max])
    % axis equal; axis([-3 3 -3 3]); hold on;
    

    % line_plot = plot([0,1],0,'k','LineWidth',2);
    dot_plot = plot(0,0, 'ro', 'MarkerFaceColor','r','MarkerSize',5);

    string_plot = plot([0,1],0,'LineWidth',2,'Color','k','MarkerFaceColor','r','markeredgecolor','r','markersize',3);
    wave_loc = string_params.L;
    wavespeed_plot = plot(wave_loc,Y_max, 'g--');
    

    xlabel('x'); ylabel('y');
    title('Strung Masses Motion');
    
    

    frame_index = 1;
    t_frame = toc(tstart);
    while t_frame <= t_end
        while t_list(frame_index) <t_frame
            frame_index = frame_index+1;
        end
        k = frame_index;
        Uk = X_list(:,k);
        Ui = 0;
        Uf_func = string_params.Uf_func;
        Uf = Uf_func(t_frame);
        all_U = [Ui;Uk(1:string_params.n);Uf];
        %short blurb showing how to find x-coord of tracking line
        %x = x-coord of tracking line, t = current time
        %c = wave speed, w = pulse width (in time), L = string length
        c = string_params.wave; L = string_params.L;
        x = L-c*t_frame+.5*w*c;
        x = mod(x,2*L);
        if x > L
            x = 2*L - x;
        end
        set(wavespeed_plot, 'xdata', [x,x], 'ydata', [-Y_max, Y_max])
        
        set(string_plot,'xdata',xs,'ydata',all_U)
        set(dot_plot,'xdata',xs,'ydata',all_U)
        % %draw lines between each
        % for line = 1:length(all_U)-1
        % 
        %     ys = [all_U(line),all_U(line+1)];
        %     % plot(xs, ys, 'k', 'LineWidth', 2)
        % 
        %     set(line_plot,'xdata',xs,'ydata',ys);
        % 
        % end


        % % draw masses as large dots
        % for dot = 1:string_params.n
        %     % plot(, , 'ro', 'MarkerFaceColor','r','MarkerSize',5);
        %     set(dot_plot,'xdata',string_params.dx*dot,'ydata', Uk(dot));
        % end
        % % draw end and beginning as little dots
        % plot(0, Ui, 'ro', 'MarkerFaceColor','r','MarkerSize',2);
        % plot(string_params.L, Uf, 'ro', 'MarkerFaceColor','r','MarkerSize',2)

        title(sprintf('vibrating box, t = %.2f s', t_list(k)));

        drawnow;
        frame = getframe(fig);
        writeVideo(v, frame);
        t_frame = rate_factor*toc(tstart);
    end

    close(v);
    close(fig);

    fprintf('AVI animation saved to %s\n', filename);
end