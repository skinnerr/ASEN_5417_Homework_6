function [] = Problem_2()

    %%%%%%
    % Solves the inhomogeneous continuity equation for 2D fluid flow within a unit square
    % domain using the SOR method.
    %
    % Ryan Skinner, November 2015
    %%%
    
    Set_Default_Plot_Properties();
    
    %%%
    % Define variables specific to the boundary-value problem.
    %%%
    
    % Solution domain: the closed interval [0,1]x[0,1]. Assume dx = dy = h.
    N = 101;
    x = linspace(0,1,N);
    y = linspace(0,1,N);
    h = x(2) - x(1);
    
    % Source term.
    xi = -1;
    
    % Boundary conditions (u and uprime) defined as cardinal directions (n, s, w, e).
    BC.us = 0;
    BC.uw = 0;
    BC.ue = 0;
    BC.un = 1;
    
    % Relaxation parameters to test.
    omega = [linspace(1.7,1.8,21), linspace(1.8,1.999,61)];
    
    n_to_converge = nan(length(omega));
    
    %%%
    % Solve problem numerically.
    %%%
    
    for i_omega = 1:length(omega)
        om = omega(i_omega);
        fprintf('Working on omega = %6.3f\n', om);
    
        % Initialize the solution, indexed by (x,y), and set BCs.
        u = zeros(N,N);
        u(:,1)   = BC.us;
        u(:,end) = BC.un;
        u(1,:)   = BC.uw;
        u(end,:) = BC.ue;
    
        % Solution norm.
        epsilon = inf;
        conv_crit = 0;

        % Iteration number
        n = 0;

        % Iterate over time steps.
        while epsilon > conv_crit
            n = n + 1;

            u_prev = u;
            % Loop over j (horizontal slices).
            for j = 2:N-1
                [diag, sub, sup, rhs] = Assemble_SOR(u_prev(:,j-1:j+1)', ...
                                                     om, h, xi, BC, 'horizontal');
                if n == 1
                    [LUj.L, LUj.U] = LU_Decompose(diag, sub, sup);
                end
                [sol] = LU_Solve(sub, LUj.L, LUj.U, rhs);
                u(:,j) = [BC.uw; sol; BC.ue];
            end

            u_half = u;
            % Loop over i (vertical slices).
            for i = 2:N-1
                [~, sub, ~, rhs] = Assemble_SOR(u_half(i-1:i+1,:), ...
                                                om, h, xi, BC, 'vertical');
                if n == 1
                    [LUi.L, LUi.U] = LU_Decompose(diag, sub, sup);
                end
                [sol] = LU_Solve(sub, LUi.L, LUi.U, rhs);
                u(i,:) = [BC.us; sol; BC.un];
            end

            epsilon = max(max(abs(u_prev - u)));
            if mod(n,50) == 0
                fprintf('Iteration: %4i, Error Norm: %7.1e\n', n, epsilon);
            end

            if n == 1
                conv_crit = 1e-3 * epsilon;
            end

        end
        
        n_to_converge(i_omega) = n;
        fprintf('Iteration: %4i, Error Norm: %7.1e\n', n, epsilon);
        
    end
    
    %%%
    % Process results.
    %%%
    
    figure();
    plot(omega,n_to_converge);
    xlabel('omega');
    ylabel('Iterations');
    
    figure();
    [C,h] = contour(x,y,u','LineWidth',2);
    clabel(C,h,'FontSize',14,'LabelSpacing',1000);
    axis('equal');
    xlabel('Z');
    ylabel('Y');
    
    figure();
    surf(x,y,u');
    xlabel('Z');
    ylabel('Y');
    
    disp('Done.');
    return
    
end













