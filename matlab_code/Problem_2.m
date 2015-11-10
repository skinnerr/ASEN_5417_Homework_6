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
    
    % Iteraton parameter as a function of interation number.
    rho = @(iter) 4 * sin(pi * iter / (2 * N))^2;
    
    % Relaxation parameters to test.
    omega = [linspace(1.7,1.8,21), linspace(1.8,1.999,61)];
    omega = 1;
    
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
                [diag, sub, sup, rhs] = Assemble_SOR(u_prev(:,j-1:j+1), ...
                                                     rho(n), h, xi, BC, 'horizontal');
                [LUhoriz.l, LUhoriz.u] = LU_Decompose(diag, sub, sup);
                [sol] = LU_Solve(sub, LUhoriz.l, LUhoriz.u, rhs);
                u(:,j) = [BC.uw; sol; BC.ue];
                u(:,j) = omega * u(:,j) + (1 - omega) * u_prev(:,j);
            end
    
%             surf(x,y,u');
%             xlabel('Z');
%             ylabel('Y');
%             input('asfd');

            u_half = u;
            % Loop over i (vertical slices).
            for i = 2:N-1
                [diag, sub, sup, rhs] = Assemble_SOR(u_half(i-1:i+1,:)', ...
                                                     rho(n), h, xi, BC, 'vertical');
                [LUvert.l, LUvert.u] = LU_Decompose(diag, sub, sup);
                [sol] = LU_Solve(sub, LUvert.l, LUvert.u, rhs);
                u(i,:) = [BC.us; sol; BC.un];
                u(i,:) = omega * u(i,:) + (1 - omega) * u_half(i,:);
            end

            epsilon = sum(sum(abs(u_prev - u)));
            fprintf('Iteration: %2i, Error Norm: %7.1e\n', n, epsilon);

            if n == 1
                conv_crit = 1e-3 * epsilon;
            end
    
%             surf(x,y,u');
%             xlabel('Z');
%             ylabel('Y');
%             input('asfd');

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













