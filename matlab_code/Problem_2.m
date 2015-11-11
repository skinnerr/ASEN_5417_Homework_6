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
    omega = [linspace(1.5,1.8,4), linspace(1.805,1.995,80)];
    
    n_to_converge = nan(length(omega));
    
    %%%
    % Solve problem numerically.
    %%%
    
    for i_omega = 1:length(omega)
        om = omega(i_omega);
        fprintf('\nWorking on omega = %6.3f\n', om);
    
        % Initialize the solution, indexed by (x,y), and set BCs.
        u = zeros(N,N);
        u(:,1)   = BC.us;
        u(:,end) = BC.un;
        u(1,:)   = BC.uw;
        u(end,:) = BC.ue;
    
        % Error and convergence measures.
        epsilon = inf;
        conv_crit = 0;

        % Iteration number
        n = 0;

        % Iterate over time steps.
        while epsilon > conv_crit
            
            n = n + 1;
            u_prev = u;
            
            for i = 2:N-1
                for j = 2:N-1
                    u(i,j) = (1 - om) * u_prev(i,j) ...
                             + (om/4) * (        u(i-1, j) ...
                                          +      u(i,   j-1) ...
                                          + u_prev(i+1, j) ...
                                          + u_prev(i,   j+1) ...
                                          - h^2 * xi );
                end
            end
            
            % Calculate error norm.
            epsilon = sum(sum(abs(u_prev - u)));
            
            % Set convergence criterion based on first iteration.
            if n == 1
                conv_crit = 1e-3 * epsilon;
            end
            
            % Abort if diverging.
            if 1e-3*epsilon/conv_crit > 1e5
                error('Solution diverged.');
            end
            
            if (mod(n,10)  == 0) fprintf('.'); end
            if (mod(n,500) == 0) fprintf('\n'); end
        end
        
        % Store and print convergence information.
        n_to_converge(i_omega) = n;
        fprintf('\nIteration: %2i, Error Norm: %7.1e\n', n, 1e-3*epsilon/conv_crit);
        
    end
    
    %%%
    % Process results.
    %%%
    
    figure();
    plot(omega,n_to_converge);
    xlabel('omega');
    ylabel('Iterations');
    xlim([1.5,2]);
    ylim([0,2000]);
    
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













