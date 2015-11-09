function [] = Problem_1()

    %%%%%%
    % Solves the inhomogeneous continuity equation for 2D fluid flow within a unit square
    % domain using the ADI method.
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
    BC.upn = 0; % This whole code assumes a value of zero.
    
    % Initialize the solution, indexed by (x,y), and set BCs.
    u = zeros(N,N);
    
    % Fixed iteraton parameter.
    rho = 0.002;
    
    %%%
    % Solve problem numerically.
    %%%
    
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
            [diag, sub, sup, rhs] = Assemble_fixJAY(u_prev(:,j-1:j+1), rho, h, xi, BC);
            if n == 1
                [LUj.l, LUj.u] = LU_Decompose(diag, sub, sup);
            end
            [sol] = LU_Solve(sub, LUj.l, LUj.u, rhs);
            u(:,j) = [BC.uw; sol; BC.ue];
        end
        
        u_half = u;
        % Loop over i (vertical slices).
        for i = 2:N-1
            [diag, sub, sup, rhs] = Assemble_fixEYE(u_half(i-1:i+1,:), rho, h, xi, BC);
            if n == 1
                [LUi.l, LUi.u] = LU_Decompose(diag, sub, sup);
            end
            [sol] = LU_Solve(sub, LUi.l, LUi.u, rhs);
            u(i,:) = [BC.us; sol; sol(end)];
        end
        
        epsilon = sum(sum(abs(u_prev - u)));
        fprintf('Iteration: %2i, Error Norm: %7.1e\n', n, epsilon);
        
        if n == 1
            conv_crit = 1e-3 * epsilon;
        end
        
    end
    
    %%%
    % Process results.
    %%%
        
    [C,h] = contour(x,y,u','LineWidth',2);
    clabel(C,h,'FontSize',14,'LabelSpacing',1000);
    axis('equal');
    xlabel('Z');
    ylabel('Y');
    
    disp('Done.');
    return
    
end













