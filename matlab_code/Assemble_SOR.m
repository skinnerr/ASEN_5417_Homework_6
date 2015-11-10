function [diag, sub, sup, rhs] = Assemble_SOR( u_slice, rho, h, xi, BC, direction )

    %%%%%%
    % Assembles the LHS matrix and the RHS vector for the g-system.
    %   diag -- diagonal
    %    sub -- sub-diagonal
    %    sup -- super-diagonal
    %    rhs -- right-hand side vector
    %
    % Ryan Skinner, November 2015
    %%%
    
    N = max(size(u_slice));
    
    diag_range = 2:N-1;
     sub_range = 3:N-1;
     sup_range = 2:N-2;
    
    diag = - (2 + rho) * ones(length(diag_range),1);
     sub =               ones(length(sub_range), 1);
     sup =               ones(length(sup_range), 1);
     rhs = -             u_slice(diag_range,3) ...
           + (2 - rho) * u_slice(diag_range,2) ...
           -             u_slice(diag_range,1) ...
           - h^2 * xi;
    
    % Account for boundary conditions.
    if strcmp(direction,'vertical')
        rhs(1)   = rhs(1)   - BC.us;
        rhs(end) = rhs(end) - BC.un;
    elseif strcmp(direction,'horizontal')
        rhs(1)   = rhs(1)   - BC.uw;
        rhs(end) = rhs(end) - BC.ue;
    else
        error('Invalid direction specified.');
    end

end