function [diag, sub, sup, rhs] = Assemble_SOR( u_slice, omega, h, xi, BC, direction )

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
    
    diag =                ones(length(diag_range),1);
     sub = (-omega / 4) * ones(length(sub_range), 1);
     sup = (-omega / 4) * ones(length(sup_range), 1);
     rhs =   (omega / 4) * u_slice(3,diag_range) ...
           + (1 - omega) * u_slice(2,diag_range) ...
           + (omega / 4) * u_slice(1,diag_range) ...
           - (omega / 4) * h^2 * xi;
    
    % Account for boundary conditions.
    if strcmp(direction, 'vertical')
        rhs(1)   = rhs(1)   + (omega / 4) * BC.us;
        rhs(end) = rhs(end) + (omega / 4) * BC.un;
    elseif strcmp(direction, 'horizontal')
        rhs(1)   = rhs(1)   + (omega / 4) * BC.uw;
        rhs(end) = rhs(end) + (omega / 4) * BC.ue;
    else
        error('Invalid direction string.');
    end

end