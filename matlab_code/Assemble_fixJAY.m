function [diag, sub, sup, rhs] = Assemble_fixJAY( uslice, rho, h, xi, BC )

    %%%%%%
    % Assembles the LHS matrix and the RHS vector for the g-system.
    %   diag -- diagonal
    %    sub -- sub-diagonal
    %    sup -- super-diagonal
    %    rhs -- right-hand side vector
    %
    % Ryan Skinner, November 2015
    %%%
    
    N = max(size(uslice));
    
    diag_range = 2:N-1;
     sub_range = 3:N-1;
     sup_range = 2:N-2;
    
    diag = - (2 + rho) * ones(length(diag_range),1);
     sub =               ones(length(sub_range), 1);
     sup =               ones(length(sup_range), 1);
     rhs = -             uslice(diag_range,3) ...
           + (2 - rho) * uslice(diag_range,2) ...
           -             uslice(diag_range,1) ...
           - h^2 * xi;
    
    % Account for boundary conditions.
    rhs(1)   = rhs(1)   - BC.uw;
    rhs(end) = rhs(end) - BC.ue;

end