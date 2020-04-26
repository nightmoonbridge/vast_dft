function xd_opt = calRootFinding(hfcn, x0, options)
% function xd_opt = calRootFinding(hfcn, x0, options)
% A function CALROOTFINDING returns the root of the given function handler.
%
if nargin < 3
    options = optimset('fzero');
    options.Display = 'off';
    options.Constraint = 'sd';

    if nargin < 2
        x0 = [0 1e4]; % initial interval
    end
end

try
    [xd_opt, fval, exitflag, output] = fzero(hfcn, x0, options);
    xd_opt_flag = 1;
catch
    xx = logspace(-10,10,500);
    
    Sddx = 10;
    iter = 1;
    itit = zeros(1e5,1);
    while (Sddx > 1e-5)
        x1 = xx(iter);
        x2 = xx(iter+1);
        Sddx = abs(hfcn(x2) - hfcn(x1));
        itit(iter) = Sddx;
        iter = iter + 1;
        if iter > 1e5
            break;
        end
    end
    itit = itit(1:iter-1);
    
    xd_opt = x1;
    xd_opt_flag = NaN;
end
% ctrfilter.cvxopt_properties.optpara(fbinidx,zidx) = xd_opt;

end