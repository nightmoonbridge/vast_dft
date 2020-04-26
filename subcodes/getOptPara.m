function ctrfilter = getOptPara(uvrb,d_part,hz,ctrfilter,idxset)
if nargin < 5
    idxset = [1 1];
end

% idxset contains an index set of [frequency bin, zone]
fbinidx = idxset(1);
zidx = idxset(2);

x0 = [1e-6 1e4]; % initial interval
options = optimset('fzero');
options.Display = 'off';
options.Constraint = ctrfilter.cvxopt_properties.const;
% Using the Newton's method or a built-in function in MATLAB
if ctrfilter.cvxopt_properties.findopt
    if strcmpi(ctrfilter.cvxopt_properties.opttype,'min_sb')
        % find the optimal para fulfills sd <= const
        ed = ctrfilter.cvxopt_properties.tarval;
        if strcmpi(ctrfilter.cvxopt_properties.const, 'sd')
            Sd = @(x) sum(abs(uvrb).^2 ./ abs(d_part + x).^2) - ed;
        else % which is for the normalized sd (nsd)
            Sd = @(x) sum(abs(uvrb).^2 ./ abs(d_part + x).^2) / hz- ed;
%             Sd = @(x) hz / sum(abs(uvrb).^2 ./ abs(d_part + x).^2) - ed;
        end
        
%         options.Constraint = 'sd';
        
        % the Newton's method
%         Sdprime = @(x) -2*sum(abs(uvrb).^2 ./ abs(d_part + x).^3);
%         mu_old = ctrfilter.cvxopt_properties.initpara;
%         ctrfilter.cvxopt_properties.optpara = ...
%             newton_root_power(Sd, Sdprime, mu_old);

        xd_opt = calRootFinding(Sd, x0, options);
        ctrfilter.cvxopt_properties.optpara(fbinidx,zidx) = xd_opt;
    else
        vd = hz'*hz;
        % find the optimal para fulfills sb <= const
        eb = ctrfilter.cvxopt_properties.tarval;
%         Sb = @(x) hz'*hz - sum((2*x + d_part).*abs(uvrb).^2 ./ abs(d_part + x).^2) - eb;
        if strcmpi(ctrfilter.cvxopt_properties.const, 'sb')
          Sb = @(x) (vd - sum((abs(uvrb).^2.*(2*x + d_part)) ./ abs(x + d_part).^2)) - eb;
            
            % this also will give the same solution
%             Sb = @(x) hz'*hz - sum((abs(uvrb).^2.*(2*x + x^2*d_part)) ./ abs(1 + x*d_part).^2) - eb;
        else % which is for the normalized sb (nsb)
            Sb = @(x) (vd - sum((abs(uvrb).^2.*(2*x + d_part)) ./ abs(x + d_part).^2))/vd - eb;
        end
        xb_opt = calRootFinding(Sb, x0, options);
        ctrfilter.cvxopt_properties.optpara(fbinidx,zidx) = xb_opt;
    end
else
    ctrfilter.cvxopt_properties.optpara = 1;
end

end