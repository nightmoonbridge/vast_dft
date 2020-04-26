function [x0, xiter] = newton_root_power(f, fprime, x0, maxIterations, tolerance, constoption)
if nargin < 5
    constoption = 'sd';
    if nargin < 4
        maxIterations = 1e4;
        tolerance = 1e-7;
    end
end
% xiter = x0;
xcumval = zeros(maxIterations,1);
for iteridx = 1 : maxIterations
    y = f(x0);
    if y <= 0
        break;
    end
    
    yprime = fprime(x0);
    
    if(abs(yprime) < eps) % Don't want to divide by too small of a number
        % denominator is too small
        break; %Leave the loop
    end
    
    if strcmpi(constoption, 'sd')
        % sd is the constraint
        x1 = x0 - y/yprime; %Do Newton's computation
    else
        % sb is the constraint
        x1 = x0 - y/yprime; %Do Newton's computation
        mu = 1;
%         x1 = -1; mu = 1;
        while x1<0
            x1 = x0 - mu*y/yprime; %Do Newton's computation
            mu=mu/2;
        end
    end

    
    % If the result is within the desired tolerance
    if (norm(x1 - x0)/norm(x0) <= tolerance) 
        xcumval(iteridx) = x1;
        break; %Done, so leave the loop
    end
    
    x0 = x1; % Update x0 to start the process again
    xcumval(iteridx) = x1;
%     xiter = [xiter;x1];
end
xiter = xcumval(1:iteridx);
end