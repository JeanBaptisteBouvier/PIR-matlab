function orbit = halo_orbit_interpolation_2(cr3bp, orbit, abacus, params, cst, varargin)

%--------------------------------------------------------------------------
% Preliminary checks
%--------------------------------------------------------------------------
if(nargin ~= 7)
    error('Wrong number of arguments');
end

if(~strcmp(cr3bp.name,'EARTH+MOON'))
    error('Wrong CRTBP. Only EARTH+MOON system is authorized.');
end

%--------------------------------------------------------------------------
% First guess from abacus
%--------------------------------------------------------------------------
switch(varargin{1})
    
    case 'Az'
        Az = varargin{2}/cr3bp.L;
        if(Az > abacus.AzLimit(1) && Az < abacus.AzLimit(2))
            %Interpolate
            fit = interpol(abacus, abacus.Az, Az);
        else
            error('The desired vertical extension is out of bounds. Current limits: [%5.5f, %5.5f] km',...
                abacus.AzLimit(1)*cr3bp.L, abacus.AzLimit(2)*cr3bp.L);
        end
    
    case 'X0'
        X0 = varargin{2};
        if(X0 > abacus.X0Limit(1) && X0 < abacus.X0Limit(2))
            %Interpolate
            fit = interpol(abacus, abacus.X0, X0);
        else
            error('The desired X0 is out of bounds. Current limits: [%5.5f, %5.5f] km',...
                abacus.X0Limit(1), abacus.X0Limit(2));
        end
        
    case 'C'
        C = varargin{2};
        if(C > abacus.CLimit(1) && C < abacus.CLimit(2))
            %Interpolate
            fit = interpol(abacus, abacus.C, C);
        else
            error('The desired Jacobi constant C is out of bounds. Current limits: [%5.5f, %5.5f] km',...
                abacus.CLimit(1), abacus.CLimit(2));
        end
        
end

% First guess
yv0 = fit.f(1:6);

% Specific case of the SOUTHERN family
if(strcmp(orbit.family,cst.orbit.family.SOUTHERN))
    yv0(3) = -yv0(3);
end


%--------------------------------------------------------------------------
% Refinement
%--------------------------------------------------------------------------
orbit = orbit_refinement(cr3bp, orbit, params, yv0, cst);

end

%--------------------------------------------------------------------------
% Interpolate the initial condition y0 to get data(y0) = value.
%--------------------------------------------------------------------------
function fit = interpol(abacus, data, value)

fit.half = 2;
[~, array_position] = min(abs(data - value));
mini = max(array_position - fit.half, 1);
maxi = min(array_position + fit.half, length(data));
fit.degree = length(mini:maxi)-1;
fit.x =  data(mini:maxi)';
for count =1:6
    fit.y(:,count) =  abacus.initialConditions(mini:maxi,count);
    %Fitting for every dimension of the state (6)
    [fit.p, ~, fit.mu] = polyfit(fit.x,fit.y(:,count),fit.degree);
    %Evaluation
    fit.f(count) = polyval(fit.p,value,[],fit.mu);
end

end