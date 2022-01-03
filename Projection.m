function out = Projection(theta,gain, thetamax, E)
% Theta - current value of adapting parameter
% gain - The value of the derivative without projection
% thetamax - the limit on the value of theta you'd like to enforce
% E - epsilon, the size of the band in which change in theta is throttled

f = (theta^2-thetamax^2)/(E*thetamax^2);
df = 2*theta/(E*thetamax^2);

if f < 0
    out = gain;
elseif df*gain<=0
    out = gain;
else
    out = gain-gain*f;
end

