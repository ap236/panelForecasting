function [yhat,bFE,bslope] = FEforecastLargeN (y,x,xhat,T)

% usage: [yhat,bFE,bslope] = FEforecastLargeN (y,x,xhat,T)
%
% fixed effects (direct) FEforecast
% where y - TxN matrix of dependent variables
%       x - T*NxK matrix of independent variables w/o intercept
%       xhat - NxK matrix of x in the forecast period
%       T - no. of observations per individual
% output: yhat - Nx1 vector of forecasts
%         bFE - Nx1 vector of estimated fixed effects
%         bslope - Kx1 vector of estimated slope coefficients

% Andreas Pick 

[TN,K] = size(x);
N = TN/T;
xdemean = nan(N*T,K);
ydemean = nan(N*T,1);
for i = 1:N
    xi = x((i-1)*T+1:i*T,:);
    yi = y((i-1)*T+1:i*T);
    nans = isnan(sum([yi xi],2));
    xdemean((i-1)*T+1:i*T,:) = xi - ones(T,1)*mean(xi(nans == 0,:));
    ydemean((i-1)*T+1:i*T) = yi - ones(T,1)*mean(yi(nans == 0));
end
nans = isnan(sum([ydemean xdemean],2));
bslope = xdemean(nans==0,:)\ydemean(nans==0);
bFE = nan(N,1);
for i = 1:N
    bFE(i) = nanmean(y((i-1)*T+1:i*T) - x((i-1)*T+1:i*T,:)*bslope);
end

yhat = bFE + xhat*bslope;

end
