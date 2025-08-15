function yhat = panelRandomCorrCoeffFreqForecasts_wOracle (y, z, zf, N, additionalMaterial)

% usage: yhat = panelRandomCorrCoeffFreqForecasts_wOracle (y, x, x_in_forecast_period, N, additionalMethods)
%
% forecasts from different forecasting methods for balanced panel with random
% coefficients
%
% y - T*Nx1 vector of dependent variable
% x - T*NxK, regressors without intercept
% x_in_forecast_period - NxK, regressors in T+h without intercept
% N - scaler, number of individuals
% additionalMethods - (optional) cell: {1} - include old weights,
%                     {2} - include oracle weights, then contains parms
%                     and disturbances
%
% yhat - Nx2, forecasts from different forecast methods

% Andreas Pick - 2025

[NT,k] = size(z);
T = NT/N;
K = k+1;

nMethods = 2; 
yhat = nan(N,nMethods);

% -----------------------------------------------------------------------
% Forecasts:

% Forecast individual
zz = zeros(K,K);
zy = zeros(K,1);
betaInd = nan(K,N);
Ti = nan(N,1);
yhatInd = nan(N,1);
for ic = 1:N
  zi = [ones(T,1) z((ic-1)*T+1:ic*T,:)];
  yi = y((ic-1)*T+1:ic*T);
  nansi = sum(isnan([zi yi]),2);
  Ti(ic) = sum(nansi==0);
  zi = zi(nansi==0,:);
  yi = yi(nansi==0);
  zz = zz + zi'*zi;
  zy = zy + zi'*yi;
  betaInd(:,ic) = pinv(zi'*zi)*(zi'*yi);

  % Forecast #3: Individual ====
  yhatInd(ic) = betaInd(1,ic) + zf(ic,:)*betaInd(2:end,ic);
end

% Forecast Pooled ====
betaPool = pinv(zz)*zy;
yhatPool = betaPool(1) + (zf*betaPool(2:end))';

% Forecast FE ====
[yhatFE, ~, ~] = FEforecastLargeN (y, z, zf, T);

% Oracle forecasts ===

w = [ones(NT,1) z];
wf = [ones(N,1) zf];

theta = [additionalMaterial{1}, additionalMaterial{2}, additionalMaterial{3}]; %NxK matrix
eta = theta - ones(N,1)*mean(theta,1);
eps = additionalMaterial{4}; % TxN matrix
M = eye(T) - ones(T,T)/T;

QN = 0;
qN = 0;
D_oracle = 0;
psi1 = 0;
psi2 = 0;
hNT_oracle = 0;

% FE
psiFE1 = 0;
psiFE2 = 0;
qFE = 0;
QFE = 0;
hFE = 0;
D_FE = 0;
cFE = 0;
cFE2 = 0;
cNT = 0;

for i = 1:N
  D_oracle = D_oracle + wf(i,:)*eta(i,:)'*eta(i,:)*wf(i,:)'/N;
  QN = QN + w((i-1)*T+1:i*T,:)'*w((i-1)*T+1:i*T,:)/(T*N);
  qN = qN + w((i-1)*T+1:i*T,:)'*w((i-1)*T+1:i*T,:)*eta(i,:)'/(T*N);
  psi1 = psi1 + T/N*eps(:,i)'*w((i-1)*T+1:i*T,:)*inv(w((i-1)*T+1:i*T,:)'*w((i-1)*T+1:i*T,:))*wf(i,:)'*wf(i,:);
  psi2 = psi2 + T/N*eps(:,i)'*w((i-1)*T+1:i*T,:)*inv(w((i-1)*T+1:i*T,:)'*w((i-1)*T+1:i*T,:))*wf(i,:)'*wf(i,:)*eta(i,:)';
  iQiT = inv(w((i-1)*T+1:i*T,:)'*w((i-1)*T+1:i*T,:)/T);
  hNT_oracle = hNT_oracle + wf(i,:)*iQiT*(w((i-1)*T+1:i*T,:)'*eps(:,i)*eps(:,i)'*w((i-1)*T+1:i*T,:)/T)*iQiT*wf(i,:)'/N;

  % FE
  xi = w((i-1)*T+1:i*T,2:end);
  xfdmi = wf(i,2:end) - mean(xi);

  D_FE = D_FE + xfdmi*eta(i,2:end)'*eta(i,2:end)*xfdmi'/N;
  qFE = qFE + xi'*M*xi*eta(i,2:end)'/(T*N);
  QFE = QFE + xi'*M*xi/(T*N);
  bfemb = eta(i,2:end) + eps(:,i)'*M*xi*inv(xi'*M*xi);
  psiFE1 = psiFE1 + T/N*bfemb*(xfdmi'*xfdmi);
  psiFE2 = psiFE2 + T/N*bfemb*(xfdmi'*xfdmi)*eta(i,2:end)';
  iQFE = inv(xi'*M*xi/T);
  hFE = hFE + xfdmi*iQFE*(xi'*M*eps(:,i)*eps(:,i)'*M*xi/T)*iQFE*xfdmi'/N;
  cFE = cFE + (eta(i,2:end)*xfdmi'*mean(eps(:,i)))/N;
  cFE2 = cFE2 + mean(xi)*mean(eps(:,i))/N;
  cNT = cNT + (eps(:,i)'*M*xi*inv(xi'*M*xi)*xfdmi')*mean(eps(:,i))/N;

end
% ind-pool
D_oracle = max([D_oracle - qN'*inv(QN)*qN, 0]);
psi_oracle = psi1*inv(QN)*qN - psi2;

% ind-FE
D_FE = max([D_FE - qFE'*inv(QFE)*qFE, 0]);
psiFE = psiFE1*inv(QFE)*qFE - psiFE2;
cFE = -cFE + qFE'*inv(QFE)*cFE2';

w_oracle(1) = (D_oracle - psi_oracle/T)/(D_oracle + hNT_oracle/T - 2*psi_oracle/T);
if w_oracle(1) < 0
  w_oracle(1) = 0;
elseif w_oracle(1) > 1
  w_oracle(1) = 1;
end
w_oracle(2) = (D_FE - psiFE/T - (cFE - cNT))/(D_FE + hFE/T - 2*psiFE/T);
if w_oracle(2) < 0
  w_oracle(2) = 0;
elseif w_oracle(2) > 1
  w_oracle(2) = 1;
end

yhat(:,1) = w_oracle(1)*yhatInd + (1-w_oracle(1))*yhatPool';
yhat(:,2) = w_oracle(2)*yhatInd + (1-w_oracle(2))*yhatFE;

end
