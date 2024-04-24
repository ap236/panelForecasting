function Fout = principalComponentsNaN (y, R, varargin)

% usage: F = principalComponentsNaN (y, R, varargin)
%
% Calculates the factors, F, for the (TxN) matrix y with missing values
% (NaN), the number of factors is R. The EM-algorithm similar to the one
% suggested by Stock and Watson (2002, JBES), however with a standard
% principal components analysis.
%
% if varargin exists, then the eigenvectors are normalised to sum to one
% (instead of being of length one, which is the Matlab default)
%
% Delivers a TxR matrix of factors.

% Andreas Pick 

[T,N] = size(y);
for i = 1:N % normalise in the presence of NaNs
    %y(:,i) = y(:,i) - mean(y(y(:,i)>-999999,i));
    %y(:,i) = y(:,i)./std(y(y(:,i)>-999999,i));
    y(:,i) = y(:,i) - nanmean(y(:,i));
    y(:,i) = y(:,i)./nanstd(y(:,i));
end
if nargin == 3
    nor = 2;
else
    nor = 1;
end

rowNaN = (mean(isnan(y),2) == 1);
Fout = zeros(T,R)*NaN;
y = y(rowNaN == 0,:);
colNaN = (mean(isnan(y),1) == 1);
y = y(:,colNaN == 0);
[T,N] = size(y);

if sum(isnan(y),1)==0 > 3
    % starting values from the variables that don't have missing values if more
    % than three
    ycomplete = y(:,sum(isnan(y),1)==0);
    covmat = cov(ycomplete); % get the error covariance matrix
    [eigvec, eigval] = eig(covmat);
    if nor == 2
        normvec = eigvec(:,end-R+1:end)./(ones(N,1)*sum(eigvec(:,end-R+1:end),1));
    else
        normvec = eigvec(:,end-R+1:end);
    end
    Fnew = ycomplete*normvec;
else
    F = zeros(T,R); % starting values
end
lam = zeros(R,N);
x = y;
x(isnan(y)) = 1; % starting values
cc = 10^(-6)/(T*R); % convergence criterium
iterLimit = 10^4; % max no of iterations
dif = 1;
ctr = 0;
converged = 0;

while converged == 0
    while (dif > cc)

        ctr = ctr+1;
        % E-step:
        ex = F*lam; % expectation of x
        xnew = y;
        xnew(isnan(y)) = ex(isnan(y)); % filling missing value with expectations

        % M-step:
        covmat = cov(x); % get the error covariance matrix
        [eigvec, eigval] = eig(covmat);
        if nor == 2
            normvec = eigvec(:,end-R+1:end)./(ones(N,1)*sum(eigvec(:,end-R+1:end),1));
        else
            normvec = eigvec(:,end-R+1:end);
        end
        Fnew = x*normvec;
        lnew = Fnew\xnew; % OLS regression of xnew on Fnew.

        dif = sum(sum(abs(F - Fnew)));

        F = Fnew;
        lam = lnew;
        x = xnew;
        if ctr == iterLimit
            break
        end
    end

    if ctr < iterLimit
        Fout(rowNaN == 0,:) = F;
        converged = 1;
    elseif ((ctr == iterLimit) && (cc == 10^(-10)/(T*R)))
        ctr = 0;
        cc = cc*100;
    else
        outstr = ['principalComponentsNaN not converged after ' num2str(iterLimit) ' iterations'];
        disp(outstr)
        Fout = NaN;
        converged = 1;
    end
end
