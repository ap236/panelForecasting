function p  = densityplotcalcs (y, x, bw, w)
  
  % usage: p  = densityplotcalcs (y, x, bw, w)
  % 
  % calculations for density plots cause the matlab stuff is utter shit
  % 
  % y - data to be evaluated, vector
  % x - grid to be evaluated on, vector
  % bw - bandwidth, scalar (here std of normal)
  % w - weights (optional), vector of same dim as y
      
  % Andreas Pick 
  
  if ~isvector(y)
    error('x must be vector')
  end
  if nargin == 4
    if ~isvector(w)
      error('x must be vector')
    end
    if length(w)~=length(y)
      error('weights are not same length as data')
    end
    w(w==Inf) = 0;
    w = w./sum(w);
  end
  if ~isvector(x)
    error('x must be vector')
  end
  if mean(isnan(y)) == 1
    error('input is vector of nans')
  end
  
  p = nan(size(x));
  if nargin == 3
    for i = 1:length(x)
      p(i) = nanmean(normpdf(y,x(i),bw));
    end
  elseif nargin == 4
    for i = 1:length(x)
      p(i) = nansum(normpdf(y,x(i),bw).*w);
    end    
  else
    error('not right no of nargins')
  end
  
  
end
  
  
  