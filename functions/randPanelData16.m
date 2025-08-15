function output = randPanelData16 (N, T, h, prms, seedValue, step, uniform_sig2i)

  % usage: output = randPanelData16 (N, T, prms, xTp1, seed, step, optional value)
  %
  % Generates dynamic panel data for Pesaran, Pick, Timmermann 2025
  % in two steps. First the parameters and then the data.
  %
  % Step 1: generate parameters
  % Step 2: read in output of step 1 as parameters and use to generate
  %         data

  % Andreas Pick - 2025

  if step == 1 % Generate parameters

    if nargin < 8
      uniform_sig2i = 0;
    end

    rho0      = prms{1};
    arho      = prms{2};
    g0        = prms{3};
    sig2alpha = prms{4};
    sig2gamma = prms{5};
    Rbx       = prms{6};
    kappa1    = prms{7};
    kappa2    = prms{8};
    sign4Variance = prms{9};

    randn('seed',seedValue);
    rand('seed',seedValue+1);

    if sig2gamma > 0
      g0i = (1:N)';
      g0i(g0i<=N/2) = 0.2/3;
      g0i(g0i>N/2) = 0.4/3;
    else
      g0i = g0;
    end

    if sig2alpha > 0
      alpha0i = (1:N)';
      alpha0i(alpha0i<=N/2) = 2/3;
      alpha0i(alpha0i>N/2) = 4/3;
    else
      alpha0i = 1;
    end

    if g0 > 0 % panel ARX
      ARX = 1;
      phi = Rbx*sqrt(sig2alpha);
      ppi = Rbx*sqrt(sig2gamma);
      sig2zeta = sig2gamma - ppi^2;
      sig2eta = sig2alpha - phi^2;

      mux = (randn(N,1).^2 - 1)/sqrt(2);
      alphai = alpha0i + phi*mux + sqrt(sig2eta)*randn(N,1);
      gammas = g0i + ppi*mux + sqrt(sig2zeta)*randn(N,1);
    else % panel AR
      ARX = 0;
      alphai = alpha0i + randn(N,1)*sqrt(sig2alpha);
      mux = 0;
      gammas = zeros(N,1);
    end

    rhox = rand(N,1)*0.95; 
    rhoy = rho0 + rand(N,1)*arho - arho/2;

    rhoy(rhoy>0.986) = 0.986;
    rhoy(rhoy<-0.986) = -0.986;

    % Generate error terms
    if uniform_sig2i == 0
      sig2e = .5 + .5*randn(N,1).^2;
    else
      sig2e = 1/(0.2 + rand(N,1)*0.6);
    end
    sig2x = .5 + .5*randn(N,1).^2;

    if ARX == 0
        wiTp1 = alphai./(1-rhoy) ...  % mean
          + sign4Variance.*sqrt(sig2e./(1-rhoy.^2));
    elseif ARX == 1
      wiTp1 = nan(N,2); % y_iT and x_iT+1

      varyterm = gammas.^2.*sig2x./(1-rhoy.^2).*(1 + 2.*rhoy.*rhox./(1-rhoy.*rhox));

      wiTp1(:,1) = (alphai + mux.*gammas)./(1-rhoy) ...
          + sign4Variance.*sqrt(sig2e./(1-rhoy.^2) + varyterm);

      wiTp1(:,2) = mux + sign4Variance.*sig2x;
    end

    output{1} = alphai;
    output{2} = mux;
    output{3} = rhoy;
    output{4} = rhox;
    output{5} = gammas;
    output{6} = sig2x;
    output{7} = sig2e;
    output{8} = kappa1;
    output{9} = kappa2;
    output{10} = wiTp1; % predictors in the forecast period

  elseif step == 2 % given prms, generate data

    alphai = prms{1};
    mux    = prms{2};
    rhoy   = prms{3};
    rhox   = prms{4};
    gammas = prms{5};
    sig2x  = prms{6};
    sig2e  = prms{7};
    kappa1 = prms{8};
    kappa2 = prms{9};
    wiTp1  = prms{10};

    randn('seed',seedValue);
    % Generate data

    T = T+h;
    if kappa1 > 0
      Tinit = 1;
    else
      Tinit = 100;
    end

    nu = randn(T+Tinit,N);
    epsilon = (randn(T+Tinit,N).^2-1)/sqrt(2);

    y = nan(T+Tinit,N);
    xi = nan(T+Tinit,N);
    x = nan(T+Tinit,N);

    % initialisation:
    if sum(abs(gammas)) > 0
      if kappa1 < 0
        y(1,:) = ((alphai + gammas.*mux)./(1-rhoy))';
        xi(1,:) = zeros(1,N);
      else
        sigma2i0 = (gammas.^2.*sig2x + sig2e)./(1-rhoy.^2);
        y(1,:) = kappa1.*((alphai + gammas.*mux)./(1-rhoy))' ...
          + kappa2.*sqrt(sigma2i0)'.*randn(1,N);
        xi(1,:) = zeros(1,N);
      end

      for i = 2:(T+Tinit)
        xi(i,:)  = rhox'.*xi(i-1,:) + sqrt(sig2x)'.*sqrt(1-rhox.^2)'.*nu(i,:);
        x(i,:)   = mux' + xi(i,:);
        y(i,:) = alphai' + rhoy'.*y(i-1,:) + gammas'.*x(i,:) + sqrt(sig2e)'.*epsilon(i,:);
      end

      y(end,:) = alphai' + rhoy'.*wiTp1(:,1)' + gammas'.*wiTp1(:,2)' + sqrt(sig2e)'.*epsilon(end,:); % make forecast period y consistent with w_i,T+1

      output{1} = y(Tinit+h:end,:);
      output{2} = y(Tinit:end-h,:);
      output{3} = x(Tinit+h:end,:);
      output{4} = (ones(T-h,1)*sqrt(sig2e)').*epsilon(Tinit+1:end-h,:);

    else % panel AR

      if kappa1 < 0
        y(1,:) = (alphai./(1-rhoy))';
      else
        sigma2i0 = sig2e./(1-rhoy.^2);
        y(1,:) = kappa1.*(alphai./(1-rhoy))' ...
          + kappa2.*sqrt(sigma2i0)'.*randn(1,N);
      end

      for i = 2:(T+Tinit)
        y(i,:) = alphai' + rhoy'.*y(i-1,:) + sqrt(sig2e)'.*epsilon(i,:);
      end

      y(end,:) = alphai' + rhoy'.*wiTp1(:,1)' + sqrt(sig2e)'.*epsilon(end,:); % make forecast period y consistent with w_i,T+1

      output{1} = y(Tinit+1:end,:);
      output{2} = y(Tinit:end-1,:);
      output{3} = [];
      output{4} = (ones(T,1)*sqrt(sig2e)').*epsilon(Tinit+1:end,:);
    end

  else
    error('input "step" needs to be either 1 or 2')
  end

end
