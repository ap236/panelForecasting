function ppt_cpi_density_plot (dataPath, doHouse)

% Makes density plots for the relative msfe for CPI application 
% for Pesaran, Pick, Timmermann (2025, QE) 

% Andreas Pick - 2025

mctr = 0;
for meanWselInd = 0:2 % 1 - close to mean, 2 - ones stdev away, 0 - all
  mctr = mctr+1;

  if meanWselInd == 1
    WselText = 'mean';
  elseif meanWselInd == 2
    WselText = 'std';
  elseif meanWselInd == 0
    WselText = 'all';
  else
    error('meanWselInd incorrectly set')
  end

  methodList = 1:8;
  nMethods = length(methodList);

  msfei = load([dataPath 'results/cpi/msfei_60_' WselText '.txt']); % N x nMethods
  noObs = load([dataPath 'results/cpi/noObs_60_' WselText '.txt']);

  R2cell = nan(size(msfei));
  lctr = 0;
  for l = 1:nMethods
    lctr = lctr+1;
    R2cell(lctr,:) = msfei(l,:)./msfei(3,:);
  end

  colors = {'#77AC30' '#7E2F8E' '#D95319' '#0072BD' '#A2142F' '#4DBEEE' '#400EEE' '#4D0232'}; % '#477777'};
  lines = {'-','-.',':','-','--','-.',':','-'}; %,'--'};

  if meanWselInd == 1
    titletext = '\kappa_i=0';
  elseif meanWselInd == 2
    titletext = '\kappa_i=\pm 1';
  elseif meanWselInd == 0
    titletext = 'all forecasts';
  end

  figure(mctr+3*doHouse)
  kctr = 0;
  for k = [1 5:8]
    kctr = kctr+1;
    data = R2cell(k,:);
    weights = 1./(noObs(k,:)');
    x = 0.5:0.005:1.5;
    y  = densityplotcalcs (data, x, 0.04);
    xlim([0.5,1.5])
    plot(x,y,'Color',colors{kctr},'LineStyle',lines{kctr},'LineWidth',2)
    hold on
  end
  xlabel('Ratio of MSFE')
  legend('Pooled','Comb (pool)','Comb (FE)','Comb (\omega^*_i)','Emp.Bayes','Location','northwest')
  xticks(0.5:0.1:1.5)
  title(titletext)
  fontsize(gcf,scale=1.5)
  saveas(mctr,[dataPath 'results/cpi/CPI_density_' WselText '.pdf'],'pdf')

  hold off

end