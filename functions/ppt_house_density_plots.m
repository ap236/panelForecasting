function ppt_house_density_plots (dataPath)

% makes density plots for the relative msfe of the house price application
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

  nMethods = 10;
  selectMethod = [1 5 6 9 10];

  msfei = load('-ascii',[dataPath 'results/house/msfei_' WselText]);

  R2cell = nan(size(msfei))';
  for i = 1:nMethods
    R2cell(:,i) = msfei(i,:)./msfei(3,:);
  end

  colors = {'#77AC30' '#7E2F8E' '#D95319' '#0072BD' '#A2142F' '#4DBEEE' '#400EEE'};
  lines = {'-','-.',':','-','--','-.',':'};

  if meanWselInd == 1
    titletext = '\kappa_i=0'; 
  elseif meanWselInd == 2
    titletext = '\kappa_i=\pm 1'; 
  elseif meanWselInd == 0
    titletext = 'all forecasts';
  end

  figure(mctr)
  kctr = 0;
  for k = selectMethod
    kctr = kctr+1;
    data = R2cell(:,k);
    if mctr == 1
      ylim([0,9])
    end
    x = 0.6:.005:1.4;
    xlim([.6,1.4])
    y  = densityplotcalcs (data, x, 0.03); %, weights);
    plot(x,y,'Color',colors{kctr},'LineStyle',lines{kctr},'LineWidth',2)

    hold on
  end
  xlabel('Ratio of MSFE')
  legend('Pooled', 'Comb (pool)', 'Comb (FE)', 'Comb (\omega^*_i)', 'Emp.Bayes')
  title([titletext])
  fontsize(gcf,scale=1.5)

  saveas(mctr,[dataPath 'results/house/house_density_' WselText '.pdf'],'pdf')

  hold off

end