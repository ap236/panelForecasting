% making density plots for the relative msfe
% data: MSFE for each individual

% Andreas Pick 

close('all')

savePlot = 0;
do_cdf = 1;

nModels = 2;
nMethods = 8;

msfei = cell(nModels,1);
nObs = cell(nModels,1);
for j = 1:nModels
  msfei{j,1} = load('-ascii',[savePath 'housePrice_msfei_model_' num2str(j) '.txt']);
end

R2cell = cell(nModels,1);
for j = 1:nModels
  for i = 1:nMethods
    R2cell{j,1}(:,i) = msfei{j,1}(i,:)./msfei{j,1}(3,:);
  end
end

colors = {'#77AC30' '#7E2F8E' '#D95319' '#0072BD' '#A2142F' '#4DBEEE' '#400EEE'};
lines = {'-','-.',':','-','--','-.',':'};

titletext{1,1} = 'SAR';
titletext{2,1} = 'SARX';

cntr = 0;
for j = 1:2 % model
  cntr = cntr+1;
  figure(cntr)
  kctr = 0;
  for k = [1 2 4:nMethods]
    kctr = kctr+1;
    data = R2cell{j,1}(:,k);
    x = 0.6:.005:1.4;
    xlim([.6,1.4])
    y  = densityplotcalcs (data, x, 0.03);
    if do_cdf == 1
      y = cumsum(y)./sum(y);
      ylim([0,1])
    end
    plot(x,y,'Color',colors{kctr},'LineStyle',lines{kctr},'LineWidth',2)
    hold on
  end
  xlabel('Ratio of MSFE')
  if do_cdf == 0
    legend('Pool', 'RE', 'FE', 'Emp.Bayes', 'Hier.Bayes', 'Comb (pool)', 'Comb (FE)')
  else
    legend('Pool', 'RE', 'FE', 'Emp.Bayes', 'Hier.Bayes', 'Comb (pool)', 'Comb (FE)','Location','northwest')
  end
  title([titletext{j,1}])

  if savePlot == 1
    if do_cdf == 0
      saveas(cntr,['housePrice_density_' titletext{j,1} '_pdf.pdf'],'pdf')
    else
      saveas(cntr,['housePrice_density_' titletext{j,1} '_cdf.pdf'],'pdf')
    end
  end
end
hold off
