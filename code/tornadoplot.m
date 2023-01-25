%% ------------------------------------------------------------------------
% This code was jointly written by A/Prof Roslyn Hickson and Dr Andrew Rawlinson
% The code was most recently updated in Januaray 2023
%% ------------------------------------------------------------------------
function tornadoplot(prcctotendvals,legendstring)
  
  % The Tornado Plot
  figure(999),clf;

  % first sort the results
  %non-negatives
  pinds = find(prcctotendvals>=0);  % indexes in prrctotendvals of the non-negs
  posres = prcctotendvals(pinds);   % the non-neg PRCC vals
  pnames = legendstring(pinds);    % the corrensponding variable names (for the tornado plots)
  [pos pnameinds] = sort(posres,'ascend');  % sorted positivies to ascending for tornado plot
  posnames = pnames(pnameinds);             % updated name order
  %negatives
  ninds = find(prcctotendvals<0);
  negres = prcctotendvals(ninds);  
  nnames = legendstring(ninds);
  [neg nnameinds] = sort(negres,'descend');
  negnames = nnames(nnameinds);
  
  tornres = zeros(1,length(pos)+length(neg));
  posind = 1;
  negind = 1;
  posit = pos(posind);
  negit = neg(negind);
  for i=1:(length(pos)+length(neg))
    if posit<abs(negit)
      tornres(i) = posit;
      names{i} = posnames{posind};
      posind = posind+1;
      if posind>length(pos)
        posit = Inf;
      else
        posit = pos(posind);
      end
    else
      tornres(i) = negit;
      names{i} = negnames{negind};
      negind = negind+1;
      if negind >length(neg)
        negit = -Inf;
      else
        negit = neg(negind);
      end
    end
  end
  
  %   figure
  h = barh(tornres);
  hold on
  
  % make sure the x-axis is even across postive and negative
  xlim([-1 1])
  set(gca,'FontSize',15);
  set(gca,'yticklabel',legendstring)
  set(gca,'Ytick',[1:length(legendstring)],'YTickLabel',[1:length(legendstring)])
  
  sep = 1:length(legendstring); % only here you setup
set(gca,'YTickLabel',[])
for n = 1:length(sep)
   text(-1.2,sep(n),names{n},...
      'Interpreter','latex','VerticalAlignment','Middle',...
      'HorizontalAlignment','Left','FontSize',18)%,'FontName','Times New Roman','FontSize',15)
end
  
  xlabel('Partial Rank Correlation Coefficient (PRCC)')

  saveas(gcf,'TornadoPlot.eps')
  saveas(gcf,'TornadoPlot.fig')

return

