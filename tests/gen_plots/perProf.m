function [getfig]=perProf(D,xlabelstr,ylabelstr,titlestr,legendCellStr, skipSolver, str)
% Input : np is no. of problems or no. of rows of D
%         ns is no. of solvers or no. of columns of D
%         D is the data matrix, Assume D has all the positive entries
%         fbestData array of order np x ns saving the fbest values for all the algorithms
%         targetfbest array of 1 x np will be used to penalize the entries in D matrix, 

tol=1e-8; % if some entry is less than this tolerance modify all the entries of matrix D
[np,ns]=size(D); 

r=D; % initialize
% (1) if any value in a row is < tol, then add a dummy variable for all the values in that row.
indices=(r<tol);
for j=1:np % 1:nrows
    if sum(indices(j,:))>0   % if sum value is small in a row
       r(j,:)=r(j,:)+1;      % add 1 to all the values of the row
    end
end

% (2) Divide each row by the min. value of that row
r=r./min(r,[],2);            % divide each row with the min. of each row
rmax=ceil(max(max(r)));      % find the upper bound for the data


% Protection:  remove the outliers, like NaN, maybe 0 as well
r(isnan(r))=rmax; % assign NaN values to rmax


r=sort(r,1); % sort data in each row from low to high

% The following colcode is consistent with boxPlotsAllDimRegimes.m to have the same color of the algo. in both boxplots and perprof plots 
colcode =[1         0         0        % red
          0.9294    0.6941    0.1255   % yellow
          0         0.4471    0.7412   % blue
          0.3922    0.8314    0.0745   % green
          0.7176    0.2745    1.0000   % purple             
          0         1         1        % cyan
          1.0000    0.4118    0.1608   % orange
          1.0000    0.0745    0.6510   % pink
          0.5020    0.5020    0.5020   % gray
          0.30      0.75      0.93     % sky blue
          0.47      0.67      0.19     % olive
          0.64      0.08      0.18     % dark red
          0         0         0.5      % navy blue
          ];

% The following block, uses different marker styles along with line styles
figure
linestyle=cell(1,12);markerstyle=cell(1,12);
idxcell=1;
linestyle{idxcell}='--';markerstyle{idxcell}='*'; idxcell=idxcell+1;
linestyle{idxcell}='-';markerstyle{idxcell}='x';  idxcell=idxcell+1;
linestyle{idxcell}=':';markerstyle{idxcell}='v';  idxcell=idxcell+1;
linestyle{idxcell}='-.';markerstyle{idxcell}='^'; idxcell=idxcell+1;
linestyle{idxcell}='--';markerstyle{idxcell}='*'; idxcell=idxcell+1;
linestyle{idxcell}='-';markerstyle{idxcell}='x';  idxcell=idxcell+1;
linestyle{idxcell}=':';markerstyle{idxcell}='v';  idxcell=idxcell+1;
linestyle{idxcell}='-.';markerstyle{idxcell}='^'; idxcell=idxcell+1;
linestyle{idxcell}='--';markerstyle{idxcell}='*'; idxcell=idxcell+1;
linestyle{idxcell}='-';markerstyle{idxcell}='x';  idxcell=idxcell+1;
linestyle{idxcell}=':';markerstyle{idxcell}='v';  idxcell=idxcell+1;
linestyle{idxcell}='-.';markerstyle{idxcell}='^'; idxcell=idxcell+1;

if ~isempty(skipSolver)
   colcode(skipSolver,:)=[];linestyle(skipSolver)=[];markerstyle(skipSolver)=[]; 
   idxcell=idxcell-1;
end

for i=1:ns
   [xpts,ypts]=stairs( r(:,i),(cumsum(r(:,i)<rmax))/np);
   xs=[xpts;rmax];ys=[ypts;ypts(end)]; % adding an extra point so that the graphs start from (1,0) point
   if ns>9
      getfig=plot(xs,ys,linestyle{i},'LineWidth',2); 
   else
      getfig=plot(xs,ys,linestyle{i},'LineWidth',2,'Color',colcode(i,:));
   end
   hold on
end
hold off

if rmax==1,rmax=1.1;
end
newxl=[1 rmax];
axis([newxl , 0  1]);
set(gca,'Units','normalized','YTick',0:0.25:1,... 
   'FontUnits','points','FontSize',13,'FontName','Times','FontWeight','normal','LineWidth',1.5);

% Latex typesetting by using interpreter as LATEX    
xlabel({sprintf('$%s$',xlabelstr)},'FontUnits','points','FontSize',15,'FontName','Times','FontWeight','normal','interpreter','latex');
ylabel({sprintf('$%s$',ylabelstr)},'FontUnits','points','FontSize',14,'FontName','Times','FontWeight','normal','interpreter','latex');

if isempty(legendCellStr)
   legend(cellfun(@(u) append('Set', num2str(u) ), num2cell(1:ns) ,'UniformOutput',false),...
    'FontSize',12,'FontName','Times','interpreter','latex','Location','southeast');
else
   legend( legendCellStr,...
    'FontSize',12,'FontName','Times','interpreter','latex','Location','southeast'); 
end

title(sprintf('%s',titlestr),'FontUnits','points','FontWeight','normal','FontSize',15,'FontName','Times',...
    'Interpreter','latex','Color','k');

saveas(gcf,sprintf('perprofEx123%s',str),'fig');
saveas(gcf,sprintf('perprofEx123%s',str),'epsc');
% figure('visible','on');
% saveas(gcf,'PerProf');

end


