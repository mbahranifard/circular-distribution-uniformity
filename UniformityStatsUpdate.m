%% Uniformity Stats
%% stats
clear all
close all
cd('C:\Users\mbahranifard3\Desktop\Uniformity\')
    numvec = {["81OD.csv","82OS.csv","85OS.csv"],["85OD.csv","82OD.csv","84OS.csv","84OD.csv"],["79OS.csv","80OD.csv","81OS.csv","69OD.csv"]};
    colormat = ["y","r","b"];
    
figure(51)
samphist = tiledlayout(3,1);
selectvec = [3,2,2];


for j = 1:3
    NoMag = figure(j+10);
    t=tiledlayout(numel(numvec{j}),1);

for i=1:numel(numvec{j})
     calccel{i,j}=readmatrix(numvec{j}(i));
    if size(calccel{i,j},2)>21;
        calccel{i,j}(:,22:end)=[];
    end
    %finding the place with most cells and setting that as reference point
    angdist = mean(calccel{i,j},1,'omitnan');
    avvec = [1:numel(angdist)];
    
    angdistnorm = angdist/sum(angdist);
    nominalavg = round(dot(angdistnorm,avvec));
%     [~,highint] = max(angdist);
    newcalc{i,j}= [calccel{i,j}(:,nominalavg:end),calccel{i,j}(:,1:nominalavg-1)];
    colnum = round(numel(angdist)/2);
    newcalc{i,j} = [newcalc{i,j}(:,colnum+1:end),newcalc{i,j}(:,1:colnum)];
%     %plot
%     nexttile((i-1)*3+1)
%     bp = bar(max(newcalc{i,j}(1:4,:)),'barwidth',0.95,'Facecolor',colormat(j));
%     nexttile((i-1)*3+2)
%     bp = bar(max(newcalc{i,j}(5:6,:)),'barwidth',0.95,'Facecolor','b');
    nexttile(i)
    bp = bar(max(newcalc{i,j}(6:end,:)),'barwidth',0.95,'Facecolor',colormat(j));
end
figure(51)
ax(j)=nexttile(j)
bp=bar(max(newcalc{selectvec(j),j}(6:end,:)),'barwidth',0.95,'Facecolor',colormat(j));

end
linkaxes([ax(1) ax(2) ax(3)],'xy')
xticklabels([ax(1),ax(2)],{})
samphist.TileSpacing = 'compact';
samphist.Padding='compact';
xlabel(samphist,["Angle(rad)"],'fontsize',18,'fontWeight','b')
ylabel(samphist,["NFPC"],'fontsize',18,'fontWeight','b')
% set(t,'fontsize',18);
% title(ax(1),'No Magnet','fontsize',10)
% title(ax(2),'Ring','fontsize',10)
% title(ax(3),'Point','fontsize',10)
% radvec = repelem(0,6,3,2,);
% textlab = sprintf("
vec = [.5:21/8:21.5];
set(ax(3),'XLim',[0.3,21.7])
for J = 1:3
ax(J).XAxis.TickLabelInterpreter= 'latex';
ax(J).YAxis.TickLabelInterpreter= 'latex';
ax(J).XAxis.FontSize =22;
% ax(J).XAxis.FontWeight ='bold';
ax(J).YAxis.FontSize =16;
ax(J).YAxis.FontWeight ='bold';
set(ax(J),'YTick',[0:0.2:0.6])
end
% t.TickLabelInterpreter= 'latex';
% ax(3).XAxis.TickLabelFormat      = '\textbf{%g}';
% set(ax(3),'XTick',vec ,'XTickLabel', {'0','\frac{\pi}{4}', '\frac{\pi}{2}$', '\boldmath$\frac{3\pi}{4}$', '\boldmath${\pi}$', '\boldmath$\frac{5\pi}{4}$', '\boldmath$\frac{3\pi}{2}$', '\boldmath$\frac{7\pi}{4}$', '\boldmath$\frac{2\pi}$'})
% set(ax(3),'XTick',vec ,'XTickLabel', {'0','\boldmath$\frac{\pi}{4}$', '\boldmath$\frac{\pi}{2}$', '\boldmath$\frac{3\pi}{4}$', '\boldmath${\pi}$', '\boldmath$\frac{5\pi}{4}$', '\boldmath$\frac{3\pi}{2}$', '\boldmath$\frac{7\pi}{4}$', '\boldmath$\frac{2\pi}$'})
set(ax(3),'XTick',vec ,'XTickLabel', {'0','$\frac{\pi}{4}$', '$\frac{\pi}{2}$', '$\frac{3\pi}{4}$', '${\pi}$', '$\frac{5\pi}{4}$', '$\frac{3\pi}{2}$', '$\frac{7\pi}{4}$', '$\frac{2\pi}$'})

ax(3).XTick
ax(3).XTickLabel

totplot = figure(5);
t=tiledlayout(3,1);
for j=1:3
for i =1:numel(numvec{j})
% numbin = 3; %number of bins to be included
% k = (j-1)*3+i;
outring = NaN([2 size(newcalc{i,j},2)]);
% calculating average of last 2 rows with removign NaN values
test = newcalc{i,j};
[~,twonan] = find(isnan(test(end-1,:)));
[~,onenan] = find(isnan(test(end,:)));
outring = test(end-1:end,:);
outring(:,onenan)=test(end-2:end-1,onenan);
outring(:,twonan)=test(end-3:end-2,twonan);
ringmat(i,:) = max(outring);
end

  Navg(j,:) = mean(ringmat,1);

%
ringthreshmat=ringmat>0.05;
ringthresh = sum(ringthreshmat,2)/size(ringthreshmat,2); 
 Kstd(:,j)=std(ringmat,[],2);
 Kstd(4,1)=NaN;
 KWavg(:,j) = mean(ringmat,2);
 KWavg(4,1)=NaN;
 Kperc(:,j)=ringthresh;
 Kperc(4,1)=NaN;
 

end
%histogram overlays

%% stats
Kcov = Kstd./KWavg;
%% COV plots
[p,table,stats]=anova1(Kcov);
c1=multcompare(stats)
meanbarstd = mean(Kcov,1,'omitnan')
repmat = repelem([1:3],4,1);
for i =1:3
figure(199)
hold all
barploto=bar(i,meanbarstd(i),'facecolor',colormat(i))
barploto.FaceAlpha = 0.3;
% jitterx = -.3 +0.6*randn(numel((Kcov(:,i))),1);
jitterx([2,4]) = 0.2*rand(2,1);
jitterx([1,3])=-1*jitterx([2,4]); 
scatter(repmat(:,i)+jitterx',Kcov(:,i),80,'filled','markerfacecolor',colormat(i),'markeredgecolor','k')
end
errorbar([1,2,3],meanbarstd,std(Kcov,[],1,'omitnan'),'linestyle','none','color','k','CapSize',12)
xtickvec = {'No Magnet','Ring','Point'};
set(gca,'Xtick',[1,2,3],'XtickLabel',xtickvec)
ylabel(['CV'],'fontsize', 20)
set(gca,'XLim',[.5,3.5],'YLim',[0,10])
ax=get(gca);
ax.XAxis.FontSize =22;
ax.XAxis.FontWeight ='bold';
ax.YAxis.FontSize =20;
ax.YAxis.FontWeight ='bold';



%% percent distribution

Kcov = Kstd./KWavg;
%% Kperc plots
Kperc= Kperc*100;
[p,table,stats]=anova1(Kperc);
c=multcompare(stats);
meanbarstd = mean(Kperc,1,'omitnan')
repmat = repelem([1:3],4,1);
for i =1:3
figure(201)
hold all
barploto=bar(i,meanbarstd(i),'facecolor',colormat(i))
barploto.FaceAlpha = 0.3;
% jitterx = -.3 +0.6*randn(numel((Kcov(:,i))),1);
jitterx([2,4]) = 0.2*rand(2,1);
jitterx([1,3])=-1*jitterx([2,4]); 
scatter(repmat(:,i)+jitterx',Kperc(:,i),80,'filled','markerfacecolor',colormat(i),'markeredgecolor','k')
end
errorbar([1,2,3],meanbarstd,std(Kperc,[],1,'omitnan'),'linestyle','none','color','k','CapSize',12)
xtickvec = {'No Magnet','Ring','Point'};
set(gca,'Xtick',[1,2,3],'XtickLabel',xtickvec)
ylabel(['Delivery Adequacy (%)'],'fontsize', 20)
set(gca,'XLim',[.5,3.5],'YLim',[0,100])
ax=get(gca);
ax.XAxis.FontSize =22;
ax.XAxis.FontWeight ='bold';
ax.YAxis.FontSize =20;
ax.YAxis.FontWeight ='bold';
%%




for i=1:3
%     i=3;
nulldist = makedist('uniform','lower',0,'upper',1);
meandist = mean(Navg(i,:));
IOPdiff = Navg(i,:)-meandist; %copied from another code that's why name is IOPdiff
% [h,p]=ttest(IOPdiff);
SEM = std(IOPdiff)/sqrt(length(IOPdiff));               % Standard Errory 
y = normcdf(mean(IOPdiff),0,SEM);
% y = normpdf((IOPdiff))

ts = tinv([ 0.975],length(IOPdiff)-1);      % T-Score
UCI =  ts*SEM;
% p = circ_rtest(Navg(i,:))
% [~,p]=adtest(IOPdiff,'Distribution',nulldist);
end

%f-test
meandist = mean(Navg,2);
SEM = std(Navg,[],2)/sqrt(size(Navg,2));               % Standard Errory 
resdist = Navg-meandist;
[~,p]=vartest2(resdist(:,1),resdist(:,2))


[p2,tbl,stats] = kruskalwallis(KWavg)
linkaxes([ax(:)],'xy')
xticklabels([ax(1),ax(2)],{})
t.TileSpacing = 'compact';
t.Padding='compact';
xlabel(t,["Angle(rad)"],'fontsize',18,'fontWeight','b')
ylabel(t,["Normalized Fluorescent Pixel Count"],'fontsize',18,'fontWeight','b')
% set(t,'fontsize',18);
% title(ax(1),'No Magnet','fontsize',10)
% title(ax(2),'Ring','fontsize',10)
% title(ax(3),'Point','fontsize',10)
% radvec = repelem(0,6,3,2,);
% textlab = sprintf("
vec = [.5:21/8:21.5];
set(ax(3),'XLim',[0.3,21.7])
set(ax(3), 'TickLabelInterpreter', 'latex','XTick',vec ,'XTickLabel', {'0','$\frac{\pi}{4}$', '$\frac{\pi}{2}$', '$\frac{3\pi}{4}$', '${\pi}$', '$\frac{5\pi}{4}$', '$\frac{3\pi}{2}$', '$\frac{7\pi}{4}$', '$\frac{2\pi}$'})
get(ax(3),'XTickLabel')
get(ax(3),'XTick')

for i=1:3
    set(ax(i),'fontsize',15)
end

ax(3).XAxis.FontSize =10;
c = multcompare(stats);

for i =1:3
[test,p] = swtest(KWavg(:,i))
end

[p,table,stats]=anova1(KWavg);
multcompare(stats)

meanbar = mean(KWavg,1,'omitnan')
repmat = repelem([1:3],4,1);

for i =1:3
figure(200)
hold all
barploto=bar(i,meanbar(i),'facecolor',colormat(i))
barploto.FaceAlpha = 0.3;
% jitterx([2,4]) = 0.2*rand(2,1);
% jitterx([1,3])=-1*jitterx([2,4]);
scatter(repmat(:,i)+jitterx',KWavg(:,i),80,'filled','markerfacecolor',colormat(i),'markeredgecolor','k')
end
errorbar([1,2,3],mean(KWavg,1,'omitnan'),std(KWavg,[],1,'omitnan'),'linestyle','none','color','k','CapSize',12)
xtickvec = {'No Magnet','Ring','Point'};
set(gca,'Xtick',[1,2,3],'XtickLabel',xtickvec)
ylabel(['NFPC'],'fontsize', 20)
set(gca,'XLim',[.5,3.5],'YLim',[-0.0001,.27])
ax=get(gca);
ax.XAxis.FontSize =22;
ax.XAxis.FontWeight ='bold';
ax.YAxis.FontSize =20;
ax.YAxis.FontWeight ='bold';
% t.YAxis.FontWeight ='bold';


