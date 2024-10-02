% Plot mesh data produced by PolarMaxPeakesMesh.py
close all; clear all;
load('matlabData/TorusCurves.mat')

labelSize =11;
ticklengthUn = 0.2;

Epsa = load('pythonData/MaxPeaksMesh_epsa.txt');
Epsp = load('pythonData/MaxPeaksMesh_epsp.txt');
PeakVals = load('pythonData/MaxPeaksMesh.txt');
FreqVals = load('pythonData/MaxPeaksMesh_freq.txt');
%sky = load('sky.txt');
cmap0 = jet(512);
cmap0(2:3,:)=repmat([1 1 1],2,1);
cmap0(1,:)=repmat([0.7 0.7 0.7],1,1);

FreqVals(PeakVals<0.01)=NaN;
FreqVals(PeakVals>0.7)=NaN;
PeakVals(PeakVals<0.01)=NaN;
PeakVals(PeakVals>0.7)=NaN;

num = size(Epsa,1)/2;
Epsa1 = Epsa(1:num,:);
Epsa2 = Epsa(num+1:end,:);
Epsp1 = Epsp(1:num,:);
Epsp2 = Epsp(num+1:end,:);
PeakVals1 = PeakVals(1:num,:);
PeakVals2 = PeakVals(num+1:end,:);
FreqVals1 = FreqVals(1:num,:);
FreqVals2 = FreqVals(num+1:end,:);

%% Variability of total intensity peak high in %
figure(1); clf; hold on;
%plot(Epsa1(:),Epsp1(:),'.k','MarkerSize',3);
scatter(Epsa2(:),Epsp2(:),[],-0.001*PeakVals2(:)./PeakVals2(:),'filled','s')
scatter(Epsa1(:),Epsp1(:),[],PeakVals1(:),'filled','s');
plot([0,0,0,-0.01],[0.295,0.355,0.415,0.355],'o','MarkerFaceColor', [0. 0. 0.],'MarkerEdgeColor', [0. 0. 0.],'linewidth',1,'MarkerSize',7)
for i=4:length(Br_tr)
    if ~isempty(Br_tr{i})
        [perv, idx] = sort(arrayfun(@(x)x.parameter(ind.epsp),Br_tr{i}.point));
        parv = arrayfun(@(x)x.parameter(ind.epsa),Br_tr{i}.point);
        parv = parv(idx);
        plot(parv,perv,'-','Color', [0. 0. 0.],'linewidth',2);
    end
end
% for i=1:length(Br_pd)
%     if ~isempty(Br_pd{i})
%         [perv, idx] = sort(arrayfun(@(x)x.parameter(ind.epsp),Br_tr{i}.point));
%         parv = arrayfun(@(x)x.parameter(ind.epsa),Br_tr{i}.point);
%         parv = parv(idx);
%         plot(parv,perv,'--','Color', [0. 0. 0.],'linewidth',2);
%     end
% end
hold off;

bd = floor(max(PeakVals1(:))*512);
cmap = cmap0(1:bd,:);
colormap(cmap);
a = colorbar();
a.Label.String = '(max-min)/max in %';
title('Variability of peak total intensity - upsweap $\varepsilon_a$','interpreter','latex');
xlabel('$\varepsilon_a$','interpreter','latex');
ylabel('$\varepsilon_p$','interpreter','latex');

box on;
xAuxTicks = [-0.03, -0.01, 0.01];
yAuxTicks = [0.2, 0.4];
xlim([-0.035, 0.015]);
ylim([0.1, 0.5]);
hold off;
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,8,5])
set(gca,'position',[0.005,0.01,0.99,0.98],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',1.5) %[0.07,0.10,0.92,0.88]
xticklabels({'',''})
yticklabels({'',''})
hgexport(gcf, 'img/Fig2a.eps', hgexport('factorystyle'), 'Format', 'eps');

figure(2); clf; hold on;
%plot(Epsa1(:),Epsp1(:),'.k','MarkerSize',3);
scatter(Epsa1(:),Epsp1(:),[],-0.001*PeakVals1(:)./PeakVals1(:),'filled','s')
scatter(Epsa2(:),Epsp2(:),[],PeakVals2(:),'filled','s');
plot([0,0,0,-0.01],[0.295,0.355,0.415,0.355],'o','MarkerFaceColor', [0. 0. 0.],'MarkerEdgeColor', [0. 0. 0.],'linewidth',1,'MarkerSize',7)
for i=1:4
    if ~isempty(Br_tr{i})
        [perv, idx] = sort(arrayfun(@(x)x.parameter(ind.epsp),Br_tr{i}.point));
        parv = arrayfun(@(x)x.parameter(ind.epsa),Br_tr{i}.point);
        parv = parv(idx);
        plot(parv,perv,'-','Color', [0. 0. 0.],'linewidth',2);
    end
end
% for i=1:length(Br_pd)
%     if ~isempty(Br_pd{i})
%         [perv, idx] = sort(arrayfun(@(x)x.parameter(ind.epsp),Br_tr{i}.point));
%         parv = arrayfun(@(x)x.parameter(ind.epsa),Br_tr{i}.point);
%         parv = parv(idx);
%         plot(parv,perv,'--','Color', [0. 0. 0.],'linewidth',2);
%     end
% end
hold off;

bd = floor(max(PeakVals1(:))*512);
cmap = cmap0(1:bd,:);
colormap(cmap);
a = colorbar();
a.Label.String = '(max-min)/max in %';
title('Variability of peak total intensity - downsweap $\varepsilon_a$','interpreter','latex');
xlabel('$\varepsilon_a$','interpreter','latex');
ylabel('$\varepsilon_p$','interpreter','latex');

box on;
xAuxTicks = [-0.03, -0.01, 0.01];
yAuxTicks = [0.2, 0.4];
xlim([-0.035, 0.015]);
ylim([0.1, 0.5]);
hold off;
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,8,5])
set(gca,'position',[0.005,0.01,0.99,0.98],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',1.5) %[0.07,0.10,0.92,0.88]
xticklabels({'',''})
yticklabels({'',''})
hgexport(gcf, 'img/Fig2c.eps', hgexport('factorystyle'), 'Format', 'eps');

%% Frequency information
% figure(3); clf; hold on;
% scatter(Epsa1(:),Epsp1(:),[],FreqVals1(:),'filled');
% scatter(Epsa2(:),Epsp2(:),[],FreqVals2(:),'Marker','x');
% hold off;
% colormap(winter);
% colorbar();

%% Plot mesh data produced by PolarMaxPeakesMeshX.py
load('TorusCurves.mat')

Epsa = load('pythonData/MaxPeaksMeshX_epsa.txt');
Epsp = load('pythonData/MaxPeaksMeshX_epsp.txt');
PeakVals = load('pythonData/MaxPeaksMeshX.txt');

sky = load('sky.txt');
cmap = jet(256);
cmap(2,:)=repmat([1 1 1],1,1);
cmap(1,:)=repmat([0.7 0.7 0.7],1,1);

num = size(Epsa,1)/2;

xinds1 = find(isnan(PeakVals(1:num,:)));
xinds2 = find(isnan(PeakVals(num+1:end,:)));
PeakVals(PeakVals<0.01)=NaN;
PeakVals(PeakVals>0.99)=NaN;

Epsa1 = Epsa(1:num,:);
Epsa2 = Epsa(num+1:end,:);
Epsp1 = Epsp(1:num,:);
Epsp2 = Epsp(num+1:end,:);
PeakVals1 = PeakVals(1:num,:);
PeakVals2 = PeakVals(num+1:end,:);
Epsa1X = Epsa1(xinds1);
Epsp1X = Epsp1(xinds1);
Epsa2X = Epsa2(xinds2);
Epsp2X = Epsp2(xinds2);
 
%% Variability of x intensity peak high in %
figure(3); clf; hold on;
%plot(Epsa1(:),Epsp1(:),'.k','MarkerSize',3);
plot(Epsa1X(:),Epsp1X(:),'xk','Color',[0.5, 0.5, 0.5],'MarkerSize',6);
scatter(Epsa2(:),Epsp2(:),[],-0.001*PeakVals2(:)./PeakVals2(:),'filled','s')
scatter(Epsa1(:),Epsp1(:),[],PeakVals1(:),'filled','s');
plot([0,0,0,-0.01],[0.295,0.355,0.415,0.355],'o','MarkerFaceColor', [0. 0. 0.],'MarkerEdgeColor', [0. 0. 0.],'linewidth',1,'MarkerSize',7)
for i=5:length(Br_tr)
    if ~isempty(Br_tr{i})   
        [perv, idx] = sort(arrayfun(@(x)x.parameter(ind.epsp),Br_tr{i}.point));
        parv = arrayfun(@(x)x.parameter(ind.epsa),Br_tr{i}.point);
        parv = parv(idx);
        plot(parv,perv,'-','Color', [0. 0. 0.],'linewidth',2);
    end
end
% for i=1:length(Br_pd)
%     if ~isempty(Br_pd{i})
%         [perv, idx] = sort(arrayfun(@(x)x.parameter(ind.epsp),Br_tr{i}.point));
%         parv = arrayfun(@(x)x.parameter(ind.epsa),Br_tr{i}.point);
%         parv = parv(idx);
%         plot(parv,perv,'--','Color', [0. 0. 0.],'linewidth',2);
%     end
% end
hold off;

bd = floor(max(PeakVals1(:))*512);
cmap = cmap0(1:bd,:);
colormap(cmap);
a = colorbar();
a.Label.String = '(max-min)/max in %';
title('Variability of peak total intensity - downsweap $\varepsilon_a$','interpreter','latex');
xlabel('$\varepsilon_a$','interpreter','latex');
ylabel('$\varepsilon_p$','interpreter','latex');

box on;
xAuxTicks = [-0.03, -0.01, 0.01];
yAuxTicks = [0.2, 0.4];
xlim([-0.035, 0.015]);
ylim([0.1, 0.5]);
hold off;
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,8,5])
set(gca,'position',[0.005,0.01,0.99,0.98],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',1.5) %[0.07,0.10,0.92,0.88]
xticklabels({'',''})
yticklabels({'',''})
hgexport(gcf, 'img/Fig2b.eps', hgexport('factorystyle'), 'Format', 'eps');

figure(4); clf; hold on;
%plot(Epsa1(:),Epsp1(:),'.k','MarkerSize',3);
plot(Epsa2X(:),Epsp2X(:),'xk','Color',[0.5, 0.5, 0.5],'MarkerSize',6);
scatter(Epsa1(:),Epsp1(:),[],-0.001*PeakVals1(:)./PeakVals1(:),'filled','s')
scatter(Epsa2(:),Epsp2(:),[],PeakVals2(:),'filled','s');
plot([0,0,0,-0.01],[0.295,0.355,0.415,0.355],'o','MarkerFaceColor', [0. 0. 0.],'MarkerEdgeColor', [0. 0. 0.],'linewidth',1,'MarkerSize',7)
for i=1:4
    if ~isempty(Br_tr{i})
        [perv, idx] = sort(arrayfun(@(x)x.parameter(ind.epsp),Br_tr{i}.point));
        parv = arrayfun(@(x)x.parameter(ind.epsa),Br_tr{i}.point);
        parv = parv(idx);
        plot(parv,perv,'-','Color', [0. 0. 0.],'linewidth',2);
    end
end
% for i=1:length(Br_pd)
%     if ~isempty(Br_pd{i})
%         parv = arrayfun(@(x)x.parameter(ind.epsa),Br_pd{i}.point);
%         perv = arrayfun(@(x)x.parameter(ind.epsp),Br_pd{i}.point);
%         plot(parv,perv,'--','Color', [0. 0. 0.],'linewidth',2);
%     end
% end
hold off;

bd = floor(max(PeakVals1(:))*512);
cmap = cmap0(1:bd,:);
colormap(cmap);
a = colorbar();
a.Label.String = '(max-min)/max in %';
title('Variability of peak total intensity - downsweap $\varepsilon_a$','interpreter','latex');
xlabel('$\varepsilon_a$','interpreter','latex');
ylabel('$\varepsilon_p$','interpreter','latex');
hgexport(gcf, 'Fig2colorbar.eps', hgexport('factorystyle'), 'Format', 'eps');

box on;
xAuxTicks = [-0.03, -0.01, 0.01];
yAuxTicks = [0.2, 0.4];
xlim([-0.035, 0.015]);
ylim([0.1, 0.5]);
hold off;
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,8,5])
set(gca,'position',[0.005,0.01,0.99,0.98],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',1.5) %[0.07,0.10,0.92,0.88]
xticklabels({'',''})
yticklabels({'',''})
hgexport(gcf, 'img/Fig2d.eps', hgexport('factorystyle'), 'Format', 'eps');
