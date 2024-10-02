% Plot slice data produced by BifDiag_A_*.py
close all; clear all;

labelSize =11;
ticklengthUn = 0.2;

DataDownX = load('pythonData/MaxPeaks_A0.355_tau2000.0_downX.txt');
DataDownY = load('pythonData/MaxPeaks_A0.355_tau2000.0_down.txt');
DataUpX = load('pythonData/MaxPeaks_A0.355_tau2000.0_upX.txt');
DataUpY = load('pythonData/MaxPeaks_A0.355_tau2000.0_up.txt');
ADown = DataDownX(:,1);
AUp = DataUpX(:,1);

%% Plot Downsweep data
figure(1); clf; hold on;
plot(ADown,DataDownX(:,2:end),'.','Color', [0 0.4470 0.7410],'MarkerSize',1);
plot(ADown,DataDownY(:,2:end),'.','Color', [0.8500 0.3250 0.0980],'MarkerSize',1);
hold off;


xlabel('$A$','interpreter','latex');
ylabel('$I$','interpreter','latex');

box on;
xAuxTicks = [2.1, 2.3, 2.5];
yAuxTicks = [0, 50, 100];
xlim([2.1, 2.5]);
ylim([-5, 115]);
hold off;
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,8/6*4,2.5])
set(gca,'position',[0.005,0.01,0.99,0.98],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',1.) %[0.07,0.10,0.92,0.88]
xticklabels({'',''})
yticklabels({'','',''})
hgexport(gcf, 'img/Fig4a.eps', hgexport('factorystyle'), 'Format', 'eps');

%% Plot Upsweep data
figure(2); clf; hold on;
plot(AUp,DataUpY(:,2:end),'.','Color', [0.8500 0.3250 0.0980],'MarkerSize',1);
plot(AUp,DataUpX(:,2:end),'.','Color', [0 0.4470 0.7410],'MarkerSize',1);
hold off;

xlabel('$A$','interpreter','latex');
ylabel('$I$','interpreter','latex');

box on;
xAuxTicks = [2.5, 2.7];
yAuxTicks = [0, 50, 100];
xlim([2.5, 2.7]);
ylim([-5, 115]);
hold off;
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,8/6*2,2.5])
set(gca,'position',[0.005,0.01,0.99,0.98],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',1.) %[0.07,0.10,0.92,0.88]
xticklabels({'',''})
yticklabels({'','',''})
hgexport(gcf, 'img/Fig4b.eps', hgexport('factorystyle'), 'Format', 'eps');

