% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for producing panels (a)--(d) of figure 3
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot Figure a 
Data = load('./pythonData/Figure3a_sol.txt');
t = Data(:,1);
Ep = Data(:,2)+1j.*Data(:,3);
Em = Data(:,4)+1j.*Data(:,5);

%%
Ix = abs(Ep + Em).^2/2;
Iy = abs(Ep - Em).^2/2;
I  = abs(Ep).^2+abs(Em).^2; % total intensity

id = find(islocalmax(I) & I>10);
idx = find(islocalmax(Ix) & Ix>10);
idy = find(islocalmax(Iy) & Iy>10);
mm = mean(Ix(idx))+mean(Iy(idy));
mx = max(I);
mxinterp1 = interp1(t(id),I(id),t);

figure(2); clf; hold on;
plot(t-t(1),Ix, 'LineWidth', 1., 'Color', [0, 0.4470, 0.7410]);
plot(t-t(1),mm-Iy, 'LineWidth', 1., 'Color', [1.,0.5,0.]);
plot(t(id)-t(1),I(id),'*','MarkerSize',5,'MarkerEdgeColor','k');
hold off;
xlim([0,50]);
ylim([0,100]);

ticklengthUn = 0.1;
labelSize = 11;
box on;
xAuxTicks = [0 20 40];
yAuxTicks = [0 50];
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,8,1.2])
set(gca,'position',[0.005,0.02,0.99,0.975],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',1.) %[0.07,0.10,0.92,0.88]
xticklabels({'','','',''})
yticklabels({'','','',''})
hgexport(gcf, 'img/Fig3a_raw.eps', hgexport('factorystyle'), 'Format', 'eps');

%% Plot Figure b 
Data = load('./pythonData/Figure3b_sol.txt');
t = Data(:,1);
Ep = Data(:,2)+1j.*Data(:,3);
Em = Data(:,4)+1j.*Data(:,5);

%%
Ix = abs(Ep + Em).^2/2;
Iy = abs(Ep - Em).^2/2;
I  = abs(Ep).^2+abs(Em).^2; % total intensity

id = find(islocalmax(I) & I>10);
idx = find(islocalmax(Ix) & Ix>10);
idy = find(islocalmax(Iy) & Iy>10);
mm = mean(Ix(idx))+mean(Iy(idy));
mx = max(I);
mxinterp1 = interp1(t(id),I(id),t);

figure(2); clf; hold on;
plot(t-t(1),Ix, 'LineWidth', 1., 'Color', [0, 0.4470, 0.7410]);
plot(t-t(1),mm-Iy, 'LineWidth', 1., 'Color', [1.,0.5,0.]);
plot(t(id)-t(1),I(id),'*','MarkerSize',5,'MarkerEdgeColor','k');
hold off;
xlim([0,50]);
ylim([0,100]);

ticklengthUn = 0.1;
labelSize = 11;
box on;
xAuxTicks = [0 20 40];
yAuxTicks = [0 50];
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,8,1.2])
set(gca,'position',[0.005,0.02,0.99,0.975],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',1.) %[0.07,0.10,0.92,0.88]
xticklabels({'','','',''})
yticklabels({'','','',''})
hgexport(gcf, 'img/Fig2b_raw.eps', hgexport('factorystyle'), 'Format', 'eps');

%% Plot Figure c 
Data = load('./pythonData/Figure3c_sol.txt');
t = Data(:,1);
Ep = Data(:,2)+1j.*Data(:,3);
Em = Data(:,4)+1j.*Data(:,5);

%%
Ix = abs(Ep + Em).^2/2;
Iy = abs(Ep - Em).^2/2;
I  = abs(Ep).^2+abs(Em).^2; % total intensity

id = find(islocalmax(I) & I>10);
idx = find(islocalmax(Ix) & Ix>10);
idy = find(islocalmax(Iy) & Iy>10);
mm = mean(Ix(idx))+mean(Iy(idy));
mx = max(I);
mxinterp1 = interp1(t(id),I(id),t);

figure(2); clf; hold on;
plot(t-t(1),Ix, 'LineWidth', 1., 'Color', [0, 0.4470, 0.7410]);
plot(t-t(1),mm-Iy, 'LineWidth', 1., 'Color', [1.,0.5,0.]);
plot(t(id)-t(1),I(id),'*','MarkerSize',5,'MarkerEdgeColor','k');
hold off;
xlim([0,50]);
ylim([0,100]);

ticklengthUn = 0.1;
labelSize = 11;
box on;
xAuxTicks = [0 20 40];
yAuxTicks = [0 50];
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,8,1.2])
set(gca,'position',[0.005,0.02,0.99,0.975],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',1.) %[0.07,0.10,0.92,0.88]
xticklabels({'','','',''})
yticklabels({'','','',''})
hgexport(gcf, 'img/Fig3c_raw.eps', hgexport('factorystyle'), 'Format', 'eps');

%% Plot Figure d 
Data = load('./pythonData/Figure3d_sol.txt');
t = Data(:,1);
Ep = Data(:,2)+1j.*Data(:,3);
Em = Data(:,4)+1j.*Data(:,5);

%%
Ix = abs(Ep + Em).^2/2;
Iy = abs(Ep - Em).^2/2;
I  = abs(Ep).^2+abs(Em).^2; % total intensity

id = find(islocalmax(I) & I>10);
idx = find(islocalmax(Ix) & Ix>10);
idy = find(islocalmax(Iy) & Iy>10);
mm = mean(I(id));
mx = max(I);
mxinterp1 = interp1(t(id),I(id),t);

figure(2); clf; hold on;
plot(t-t(1),Ix, 'LineWidth', 1., 'Color', [0, 0.4470, 0.7410]);
plot(t-t(1),mm-Iy, 'LineWidth', 1., 'Color', [1.,0.5,0.]);
plot(t(id)-t(1),I(id),'*','MarkerSize',5,'MarkerEdgeColor','k');
hold off;
xlim([0,50]);
ylim([0,100]);

ticklengthUn = 0.1;
labelSize = 11;
box on;
xAuxTicks = [0 20 40];
yAuxTicks = [0 50 100];
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,8,1.2])
set(gca,'position',[0.005,0.02,0.99,0.975],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',1.) %[0.07,0.10,0.92,0.88]
xticklabels({'','','',''})
yticklabels({'','','',''})
hgexport(gcf, 'img/Fig3d_raw.eps', hgexport('factorystyle'), 'Format', 'eps');

