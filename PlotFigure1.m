% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for producing panels (b) and (c) of figure 1
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["e06", "e05"];
opts.VariableTypes = ["double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Import the data and clear opts
%Data = readtable("/home/stefan/Dropbox/Work/Research/MATLAB/ddebiftool-git/projects/POLAR/ExpData/trc06-test.txt", opts);
Data = readtable("./expData/trc06-test.txt", opts);
clear opts

% Convert to array
Data = table2array(Data);
PeakData = load('./expData/PeakData.txt');

%% Call data
loc = find(isnan(Data(:,1)));

ty = Data(loc(1)+1:loc(2)-1,1);
Iy = Data(loc(1)+1:loc(2)-1,2);

tx = Data(loc(2)+1:loc(3)-1,1);
Ix = Data(loc(2)+1:loc(3)-1,2);

tp = PeakData(:,1);
Ip = PeakData(:,2);
mm = mean(PeakData(2:end,2));

%% Make plot

figure(1); clf; hold on;
plot(tx,Ix, 'LineWidth', 1., 'Color', [0, 0.4470, 0.7410]);
plot(ty,mm - Iy, 'LineWidth', 1., 'Color', [1.,0.5,0.]);
plot(tp,Ip,'*','MarkerSize',5,'MarkerEdgeColor','k');
xlim([3e-7,8e-7]);
ylim([-0.1,0.3]);
box on;

ticklengthUn = 0.1;
labelSize = 11;
box on;
xAuxTicks = [3e-7 5e-7 7e-7];
yAuxTicks = [0 0.2];
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,8,1.5])
set(gca,'position',[0.005,0.02,0.99,0.975],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',1.) %[0.07,0.10,0.92,0.88]
xticklabels({'','',''})
yticklabels({'',''})
%legend('C4','$\bar I$ -C2','interpreter','latex')
hgexport(gcf, 'img/Fig1b_raw.eps', hgexport('factorystyle'), 'Format', 'eps');

