%% Load DDE-Biftool and extension into Path
clear all; close all;
addpath('../../DDE-Biftool/ddebiftool/',...
    '../../DDE-Biftool/ddebiftool_extra_psol/',...
    '../../DDE-Biftool/ddebiftool_utilities/',...
    '../../DDE-Biftool/ddebiftool_extra_rotsym');

% Problem definition using |set_rotfuncs|
A=[0,-1,0,0,0,0,0,0; 1,0,0,0,0,0,0,0; 0,0,0,-1,0,0,0,0; 0,0,1,0,0,0,0,0; ...
   0,0,0,0,0,0,0,0;  0,0,0,0,0,0,0,0; 0,0,0,0,0,0,0,0; 0,0,0,0,0,0,0,0];
expA=@(phi) [cos(phi), -sin(phi),0,0,0,0,0,0; sin(phi), cos(phi),0,0,0,0,0,0; ...
             0,0,cos(phi), -sin(phi),0,0,0,0; 0,0, sin(phi), cos(phi),0,0,0,0; ...
             0,0,0,0,1,0,0,0; 0,0,0,0,0,1,0,0; 0,0,0,0,0,0,1,0; 0,0,0,0,0,0,0,1];
% Initial values of parameters and parameter indices
parnames={...
    'alpha',...   % alpha factor
    'beta',...    % alpha factor for losses
    'epsa',...    % amplitude anisoptropy
    'epsp',...    % phase anisotropy
    'cp',...      % coupling strength
    'phi',...     % feedback phase
    'tau',...     % feedback delay
    'gammaG',...  % gain relaxation rate
    'A',...       % pump strength
    'gammaGSF',...% gain spin flip rate
    'gammaQ',...  % gain relaxation rate
    'B',...       % 
    'a',...       % ratio gammaG/gammaQ*sigma
    'gammaQSF',...% losses' spin flip rate
    'T',...       % placeholder for period
    'omega'};     % rotation velocity
cs=[parnames;num2cell(1:length(parnames))];
ind=struct(cs{:});
par=zeros(1,numel(fieldnames(ind)));
% placeholder for parameter values
par([ind.alpha,ind.beta,ind.epsa,ind.epsp,ind.cp,ind.phi,ind.tau,ind.gammaG])=...
    [       2.,      0.5,     0.0,     0.1,    .2,      0.,  100.,      0.01];
par([ind.A,ind.gammaGSF,ind.gammaQ,ind.B,ind.a,ind.gammaQSF,ind.T,ind.omega])=...
    [  2.5,     10*0.01,      0.01,  2.0,  10.,     10*0.01,   0.,       0.];

% Right-hand side and call to |set_rotfuncs|
f=@(x,p)SFMSA(x(1,1,:)+1i*x(2,1,:),x(1,2,:)+1i*x(2,2,:),...
              x(3,1,:)+1i*x(4,1,:),x(3,2,:)+1i*x(4,2,:),...
              x(5,1,:),x(6,1,:),x(7,1,:),x(8,1,:),p,ind);
rfuncs=set_rotfuncs('sys_rhs',f,'rotation',A,'exp_rotation',expA,...
    'sys_tau',@()ind.tau,'x_vectorized',true);
bd={'max_bound',[],'min_bound',[]};
opt_inputs=[bd,{'extra_condition',1,'print_residual_info',1},...
    'matrix','sparse','eigmatrix','sparse']; 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Set up branch from file                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nameFile = '../pythonData/PO_p1_l0_epsa_00_epsp01';
k = 1; % branch identifier. k=1:PO_epsa_001_epsp01_sol
       
% single polarization pulses are badly scaled in general
Br_psol{k}=RPObr_from_file(nameFile,rfuncs,[ind.tau ind.omega],ind.omega,ind.T, 'step', 1e-3, 'degree', 10, 'intervals', 101, opt_inputs{:});
Br_psol{k}.parameter.max_step=[0 100.; ind.epsa 0.01; ind.omega 0.01];
Br_psol{k}.parameter.max_bound=[ind.tau, 2000];
% set up from existing branch
% Br_psol{1}=RPObr_from_point(Br_psol{1},260,rfuncs,[ind.tau ind.omega],ind.tau);
        %% continue branch
        k=1;
        figure(2);ax=gca; hold on;
        Br_psol{k}=br_contn(rfuncs,Br_psol{k},2000,'plotaxis',ax);
        hold off;
        %% switch to continuation in epsa and continue branch
        k=2; 
        m=1; % setup from Br_psol{m}
        Br_psol{k} = SetupPsolFrom_psol(rfuncs,Br_psol{m},length(Br_psol{m}.point),...
                        'contpar',[ind.epsa ind.omega]);
        Br_psol{k}.parameter.max_step=[0 1.; ind.epsa 0.005; ind.omega 0.01];            
        Br_psol{k}.parameter.max_bound=[ind.epsa, 0.01];
        Br_psol{k}.parameter.min_bound=[ind.epsa, -0.03];
        figure(3);ax=gca; hold on;
        Br_psol{k}=br_contn(rfuncs,Br_psol{k},1000,'plotaxis',ax);
        Br_psol{k}.point(1:5)=[];
        Br_psol{k}=br_rvers(Br_psol{k});
        Br_psol{k}=br_contn(rfuncs,Br_psol{k},1000,'plotaxis',ax);
        hold off;            
        %% compute stability
        k=2;
        [Br_psol_nunst{k},Br_psol_dom{k},Br_psol_defect{k},Br_psol{k}.point]=GetStability(Br_psol{k},...
        'exclude_trivial',true,'locate_trivial',@(p)[1,1],'funcs',rfuncs);
        beep;
        
        %% switch to continuation in epsp and continue branch tp find PD
        k=3; 
        m=1; % setup from Br_psol{m}
        Br_psol{k} = SetupPsolFrom_psol(rfuncs,Br_psol{m},length(Br_psol{m}.point),...
                        'contpar',[ind.epsp ind.omega]);
        Br_psol{k}.parameter.max_step=[0 1.; ind.epsp 0.005; ind.omega 0.001];            
        Br_psol{k}.parameter.max_bound=[ind.epsp, 0.5];
        Br_psol{k}.parameter.min_bound=[ind.epsp, 0.08];
        figure(3);ax=gca; hold on;
        Br_psol{k}=br_contn(rfuncs,Br_psol{k},1000,'plotaxis',ax);
        Br_psol{k}.point(1:5)=[];
        Br_psol{k}=br_rvers(Br_psol{k});
        Br_psol{k}=br_contn(rfuncs,Br_psol{k},1000,'plotaxis',ax);
        hold off;            
        %% compute stability
        k=3;
        [Br_psol_nunst{k},Br_psol_dom{k},Br_psol_defect{k},Br_psol{k}.point]=GetStability(Br_psol{k},...
        'exclude_trivial',true,'locate_trivial',@(p)[1,1],'funcs',rfuncs);
        beep;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Set up branch from file                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nameFile = './pythonData/PO_p1_lpi_epsa_m003_epsp01';
k = 4; % branch identifier. k=4:PO_p1_lpi_epsa_m003_epsp01
       
% single polarization pulses are badly scaled in general
Br_psol{k}=RPObr_from_file(nameFile,rfuncs,[ind.tau ind.omega],ind.omega,ind.T, 'step', 1e-3, 'degree', 10, 'intervals', 100, opt_inputs{:});
Br_psol{k}.parameter.max_step=[0 100.; ind.epsa 0.01; ind.omega 0.01];
Br_psol{k}.parameter.max_bound=[ind.tau, 2000];
% set up from existing branch
        %% continue branch
        k=4;
        figure(2);ax=gca; hold on;
        Br_psol{k}=br_contn(rfuncs,Br_psol{k},2000,'plotaxis',ax);
        hold off;
        %% switch to continuation in epsa and continue branch
        k=5; 
        m=4; % setup from Br_psol{m}
        Br_psol{k} = SetupPsolFrom_psol(rfuncs,Br_psol{m},length(Br_psol{m}.point),...
                        'contpar',[ind.epsa ind.omega]);
        Br_psol{k}.parameter.max_step=[0 1.; ind.epsa 0.005; ind.omega 0.01];            
        Br_psol{k}.parameter.max_bound=[ind.epsa, 0.01];
        Br_psol{k}.parameter.min_bound=[ind.epsa, -0.031];
        figure(3);ax=gca; hold on;
        Br_psol{k}=br_contn(rfuncs,Br_psol{k},1000,'plotaxis',ax);
        Br_psol{k}.point(1:5)=[];
        Br_psol{k}=br_rvers(Br_psol{k});
        Br_psol{k}=br_contn(rfuncs,Br_psol{k},1000,'plotaxis',ax);
        hold off;            
        %% compute stability
        k=5;
        [Br_psol_nunst{k},Br_psol_dom{k},Br_psol_defect{k},Br_psol{k}.point]=GetStability(Br_psol{k},...
        'exclude_trivial',true,'locate_trivial',@(p)[1,1],'funcs',rfuncs);
        ind_tr=find(abs(diff(Br_psol_nunst{k}))==2);
        %% switch to continuation in epsp and continue branch
        k=6; 
        m=5; % setup from Br_psol{m}
        Br_psol{k} = SetupPsolFrom_psol(rfuncs,Br_psol{m},7,...
                        'contpar',[ind.epsp ind.omega]);
        Br_psol{k}.parameter.max_step=[0 1.; ind.epsp 0.01; ind.omega 0.01];            
        Br_psol{k}.parameter.max_bound=[ind.epsa, 0.01; ind.epsp, 0.52];
        Br_psol{k}.parameter.min_bound=[ind.epsa, -0.03; ind.epsp, 0.1];
        figure(3);ax=gca; hold on;
        Br_psol{k}=br_contn(rfuncs,Br_psol{k},1000,'plotaxis',ax);
        Br_psol{k}.point(1:5)=[];
        Br_psol{k}=br_rvers(Br_psol{k});
        Br_psol{k}=br_contn(rfuncs,Br_psol{k},1000,'plotaxis',ax);
        hold off;
        %% compute stability
        k=6;
        [Br_psol_nunst{k},Br_psol_dom{k},Br_psol_defect{k},Br_psol{k}.point]=GetStability(Br_psol{k},...
        'exclude_trivial',true,'locate_trivial',@(p)[1,1],'funcs',rfuncs);
        
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Torus of relative POs
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Br_tr{k,l}=cell(1,1); Br_tr_funcs{k,l}=cell(1,1);
        % Br_tr{k,l},   k=2, l=1: Br_psol{2}closeToZero    
        %               k=3, l=1: same TB from Br_psol{3} 
        %               k=5, l=1: TB from Br_psol{5} 
        %               k=6, l=1: same TB from Br_psol{5} but higher epsp
        k=6; l=1;
        ind_tr=find(abs(diff(Br_psol_nunst{k}))==2);
        [Br_tr_funcs{k,l},Br_tr{k,l}]=SetupMWTorusBifurcation(rfuncs,Br_psol{k},ind_tr(1),...
            'contpar',[ind.epsa,ind.epsp,ind.omega],'step',1e-3,opt_inputs{:},...
            'print_residual_info',1,'dir',ind.epsa,'newton_max_iterations',10);
        Br_tr{k,l}.parameter.max_step=[0 10.; ind.epsa 0.001; ind.epsp, 0.01; ind.omega 0.001];            
        Br_tr{k,l}.parameter.max_bound=[ind.epsa, 0.01; ind.epsp, 0.52];
        Br_tr{k,l}.parameter.min_bound=[ind.epsa, -0.03; ind.epsp, 0.1];
        %% continuation
        k=2; l=1;
        figure(20); clf; hold on;
        Br_tr{k,l}=br_rvers(Br_tr{k,l});
        Br_tr{k,l}=br_contn(Br_tr_funcs{k,l},Br_tr{k,l},200);
        hold off;
        beep
        %% continuation
        k=3; l=1;
        figure(20); clf; hold on;
        Br_tr{k,l}=br_contn(Br_tr_funcs{k,l},Br_tr{k,l},50);
        Br_tr{k,l}=br_rvers(Br_tr{k,l});
        Br_tr{k,l}=br_contn(Br_tr_funcs{k,l},Br_tr{k,l},50);
        hold off;
        beep
        %% continuation
        k=5; l=1;
        figure(20); clf; hold on;
        Br_tr{k,l}=br_contn(Br_tr_funcs{k,l},Br_tr{k,l},100);
        Br_tr{k,l}=br_rvers(Br_tr{k,l});
        Br_tr{k,l}=br_contn(Br_tr_funcs{k,l},Br_tr{k,l},100);
        hold off;
        beep
        %% continuation
        k=6; l=1;
        figure(20); clf; hold on;
        Br_tr{k,l}=br_contn(Br_tr_funcs{k,l},Br_tr{k,l},100);
        Br_tr{k,l}=br_rvers(Br_tr{k,l});
        Br_tr{k,l}=br_contn(Br_tr_funcs{k,l},Br_tr{k,l},100);
        hold off;
        beep
        
%% Quick plotting        
figure(6); clf; hold on;
for i=1:length(Br_tr)
    if ~isempty(Br_tr{i})
        parv = arrayfun(@(x)x.parameter(ind.epsa),Br_tr{i}.point);
        perv = arrayfun(@(x)x.parameter(ind.epsp),Br_tr{i}.point);
        plot(parv,perv,'-o','Color', [0. 1. 0.],'linewidth',2);
    end
end
xlabel('$\varepsilon_a$','interpreter','latex');
ylabel('$\varepsilon_p$','interpreter','latex');

%% Save data
save('TorusCurves.mat')

