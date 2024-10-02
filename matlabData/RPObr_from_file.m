%% create branch of relative psols from python raw data 
%
% Modified version of branch_from_sol. 
% Stefan Ruschel 19.10.2020
%
% inputs:
%
% * |nameFile|: address of nameFile_sol nameFile_par in folder pythonData
% * |funcs|: functions defining right-hand side (created with set_funcs)
% * |free_par_ind|: index of parameter to be varied along branch initially
%
% Important optional inputs (name-value pairs, all passed on to )
%
% * |'degree'|: (default 3) degree used for collocation
% * |'intervals'|: (default 100) number of collocation intervals (overall mesh size is
% degree x intervals + 1
% * |'stepcond'|: (default []) passed on to p_correc
% * |'step'|: (default 1e-2) initial step between first two points
% * |'corpar'|: parameters left free for initial correction (if different
% from |free_par_ind|)

function branch=RPObr_from_file(nameFile,rfuncs,free_par_ind,cor_par_ind, period_ind, varargin)

default ={'degree',3,'intervals',301,'corpar',cor_par_ind,'step',1e-3, 'indperiod', [], 'stepcond',[]};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');

% set up branch 
branch=df_brnch(rfuncs,free_par_ind,'psol');
branch=replace_branch_pars(branch,free_par_ind,pass_on);
corpar=union(options.corpar,options.indperiod);

% load solution and parameters from python data
solPy   =  load([nameFile '_sol.txt']);
parPy   =  load([nameFile '_par.txt']);
dim     =  length(solPy(1,:));

% create point struct, remesh, correct, 2nd point and set points
pt      =  dde_psol_create('parameter',parPy.','mesh',solPy(:,1).','degree',1,...
    'profile',solPy(:,(2:dim)).','period',parPy(period_ind));
pt      =  p_remesh(pt,options.degree,options.intervals);
method  =  df_mthod(rfuncs,'psol');
[pt1,suc]=   p_correc(rfuncs,pt,corpar,[],branch.method.point,0,pt);
pt2     =  pt1;
pt2.parameter(free_par_ind(1))=pt2.parameter(free_par_ind(1))+options.step;
[pt2,suc]     =   p_correc(rfuncs,pt2,corpar,[],branch.method.point,0,pt);
branch.point=[pt1,pt2];
end

