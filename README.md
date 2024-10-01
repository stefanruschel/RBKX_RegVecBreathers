# RBKX_RegVecBreathers

Source code by Stefan Ruschel to generate the data and figures for the paper:

	Regenerative vectorial breathers in a delay-coupled neuromorphic microlaser with integrated saturable absorber (https://arxiv.org/abs/2409.20177)
	by Stefan Ruschel, Venkata A. Pammi, Rémy Braive, Isabelle Sagnes, Grégoire Beaudoin, Neil G. R. Broderick, Bernd Krauskopf, Sylvain Barbay

Folder expData contains experimentally obtained data. Folder pythonData contains time series of a corresponding mathametical model which are computed numerically using python package pydelay v 0.1.1 (https://pydelay.sourceforge.net/). Folder matlabData contains source code for numerical continuation of torus bifurcation curves; it is written for Matlab2020a, and uses the latest version of the software package "DDE-Biftool" (https://github.com/DDE-BifTool/DDE-Biftool).

In order to reproduce numerical simulation and pathfollowing results please execute: 

	./pythonData/Figure3*.py: these 4 scripts simulate the model to obtain time series shown in Figure 3.

Currently, all the data must be computed once due to the size constraints of GitHub. Please change the scripts to load the correct address where your installation of DDE-Biftool is saved.

To produce panels of Figures 1-4 please run:

	PlotFigure1.m: Loads experimental data and data for main peak intensity and plots panel (b) of Fig. 1.
	PlotFigure2.m: Loads maximum data for up and down sweep together with data for the torus curves and plots Fig.2(a)-(d). (TODO!!!!!!!!!!!!!!!)
 	PlotFigure3.m: Load numerically computed time series data and plots Fig.3(a)-(d). 
 	PlotFigure4.m: Loads maximum data for up and down seep from A=2.5 and plots Fig.4(a)-(b). (TODO!!!!!!!!!!!!!!!)

Please cite https://arxiv.org/abs/2409.20177 if you use (part of) this code in your research.


 

 
