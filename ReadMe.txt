This folder contains code to reproduce all tables and figures in
"The Elliptical Ornstein-Uhlenbeck Process"
Sykulski A.M., Olhede S.C., Sykulska-Lawrence H. (2021)
To appear in Statistics and Its Interface

This folder contains 15 additional files:

Data:
- Earthwobble.mat contains the polar motion data studied

Core files:
- Fig1LeftColumn.m regenerates Figure 1 Left
- Fig1RightColumn.m regenerates Figure 1 Right
- Fig2andTable2.m regenerates Figure 2 and Table 2
- Figs3and4.m regenerates Figures 3 and 4
- Figs5and6andTable4.m regenerates Figures 5 and 6 and Table 4
- Table3.m regenerates Table 3

Supporting files:
- WILCOUmodelFull2.m is the Whittle likelihood from equation (17) for the elliptical OU
- WILCOUmodelFullF.m is the Whittle likelihood from equation (15) for the elliptical OU
- WILCOUmodelRange2.m is the Whittle likelihood from equation (17) for the elliptical OU over a reduced Frequency range
- WILCOUmodelRangeF.m is the Whittle likelihood from equation (15) for the elliptical OU over a reduced Frequency range
- WILCOUmodelRangeP is the Whittle likelihood for equation (17) for the circular OU
- fminsearchbnd.m is an extension to fminsearch which allows parameters to be bounded by transformations in such a way that unconstrained optimisation via fminsearch can still be performed
- boundedline.m produces the confidence regions for the figure
- kde.m produces kernel density estimates used in Figure 2

Instructions:
Run Core files with all accompanying files in the same directory