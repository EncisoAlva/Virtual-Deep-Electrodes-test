surface = inner_skull;

nPA = 4;
nLR = 1;
PA0 = 10/1000; % 10 mm
LR0 = 10/1000;
dPA = 10/1000;
dLR = 10/1000;

TH = 0;

labels = {'L', 'R'};

params.Side       = 'Both';
params.vis_width  = 0.001;
params.vis_tol    = params.vis_width/20;
params.max_iter   = 100;
params.cutZ1      = 10;
params.Tol1       = 0.0175;

debug.figs        = true;
debug.figs_detail = false;

[ElectrodeLocations] = PlaceSurfaceElectrodesIdeal( ...
  surface, nPA, nLR, PA0, LR0, dPA, dLR, TH, labels, ...
  params, debug );

%%
surface = inner_skull;

nIS = 8;
PA0 = 25/1000; % 10 +10 +10/2 mm
LR0 = 10/1000;
IS0 =  2/1000;
dIS =  5/1000;

TH = 0;

labels = {'LD', 'RD'};

params.Side       = 'Both';
params.vis_width  = 0.001;
params.vis_tol    = params.vis_width/20;
params.max_iter   = 100;
params.cutZ1      = 10;
params.Tol1       = 0.0175;

debug.figs        = true;
debug.figs_detail = false;

PlaceInsertedElectrodesIdeal( ...
  surface, nIS, PA0, LR0, IS0, dIS, TH, labels, ...
  params, debug )