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



none=0;