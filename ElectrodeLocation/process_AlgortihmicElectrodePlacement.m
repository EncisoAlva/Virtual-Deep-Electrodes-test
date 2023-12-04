function none = process_AlgortihmicElectrodePlacement( vars )
%placekeeper function
nPA = 4;
nLR = 1;
PA0 = 10; % 10 mm
LR0 = 10;
dPA = 10;
dLR = 10;

params.Side      = 'Left';
params.vis_width = 0.001;
params.vis_tol   = params.vis_width/20;
params.max_iter  = 100;
params.cutZ1     = 10;
params.Tol1      = 0.0075;

none=0;
end

%% SURFACE ELECTRODES
function [ElecLocs] = PlaceSurfaceElectrodesIdeal( ...
  surface, nPA, nLR, PA0, LR0, dPA, dLR, TH, ...
  params, debug )
% This script determine locations of an NxM surface electordes following 
% this simple protocol:
%  1. Identify the Central Line: intersection of the plane between brain 
%     hempispheres and the envelope of brain cortex.
%  2. Identify Anterior Edge: the most anterior points of the brain, 
%     approximately ortogonal to both the Central Line and the 
%     Posterior-Anterior line.
%  3. Set the first electrode at 
%       PA0 cm from Anterior Edge on the Posterior-Anterior direction,
%       LR0 cm from Central Line on the Left-Right direction.
%  4. Set the other electrodes following a rectangular grid with one side 
%     at a slope angle of TH degrees with respect to the Central Line.
%     The grid has center-to-center distances of dPA cm in the 
%     Posterior-Anterior direction, and dLR cm in the Left-Right direction. 
%     The grid is assumed to follow the curvature of the cortex envelope.
% 
% This function makes assumptions of ideal placing of electrodes.
%
%-------------------------------------------------------------------------
% INPUT
%    V  surface  Surface triangulation with these fields
%    |.Vertices  Location of triangles' edges, (# vertices)x3
%    |   .Faces  Triangulation, (# triangles)x3
%           nPA  Size of rectangle grid in the Posterior-Anterior
%           nLR  Size of rectangle grid in the Left-Right direction
%           PA0  Distance [mm] from Anterior Edge to the first electrode
%           LR0  Distance [mm] from Central Line to the first electrode
%           dPA  Center-to-center distance between electrodes in the
%                Posterior-Anterior direction
%           dLR  Center-to-center distance between electrodes in the
%                Left-Rigt direction
%      V params  Parameters
%      |  .Side  Which hemisphere: Left, Right, Both
%      | .Strip  Width of 
%
%-------------------------------------------------------------------------
% OUTPUT
%      ElecLocs  Locations of electrodes, (nAP)x(nLR)x3
%
%-------------------------------------------------------------------------
% Author: Julio Cesar Enciso-Alva, 2023
%         juliocesar.encisoalva@uta.edu
%

% rename of parameters
coord = surface.Vertices;

switch params.Side
  case 'Left'
    ElectSide = {'Left'};
  case 'Right'
    ElectSide = {'Right'};
  case 'Both'
    ElectSide = {'Left','Right'};
  otherwise
    ElectSide = {};
    disp('WARNING: Select a side to place electrodes.')
end

ElecLocs      = zeros(nLR, nPA, 3); % cntainer for results

% notes to self of coordinates
%  1  x  Posterior-Anterior (reverse Anterior-Posterior)
%  2  y  Left-Right
%  3  z  Inferior-Superior
dir.PA = 1;
dir.LR = 2;
dir.IS = 3;

%% IDENTIFICATION OF ANTERIOR EDGE
% defined by it's visual location, computed iteratively for robustness
iter = 1; flag=true; ptEDGE_old = [Inf,Inf,Inf];
while (iter < params.MaxIter) && flag
  if iter == 1
    [maxPA,~] = max( coord(:,dir.PA) );
    %[maxPA,idxPA] = max( coord(:,dir.PA) );
    %ptPA = coord(idxPA,:);
    region = coord( abs(coord(:,dir.PA)-maxPA)<params.Tol1,:);
    [~,idxIS] = max( region(:,dir.IS) );
    ptIS = region(idxIS,:);
    [~,idxPA] = min( region(:,dir.IS) );
    ptPA = region(idxPA,:);
  else
    region = coord( vecnorm( ptEDGE([dir.PA,dir.IS])-coord(:,[dir.PA,dir.IS]),2,2 )<params.Tol1,: );
    [~,idxPA] = max( region(:,dir.PA) );
    ptPA = region(idxPA,:);
    [~,idxIS] = max( region(:,dir.IS) );
    ptIS = region(idxIS,:);
  end
  % dist from (X,Y) to the line passin through (x1,y1) and (x2,y2)
  %   |(x2-x1)(y1-Y) - (x1-X)(y2-y1)|/norm( (x1,y1) - (x2,y2) )
  x1 = ptPA(dir.PA);  y1 = ptPA(dir.IS);
  x2 = ptIS(dir.PA);  y2 = ptIS(dir.IS);
  [~,idxEDGE] = max( (x2-x1)*(y1-coord(:,dir.IS)) - (y2-y1)*(x1-coord(:,dir.PA)) );
  ptEDGE = coord(idxEDGE,:);
  if isequal( ptEDGE, ptEDGE_old)
    flag = false;
  else
    ptEDGE_old = ptEDGE;
    iter = iter + 1;
  end
end
if debug.figs
  figure()
  scatter(coord(:,1),coord(:,3),'filled')
  xlabel('x')
  ylabel('z')
  title('Cortex Surface (Mesh Vertices)')
  subtitle('Sagittal View')
  hold on 
  scatter(coord(idxEDGE,1),coord(idxEDGE,3),200,'r','filled')
  legend({[],'Anterior Edge'})
  scatter(ptPA(1),ptPA(3),200,'k','filled')
  scatter(ptIS(1),ptIS(3),200,'k','filled')
end

coord = coord( coord(:,dir.IS) >= ptEDGE(dir.IS) ,:); % only upper part is relevant

%% LOCATION OF ELECTRODES
% 1. draw medial line, for reference
% 2. mark B regular dots in medial line
% 3. draw lines parallel to medial line
% 4. mark A regular dots on each one
% 5. TODO correction for robustness

for lazy_counter = 1:length(ElectSide)
CurrentSide = ElectSide{lazy_counter};
switch CurrentSide
  case 'Left'
    dLR = -abs(dLR);
    TH  = -abs(TH );
    sg  = -1;
  case 'Rigth'
    dLR =  abs(dLR);
    TH  =  abs(TH );
    sg  =  1;
end

% find the first point by first looking at the PA distance and ten the 
CtrlStrip = coord( coord(:,dis.LR) < params.vis_width, : );
CtrlLine  = SmoothCurveInterpolation( coord, ...
  min( CtrlStrip(:,dir.PA) ), 0, 180, params.vis_width, 0.0001 );

tmpPt = PointsInCurve( CtrlLine, PA0, dPA, 1 );
auxLR = SmoothCurveInterpolation( coord, ...
  tmpPt(dir.PA), tmpPt(dir.PA), sg*90, params.vis_width, 0.0001 );

FirstPt = PointsInCurve( auxLR, LR0, dLR, 1 );
ElecLocs(1,1,:) = FirstPt;

% iterate over Posterior-Anterio, then over Left-Right
auxPA = SmoothCurveInterpolation( coord, ...
  FirstPt(dir.PA), FirstPt(dir.LR), TH, params.vis_width, 0.0001 );
ElecLocs(1:nPA,1,:) = PointsInCurve( auxPA, 0, dPA, nPA );

for ii = 1:nPA
  refPt = ElecLocs(ii,1,:);
  auxLR = SmoothCurveInterpolation( coord, ...
    refPt(dir.PA), refPt(dir.LR), TH+sg*90, params.vis_width, 0.0001 );
  ElecLocs(ii,1:nLR,:) = PointsInCurve( auxLR, 0, dLR, nLR );
end

end

if(debug.figs)
    figure()
    trisurf(surface.Faces, coord(:,1),coord(:,2),coord(:,3),...
        'FaceColor',[.85 .85 .85],'FaceAlpha',.4)
    xlabel('Posterior-Anterior')
    ylabel('Left-Right')
    zlabel('Inferior-Superior')
    title('Cortex surface')
    hold on
    for ii = 1:nPA
      for jj = 1:nLR
        scatter3(ElecLocs(ii,jj,1),ElecLocs(ii,jj,2),ElecLocs(ii,jj,3),'filled')
      end
    end
end

writematrix(ElecLocs)
end

%% INSERTED ELECTRODES

% This script determine locations of N electordes on a stylet following 
% this simple protocol:
%  1. Identify the Central Line: intersection of the plane between brain 
%     hempispheres and the envelope of brain cortex.
%  2. Identify Anterior Edge: the most anterior points of the brain, 
%     approximately ortogonal to both the Central Line and the 
%     Posterior-Anterior line.
%  3. Set the entry point for the stylet at
%       P cm from Anterior Edge on the Anterior-Posterior direction
%       R cm from Central Line on the Left-Right direction
%  4. Define the stylet trajectory parallel to the Dorsal-Ventral
%     direction, with a deviation of 
%  4. Set the other electrodes following the grid, with center-to-center
%     distances of dP cm in the Anterior-Posterior direction, and dR cm in
%     the Left-Right direction.
% 
%-------------------------------------------------------------------------
% INPUT
%         W  Wiener Kernel, GPxSE
%         Y  Recordings from surface electrodes, SExTT
%      time  Vector of timestamps, 1xTT
%      dist  Max distance dipole-to-electrode in meters, 1x1
%   ElecPos  Positions of intended depth electrodes, DEx3
%   GridPos  Positions of distributed dipoles, GPx3
%
%-------------------------------------------------------------------------
% INPUT (OPTIONAL)
%    constr  Constrained optimization, default = F
%            Fixed orientation (1 column/dipole) or not (3 col/dipole)
%     magn   Magnitude as final result, default = T
%            If  constrained: T=magnitude, F=PCA 1st coordinate
%            If ~constrained: T=abs val, F=no changes
%    debug   If true, and if constrained optimization, give also magnitudes
%
%-------------------------------------------------------------------------
% OUTPUT
%    EstChans  Estimated channels, DExTT
% smallKernel  ...
%
%-------------------------------------------------------------------------
% Author: Julio Cesar Enciso-Alva, 2023
%         juliocesar.encisoalva@uta.edu
%


%% ADDITIONAL FUNCTIONS
function [Curve] = SmoothCurveInterpolation( ...
  surf, x0,y0, theta, w, dt )
% 1. Create a line passing trough a given point (x0,y0) with an angle (th)
%    with stps of (dt) mm.
% 2. Select the points on the given 3D triangulated surface whose 
%    projection is within (w) mm of the line.
% 3. Interpolate the a smooth curve over the surface, and evaluate over the
%    line.
% 
%-------------------------------------------------------------------------
% INPUT
%      surf   Triangulation vertices of 3d surface, (# vertices)x3.
%    (x0,y0)  Starting point.
%     theta   Angle with respect to the y-axis, in degrees.
%         w   Half width of the strip to be taken from the surface.
%        dt   Point-to-point distance on the smoothened line.
%
%-------------------------------------------------------------------------
% OUTPUT
%     Curve  Smoothened curve interpolating the strip obtained from the
%            surface.
%
%-------------------------------------------------------------------------
% Author: Julio Cesar Enciso-Alva, 2023
%         juliocesar.encisoalva@uta.edu
%

% conversion to radians
th = (90-theta)*pi/180;

% lazy variables
x=1; y=2; %z=3;

% dist from (X,Y) to line passing trough (x0,y0) at angle th
%   | cos(th)(y0-Y) - sin(th)(x0-X) |
Strip = surf(abs(cos(th)*(surf(:,y)-y0) - sin(th)*(surf(:,x)-x0)) < w, :);

% curve fitting using loess
SmoothCurve = fit( Strip(:,[x,y]), Strip(:,y), 'lowess' );

npts  = max( vecnorm(Strip(:,[x,y])-[x0,y0],2,2) )/dt;
xRan  = x0 + cos(th)*dt*(0:npts);
yRan  = y0 + sin(th)*dt*(0:npts);
zRan  = feval( SmoothCurve, [xRan', yRan'] );

Curve = [xRan', yRan', zRan];
end

function [PtsLocs] = PointsInCurve( curve, d0, dt, Npts )
% Locate points with given point-to-point distance within a curve.
% 
%-------------------------------------------------------------------------
% INPUT
%     curve  Smoothened curve interpolating the strip obtained from the
%            surface, 1x(#points)
%        d0  Distance from beginning of curve to the first point, in mm.
%        dt  Point-to-point distance, in mm.
%      Npts  Number of points to locate, integer greater than 1.
%
%-------------------------------------------------------------------------
% OUTPUT
%    PtsLoc  Coordinates of points, (Npts)x 3
%
%-------------------------------------------------------------------------
% Author: Julio Cesar Enciso-Alva, 2023
%         juliocesar.encisoalva@uta.edu
%

PtsLocs = zeros( Npts, 3 );

% cumulative length of curve
CumLen = zeros(size( curve ));
for ii = 2:length(CumLen)
  CumLen(ii) = CumLen(ii-1) + norm(curve(ii,:) - curve(ii-1,:));
end

% iteration to get points
cursor = d0;
for ii = 1:Npts
  [~,idx] = min(abs( CumLen - cursor ));
  PtsLocs(ii,:) = curve(idx,:);
  cursor = cursor + dt;
end

end