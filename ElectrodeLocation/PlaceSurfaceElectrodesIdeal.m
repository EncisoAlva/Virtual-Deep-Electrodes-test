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