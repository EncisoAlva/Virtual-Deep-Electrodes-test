function [ElectrodeLocations] = PlaceInsertedElectrodesIdeal( ...
  surface, nPA, nLR, PA0, LR0, dPA, dLR, TH, labels, ...
  params, debug )
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

% rename of parameters
coord_og = surface.Vertices;
range_old  = range(coord_og,1);
center_old = min(coord_og,[],1) + range_old/2;
coord_og = coord_og - center_old;
coord = coord_og;

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
while (iter < params.max_iter) && flag
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
  scatter(ptPA(1),ptPA(3),200,'k','filled')
  scatter(ptIS(1),ptIS(3),200,'k','filled')
  legend({'','Anterior Edge','',''}, 'Location','southwest')
end

coord = coord( coord(:,dir.IS) >= ptEDGE(dir.IS) ,:); % only upper part is relevant

%% LOCATION OF ELECTRODES
% 1. draw medial line, for reference
% 2. mark B regular dots in medial line
% 3. draw lines parallel to medial line
% 4. mark A regular dots on each one
% 5. TODO correction for robustness

for lazy_counter = 1:length(ElectSide)
CurrentSide  = ElectSide{lazy_counter};
CurrentLabel = labels{lazy_counter};
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
CtrlStrip = coord( abs(coord(:,dir.LR)) < params.vis_width, : );
CtrlLine  = SmoothCurveInterpolation( coord, ...
  min( CtrlStrip(:,dir.PA) ), 0, 0, params.vis_width, 0.0001, debug.figs_detail );

tmpPt = PointsInCurve( CtrlLine, PA0, dPA, 1 );
auxLR = SmoothCurveInterpolation( coord, ...
  tmpPt(dir.PA), tmpPt(dir.LR), sg*90, params.vis_width, 0.0001, debug.figs_detail );

FirstPt = PointsInCurve( auxLR, LR0, dLR, 1 );
ElecLocs(1,1,:) = FirstPt;

% iterate over Posterior-Anterio, then over Left-Right
auxPA = SmoothCurveInterpolation( coord, ...
  FirstPt(dir.PA), FirstPt(dir.LR), TH, params.vis_width, 0.0001, debug.figs_detail );
ElecLocs(1,1:nPA,:) = PointsInCurve( auxPA, 0, dPA, nPA );

for ii = 1:nPA
  refPt = ElecLocs(1,ii,:);
  auxLR = SmoothCurveInterpolation( coord, ...
    refPt(dir.PA), refPt(dir.LR), 180-TH+sg*90, params.vis_width, 0.0001, debug.figs_detail );
  ElecLocs(1:nLR,ii,:) = PointsInCurve( auxLR, 0, dLR, nLR );
end

if(debug.figs)
    figure()
    trisurf(surface.Faces, coord_og(:,1),coord_og(:,2),coord_og(:,3),...
        'FaceColor',[.85 .85 .85],'FaceAlpha',.4)
    xlabel('Posterior-Anterior')
    ylabel('Left-Right')
    zlabel('Inferior-Superior')
    title('Cortex surface')
    hold on
    for ii = 1:nPA
      for jj = 1:nLR
        scatter3(ElecLocs(jj,ii,1),ElecLocs(jj,ii,2),ElecLocs(jj,ii,3),...
          200,'r','filled')
      end
    end
end


ndigits   = ceil(log10(nPA*nLR));
formatStr = strcat("%",num2str(ndigits),"d");

ElecTable_pre  = zeros(nPA*nLR,3); 
ElecLabels     = cell(nPA*nLR,1);

counter   = 1;
for lr = 1:nLR
for pa = 1:nPA
  ElecTable_pre(counter,:) = reshape(ElecLocs(lr,pa,:),[1,3]) ...
    + center_old; % undo the re-center
  ElecLabels{counter}  = strcat(CurrentLabel,num2str(counter,formatStr));
  counter = counter+1;
end
end

ElecTable = array2table(ElecTable_pre, 'VariableNames',{'PA','LR','IS'});
ElecTable.label = ElecLabels;

if lazy_counter==1
  ElectrodeLocations = ElecTable;
else
  ElectrodeLocations = [ElectrodeLocations; ElecTable];
end

end

writetable(ElectrodeLocations)
end