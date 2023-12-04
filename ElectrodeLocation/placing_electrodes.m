function varargout = process_PlaceElectrodes( varargin )
%
end

%% INTRO

% Electrode positions in some animal models --such as minipig, Sus scrofa-- 
% are not yet fully standardized. Thus, the electrode positions must be 
% determined manually followingsome protocol.

%% SURFACE ELECTRODES
function locs = PlaceSurfaceElectrodesIdeal()
% This script determine locations of an NxM surface electordes following 
% this simple protocol:
%  1. Identify the Central Line: intersection of the plane between brain 
%     hempispheres and the envelope of brain cortex.
%  2. Identify Anterior Edge: the most anterior points of the brain, 
%     approximately ortogonal to both the Central Line and the 
%     Posterior-Anterior line.
%  3. Set the first electrode at 
%       P cm from Anterior Edge on the Anterior-Posterior direction
%       R cm from Central Line on the Left-Right direction
%  4. Set the other electrodes following a rectangular grid with one side 
%     at a slope angle of T degrees with respect to the Central Line.
%     The grid has center-to-center distances of dP cm in the 
%     Anterior-Posterior direction, and dR cm in the Left-Right direction. 
%     The grid is assumed to follow the curvature of the cortex envelope.
% 
% This function makes assumptions of ideal placing of electrodes. Other
% functions are provided for further adjustments based on observations.
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





%addpath 'C:\brainstorm3'
%brainstorm



%%
% PARAMETERS

%cortex_conn = example1.Faces;
%cortex_vert = example1.Vertices;

%cortex_conn = cortex.Faces;
%cortex_vert = cortex.Vertices;

dev = true;
vis_width = 0.001;
vis_tol = vis_width/20;
max_iter = 100;

% strip of AxB electrodes is --for now-- assumed parallel to medial line
% AxB means each B electrode line is perpendicular
ElecGridNum = [1,4];
ElecGridSep = [1,0.01];
ElectRefStart = [0.01,0.01];
ElectSide = 'left';

if(dev)
    % DELETE
    figure()
    trisurf(cortex_conn, cortex_vert(:,1),cortex_vert(:,2),cortex_vert(:,3),...
        'FaceColor',[.85 .85 .85])
    title('Only Cortex Surface')
    xlabel('x')
    ylabel('y')
    zlabel('z')
end

%%
% IDENTIFICATION OF ANTERIOR EDGE

% finding anterior terminal
% 1. find max
% 2. select pts close to the max
% 3. find quasi-max (~100th percentile)
% 4. repeat 2,3 until convergence
[M,~] = max( cortex_vert(:,1) );
old_M = M+vis_tol*2;
iter = 0;
while ( abs(M-old_M) >= vis_tol ) && ( iter < max_iter )
    strip_idx = abs( cortex_vert(:,1)-M ) < vis_width;
    old_M = M;
    M = quantile(cortex_vert(strip_idx,1), 0.98);
    iter = iter +1;
end
if dev
    % DELETE LATER
    figure()
    scatter(cortex_vert(:,1),cortex_vert(:,3),'filled')
    xlabel('x')
    ylabel('z')
    title('Only cortex surface (mesh vertices)')
    subtitle('Sagittal view')
    hold on 
    xline(M,'-','Anterior border','LineWidth',1)
end

% finding the anterior edge
% 1. find most superior point in anterior terminal
% 2. use the xz coordinates as anterior edge
strip_idx = abs( cortex_vert(:,1)-M ) < vis_width;
tmp_vec = cortex_vert(strip_idx,:);
m = quantile( tmp_vec(:,3)-abs(tmp_vec(:,1)-M), 0.98 );
[~,m_idx] = min( abs(tmp_vec(:,1)-M)+abs(tmp_vec(:,3)-m) );
AntEdge_x = tmp_vec(m_idx,1);
AntEdge_z = tmp_vec(m_idx,3);
if dev
    % Delete
    figure()
    scatter(cortex_vert(:,1),cortex_vert(:,3),'filled')
    xlabel('x')
    ylabel('z')
    title('Cortex Surface (Mesh Vertices)')
    subtitle('Sagittal View')
    hold on 
    scatter(AntEdge_x,AntEdge_z,200,'k','filled')
    xline(M,'--','LineWidth',1)
    legend('','Anterior Edge','Anterior Extreme')
end
if dev
    vis_idx = vecnorm( cortex_vert(:,[1,3])-[AntEdge_x,AntEdge_z], 2, 2 ) < 2*vis_width;
    figure()
    trisurf(cortex_conn, cortex_vert(:,1),cortex_vert(:,2),cortex_vert(:,3),...
        'FaceColor',[.85 .85 .85])
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('Cortex Surface')
    hold on 
    scatter3(cortex_vert(vis_idx,1),cortex_vert(vis_idx,2),cortex_vert(vis_idx,3),'r','filled')
    legend('','Anterior Edge')
end

%%

% draw the AxB grid electrodes algorithmicalli
% 1. draw medial line, for reference
% 2. mark B regular dots in medial line
% 3. draw lines parallel to medial line
% 4. mark A regular dots on each one
% 5. TODO correction for robustness

% 1.
med_line = cortex_vert( abs(cortex_vert(:,2))<vis_width, :);
med_line = med_line( med_line(:,3)>AntEdge_z, :); % only upper part is relevant
if dev
    % DELETE
    figure()
    trisurf(cortex_conn, cortex_vert(:,1),cortex_vert(:,2),cortex_vert(:,3),...
        'FaceColor',[.85 .85 .85])
    xlabel('x')
    ylabel('y')
    zlabel('y')
    hold on
    scatter3(med_line(:,1),med_line(:,2),med_line(:,3),'filled')
    scatter3(cortex_vert(vis_idx,1),cortex_vert(vis_idx,2),cortex_vert(vis_idx,3),'r','filled')
end

x_0 = min( med_line(:,1) );
x_f = max( med_line(:,1) );
x_lo = x_0:0.0001:x_f; 
y_lo = x_lo*0 + 0;
model = fit( med_line(:,[1,2]), med_line(:,3), 'lowess' );
z_up = feval( model, [x_lo',y_lo'] );

MLcurve = [x_lo', y_lo', z_up]; % medial line, stilyzed curve
guides_idx = MLcurve;

% 2.
lin_int = zeros( 1, size(x_lo,2) );
for ii = 2:size(x_lo,2)
    lin_int(ii) = lin_int(ii-1) + norm( MLcurve(ii,:)-MLcurve(ii-1,:) );
end
lin_int = lin_int(end) - lin_int;

ctr_mark = zeros( ElecGridNum(2) ,3 );
cursor = ElectRefStart(2);
for iter = 1:ElecGridNum(2)
    [~,idx] = min(abs( lin_int-cursor ));
    ctr_mark(iter,:) = MLcurve(idx,:);
    cursor = cursor + ElecGridSep(2);
end

%%%

ElectrodeCoords = zeros( ElecGridNum(1)*ElecGridNum(2), 3 );
ElectrodePositions = zeros( ElecGridNum(1)*ElecGridNum(2), 2 );
lazy_counter = 1;
pos = [1,1];

for iter2 = 1:ElecGridNum(2)
    % 3.
    aux_line = cortex_vert( abs(cortex_vert(:,1)-ctr_mark(iter2,1))<vis_width, :);
    aux_line = aux_line( aux_line(:,3)>AntEdge_z, :); % only upper part is relevant
    
    if strcmp( ElectSide, 'left' )
        y_0 = ctr_mark(iter2,2);
        y_f = max( aux_line(:,2) );
    else
        y_0 = min( aux_line(:,2) );
        y_f = ctr_mark(iter2,2);
    end
    y_lo = y_0:0.0001:y_f; 
    x_lo = y_lo*0 + ctr_mark(iter2,1);
    model = fit( aux_line(:,[1,2]), aux_line(:,3), 'lowess' );
    z_up = feval( model, [x_lo',y_lo'] );

    MLcurve = [x_lo', y_lo', z_up]; % medial line, stilyzed curve
    guides_idx = [guides_idx; MLcurve];

    % 4.
    lin_int = zeros( 1, size(x_lo,2) );
    for ii = 2:size(x_lo,2)
        lin_int(ii) = lin_int(ii-1) + norm( MLcurve(ii,:)-MLcurve(ii-1,:) );
    end
    %lin_int = lin_int(end) - lin_int;

    cursor = ElectRefStart(1);
    pos(2) = 1;
    for iter = 1:ElecGridNum(1)
        [~,idx] = min(abs( lin_int-cursor ));
        ElectrodeCoords(lazy_counter,:) = MLcurve(idx,:);
        ElectrodePositions(lazy_counter,:) = pos;
        cursor = cursor + ElecGridSep(1);
        lazy_counter = lazy_counter+1;
        pos(2) = pos(2)+1;
    end
    pos(1) = pos(1)+1;
end

if(dev)
    figure()
    trisurf(cortex_conn, cortex_vert(:,1),cortex_vert(:,2),cortex_vert(:,3),...
        'FaceColor',[.85 .85 .85],'FaceAlpha',.4)
    xlabel('x')
    ylabel('y')
    zlabel('y')
    title('Cortex surface')
    hold on
    scatter3(guides_idx(:,1),guides_idx(:,2),guides_idx(:,3),'filled')
    scatter3(cortex_vert(vis_idx,1),cortex_vert(vis_idx,2),cortex_vert(vis_idx,3),'r','filled')
    scatter3(ElectrodeCoords(:,1),ElectrodeCoords(:,2),ElectrodeCoords(:,3),200,'filled')
    legend('','Reference lines','Anterior Edge','Electrodes')
end

%%
% labels

writematrix(ElectrodeCoords)