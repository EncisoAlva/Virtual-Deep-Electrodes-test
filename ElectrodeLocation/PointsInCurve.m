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
CumLen = zeros(size( curve,1 ),1);
for ii = 2:length(CumLen)
  CumLen(ii) = CumLen(ii-1) + norm(curve(ii,:) - curve(ii-1,:));
end

% iteration to get points
for ii = 1:Npts
  [~,idx] = min(abs( CumLen - ( d0+(dt)*(ii-1) ) ));
  PtsLocs(ii,:) = curve(idx,:);
end

end