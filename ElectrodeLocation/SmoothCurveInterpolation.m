function [Curve] = SmoothCurveInterpolation( ...
  surf, x0,y0, theta, w, dt, imshow )
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
th = (theta)*pi/180;

% lazy variables
x=1; y=2; z=3;

% dist from (X,Y) to line passing trough (x0,y0) at angle th
%   | cos(th)(y0-Y) - sin(th)(x0-X) |
Strip = surf(abs(cos(th)*(surf(:,y)-y0) - sin(th)*(surf(:,x)-x0)) < w, :);

% curve fitting using loess
SmoothCurve = fit( Strip(:,[x,y]), Strip(:,z), 'lowess' );

npts  = max( vecnorm(Strip(:,[x,y])-[x0,y0],2,2) )/dt;
xRan  = x0 + cos(th)*dt*(0:npts);
yRan  = y0 + sin(th)*dt*(0:npts);
zRan  = feval( SmoothCurve, [xRan', yRan'] );

Curve = [xRan', yRan', zRan];

if(imshow)
    figure()
    scatter3(surf(:,1),surf(:,2),surf(:,3),'filled')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    hold on
    scatter3(Strip(:,1),Strip(:,2),Strip(:,3),100,'k','filled')
    scatter3(Curve(:,1),Curve(:,2),Curve(:,3),100,'g','filled')
end
end