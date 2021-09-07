function [] = drawCircle(r, centreX, centreY, varargin)
x = 0:0.01:2*pi;
xx = r * cos(x);
yy = r * sin(x);

if nargin >= 4
    plot(centreX + xx, centreY + yy, varargin{:});
else
    plot(centreX + xx, centreY + yy);
end

end
