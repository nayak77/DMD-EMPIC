function [w1,w2,w3] = bary_coord(x1,y1,x2,y2,x3,y3,x,y)
% Input: Triangle vertices coordinates (x1,y1), (x2,y2), (x3,y3) and
% coordinate of the desired point (x,y).
% Output: barucentric coordinates or weights for (x,y)

w1 = ((y2-y3)*(x-x3) + (x3-x2)*(y-y3))/((y2-y3)*(x1-x3) + (x3-x2)*(y1-y3));
w2 = ((y3-y1)*(x-x3) + (x1-x3)*(y-y3))/((y2-y3)*(x1-x3) + (x3-x2)*(y1-y3));
w3 = 1 - w1 - w2;
end