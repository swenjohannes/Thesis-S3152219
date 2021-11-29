function varargout = sfunc(varargin)
   [varargout{1:nargout}] = feval(varargin{:});
end
function [] = gravitationalForce(m1,m2,d)
velocityCalc(m1, m2, d);
end
function [] = velocityCalc(u,a,s)
u + a + s;
end
function [] = distanceTraveled(u,a,t)
...
end