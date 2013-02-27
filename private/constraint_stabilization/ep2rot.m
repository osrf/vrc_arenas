function [varargout] = ep2rot(ep)
% ep2rot   
%-------
% Given the Euler parameter vector, this function returns the rotation,
% matrix, R, and optionally the velocity kinematics matrix, B, defined
% below.
%
% Syntax:       ep2rot ([ep0 ep1 ep2 ep3])       returns R
%    or     R = ep2rot ([ep0 ep1 ep2 ep3])       returns R
%    or [R,B] = ep2rot ([ep0 ep1 ep2 ep3])       returns R and B
%
%  ******************** ARGUMENT DEFINITIONS ***********************
%
%   ep    Euler parameter vector of length four; ep = [ep0 ep1 ep2 ep3].
%   R     The rotation matrix correpsonding to the Euler paramters.
%   B     The velocity kinematic map for the Euler parameters.
%
% *********************** END ARGUMENT DEFINITIONS **********************
%
% Euler Parameter Definition
% --------------------------
% The parameters define an axis direction, e, and an angle of rotation, Phi,
% about that axis.  Orientation of a frame w.r.t. a base frame is obtained
% by rotating the base frame about e by the angle Phi.
% For mathematical formulations, see Junkins, Schaub, "Analytical Mechanics
% of Aerospace Systems," chapter 3, section 4, pp 75-80.
%
%   The Euler parameters are defined by (in Junkins/Schaub notation):
%        ep0 = cos(Phi/2)           = epsilon4 (Craig's notation)
%	     ep1 = e1*sin(Phi/2)        = epsilon1 (        "       )
%	     ep2 = e2*sin(Phi/2)        = epsilon2 (        "       )
%	     ep3 = e3*sin(Phi/2)        = epsilon3 (        "       )
%   where e = [ep1 ep2 ep3]^T.
%
%   The rotation matrix R describing the body frame orientation w.r.t. the 
%   base frame orientation is given by:
% R = [  1-2*(ep2^2 + ep3^2)   2*(ep1*ep2 - ep3*ep0)   2*(ep1*ep3 + ep2*ep0);
%       2*(ep1*ep2 + ep3*ep0)   1-2*(ep1^2 + ep3^2)    2*(ep2*ep3 - ep1*ep0);
%       2*(ep1*ep3 - ep2*ep0)  2*(ep2*ep3 + ep1*ep0)    1-2*(ep1^2 + ep2^2)  ];
%
%   The matrix mapping the angular velocity of the body (in body coordinates)
%   to the time derivative of the Euler parameter vector is:
%        B = [-ep1 -ep2 -ep3;
%              ep0 -ep3  ep2;
%              ep3  ep0 -ep1;
%             -ep2  ep1  ep0] / 2;
%   where d(ep)/dt = B*omega_B.

% written by: Jeff Trinkle 10/22/99
%     Added option to return the velocity kinematic map of the Euler
%        parameter vector.  12/12/99.
%---------------------------------------------------------------------------

%------------------------- Begin Input Checking ----------------------------
%
% Check number of arguments passed in
%
tolerance = eps * 1e9;
if (nargin == 1)
  elength = length(ep);
  if elength == 4,
    ep0 = ep(1);
    ep1 = ep(2);
    ep2 = ep(3);
    ep3 = ep(4);
  else
    disp (' ');
    disp ('****************** Error in ep2rot()!!!! *******************');
    disp ('*****************************************************************');
    disp ('ep2tr(ep) expects an Euler parameter vector of length 4, but');
    string1 = ['it has length ' num2str(elength) '!!!!!!!'];
    disp (string1);
    disp ('*****************************************************************');
    return;
  end
else
  disp (' ');
  disp ('****************** Error in ep2rot()!!!! *******************');
  disp ('*****************************************************************');
  disp ('ep2tr(ep) expects 1 input argument, but you have passed in');
  string1 = [num2str(nargin) ' arguments!!!!!!!'];
  disp (string1);
  disp ('*****************************************************************');
  return;
end

ep_norm = norm(ep);
if (abs(ep_norm - 1) > tolerance)
  disp (' ');
  disp ('****************** Warning in ep2rot()!!!! *******************');
  disp ('*****************************************************************');
  disp ('The Euler parameter vector passed to ep2rot() should have Euclidean norm of 1, but');
  string1 = ['its norm is ' num2str(ep_norm) '.  Normalizing and continuing!'];
  disp (string1);
  disp ('*****************************************************************');
  ep0 = ep0/ep_norm;
  ep1 = ep1/ep_norm;
  ep2 = ep2/ep_norm;
  ep3 = ep3/ep_norm;
end

%
%-------------------------- End input checking -----------------------------------

% Define the rotation matrix
R = [ 1-2*(ep2^2 + ep3^2)   2*(ep1*ep2 - ep3*ep0)   2*(ep1*ep3 + ep2*ep0);
     2*(ep1*ep2 + ep3*ep0)   1-2*(ep1^2 + ep3^2)    2*(ep2*ep3 - ep1*ep0);
     2*(ep1*ep3 - ep2*ep0)  2*(ep2*ep3 + ep1*ep0)    1-2*(ep1^2 + ep2^2)  ];
varargout{1} = R;

if nargout == 2
   % The transform that relates the angular velocity of the body (expressed in
   % the body frame) to the time derivative of the Euler parameters
   B = [-ep1 -ep2 -ep3;
         ep0 -ep3  ep2;
         ep3  ep0 -ep1;
        -ep2  ep1  ep0] / 2; 
   varargout{2} = B;
end

if nargout > 2
   fprintf('\n\n ******************************************************');
   fprintf('\n Error using ep2rt()!!!!!!  It has been called with %d', nargout);
   fprintf('\n return arguments, but there must be two or fewer.');
   fprintf('\n ******************************************************');
end
