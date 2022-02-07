function dpsit = qgrhs(psi, L, Dx, Dy, ymat, Re, Ro)
% Here psi represents the stream function, and q represents the vorticity
% They are related by the poisson equation
%  (1)  -L*psi = q
% The QG equations are
%  (2)    dq/dt + J(psi, q) = 1/Ro * Dx*psi + 1/Re * L*q + 1/Ro * F
% Here J is the Jacobian term (a quadrativ term)
% Ro is the Rossby number
% Re is the Reynolds number
% L is the Laplacian
% Dx is the derivative with respect to x,
% and F is the forcing term (here double gyre)
%
% We first solve for dq/dt by use of (2)
% then solve for dpsi/dt by solving the Poisson equation (1)
% Hope this helps!


% Calculate the vorticity
q = -L*psi;

% calculate Arakawa approximation to the Jacobian term
% see https://www.sciencedirect.com/science/article/pii/0021999166900155
dpsix = Dx*psi;
dpsiy = Dy*psi;
dqx = Dx*q;
dqy = Dy*q;
J1 = dpsix.*dqy     - dpsiy.*dqx;
J2 = Dx*(psi.*dqy) - Dy*(psi.*dqx);
J3 = Dy*(q.*dpsix) - Dx*(q.*dpsiy);

% This is the full Jacobian term
J = -(J1 + J2 + J3)/3;

% double gyre forcing term
F = sin(pi*(ymat - 1));

% this is the derivative of the vorticity with respect to time
dqt = -J + (1/Ro)*(dpsix) + (1/Re)*(L*q) + (1/Ro)*F;

% solve the poisson equation for the stream function
dpsit = -L\dqt;
end
