% In my spatial discretization, I do not explicitly account for the boundary points
% as they are assumed to be zero, therefore the number of points in the x and y domains
% only account for the INTERIOR points.

% defining the domain the in the x and y directions
xdomain = [0, 1];
ydomain = [0, 2];

% number of points in the interior of the domain in the x direction
nx = 63;

% number of points in the interior of the domain in the y direction
ny = 127;

% calculate the distance between the interior points for the finite difference
% we add one, as there are nx + 1 spaces between ALL the points (including the boundaries)
hx = (xdomain(2) - xdomain(1))/(nx + 1);
hy = (ydomain(2) - ydomain(1))/(ny + 1);

% Calculate the matrix of y value where the forcing term will be evaluated.
% In the rhs function, we
ys = linspace(ydomain(1), ydomain(end), ny + 2);
ys = ys(2:end-1);
ymat = repmat(ys.', 1, nx);
ymat = reshape(ymat.',  nx*ny, 1);

% The full Laplacian can be written as a sum of kronecker products in the form
% L = I_{ny} (x) Lx + Ly (x) I_{nx}
% where Lx is the 1 dimensional laplacian in the x direction
% and Ly is the 1 dimensional laplacian in the y direction

% The Laplacian in the x direction using a second order finite difference
Lx = spdiags(repmat([1 -2 1], nx, 1), [-1 0, 1], nx, nx)/(hx^2);
% The Laplacian in the y direction using a second order finite difference
Ly = spdiags(repmat([1 -2 1], ny, 1), [-1 0, 1], ny, ny)/(hy^2);
L = kron(speye(ny), Lx) + kron(Ly, speye(nx));

% We do a similar thing with the derivative operators in the x and y directions
Dx1D = spdiags(repmat([-1 0 1], nx, 1), [-1 0, 1], nx, nx)/(2*hx);
Dx = kron(speye(ny), Dx1D);
Dy1D = spdiags(repmat([-1 0 1], ny, 1), [-1 0, 1], ny, ny)/(2*hy);
Dy = kron(Dy1D, speye(nx));

Re = 450;
Ro = 0.0036;

% We are now really do define the right hand side function of the ordinary differential equation
f = @(~, psi) qgrhs(psi, L, Dx, Dy, ymat, Re, Ro);

% We set the timespan to 50 days
oneday = 0.0109;
timespan = [0, 50*oneday];

% The initial condition needs to be set manually. This is a temporary one.
y0 = zeros(nx*ny, 1);

% Now we can solve the ordinary differential equation in time with ode45
sol = ode45(f, timespan, y0);

% We can plot the solution after one day
% Make sure to pay attentio to the transpose!
imagesc(reshape(sol.y(:, end), nx, ny).');



