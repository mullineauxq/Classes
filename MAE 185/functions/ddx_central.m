function d = ddx_central(f,dx,periodic)
% ddx_central takes an ngrid matrix f and x-spacing of points 'dx' and computes the matrix of second order central 1st differences 
% ddx_central takes an n x m input ngrid matrix f and returns an n x m ngrid matrix d which
% contains the numerical first order 1st derivative of f across the
% columns. At the leftmost boundaries (n,:) a forward difference is taken
% instead, and at the rightmost boundary, a backwards difference
% if a periodic boundary is used, pass char 'p' in the periodic field, and
% a central periodic difference will be taken at both boundaries

% matlab works in meshgrid but we use ngrid, therefore n should be #
% columns and m should be # of rows, we'll rename n to be nx and m to be ny
% for convenience

%/////////////////////////////Begin Function//////////////////////////////////////////

arguments 
    f
    dx
    periodic = 'n';
end
[nx,ny] = size(f);
d=zeros(size(f));

% loop through all values except right boundary, take the equation for 
% central diff at that point

for j = 1:ny
    for i = 2:(nx-1)
        d(i,j) = (f(i+1,j)-f(i-1,j))/(2*dx) ;% converted to ngrid from meshgrid then performed forward diff
    end
end

% take the 2nd order forward 1st diff at the right boundary  (X=1)
for  k = 1:ny 
    d(1,k) = (-3*f(1,k)+4*f(2,k)-f(3,k))/(2*dx); 

end

% take the 2nd order backward 1st diff at the left boundary  (X=nx)
for  k = 1:ny 
    d(nx,k) = (3*f(nx,k)-4*f(nx-1,k)+f(nx-2,k))/(2*dx); 

end

if periodic == 'p' % user indicates periodic boundaries

    %Loop through all y values at left and right boundary
    for j = 1:ny
    d(1,j) = (f(2,j)-f(nx,j))/(2*dx);
    d(nx,j) = (f(1,j)-f(nx-1,j))/(2*dx);
    end
end
end