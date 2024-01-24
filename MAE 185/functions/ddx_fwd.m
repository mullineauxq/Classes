function d = ddx_fwd(f,dx)
% ddx_fwd takes an ngrid matrix f and x-spacing of points 'dx' and computes the matrix of first order forward 1st differences 
% ddx_fwd takes an n x m input matrix f and returns an n x m matrix d which
% contains the numerical first order 1st derivative of f across the
% columns. At the rightmost boundaries (n,:) a forward difference is taken
% instead

% matlab works in meshgrid but we use ngrid, therefore n should be #
% columns and m should be # of rows, we'll rename n to be nx and m to be ny
% for convenience


%/////////////////////////////Begin Function//////////////////////////////////////////

[nx,ny] = size(f);
d=zeros(size(f));
% loop through all values except right boundary, take the equation for 
% forward difference at that point

for j = 1:ny
    for i = 1:(nx-1)
        d(i,j) = (f(i+1,j)-f(i,j))/dx ;% converted to ngrid from meshgrid then performed forward diff
    end
end

%finally loop through the rightmost column and perform first order backwards 1st difference 
for  k = 1:ny 
    d(nx,k) = (f(nx,k)-f(nx-1,k))/dx;
end

end