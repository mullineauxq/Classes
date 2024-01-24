function d = ddx_bwd(f,dx)
% ddx_bwd takes an ngrid matrix f and x-spacing of points 'dx' and computes the matrix of first order backwards 1st differences 
% ddx_bwd takes an n x m input matrix f and returns an n x m matrix d which
% contains the numerical first order 1st derivative of f across the
% columns. At the leftmost boundaries (n,:) a forward difference is taken
% instead

% matlab works in meshgrid but we use ngrid, therefore n should be #
% columns and m should be # of rows, we'll rename n to be nx and m to be ny
% for convenience

%/////////////////////////////Begin Function//////////////////////////////////////////


[nx,ny] = size(f);
d=zeros(size(f));


% loop through all values except right boundary, take the equation for 
% backward diff at that point

for j = 1:ny
    for i = 2:(nx)
        d(i,j) = (f(i,j)-f(i-1,j))/dx ;% converted to ngrid from meshgrid then performed forward diff
    end
end


for  k = 1:ny 
    d(1,k) = (f(2,k)-f(1,k))/dx; 

end
end