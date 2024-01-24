function d = d2dy2(f,dy)
% d2dy takes an ngrid matrix f and y-spacing of points 'dy' and computes the matrix of second order central 2nd differences 
% ddx_central takes an n x m input ngrid matrix f and returns an n x m ngrid matrix d which
% contains the numerical second order 2nd derivative of f across the
% columns. At the bottom boundary (n,:) a 2nd order 2nd forward difference is taken
% instead, and at the top boundary, a backwards difference

% matlab works in meshgrid but we use ngrid, therefore n should be #
% columns and m should be # of rows, we'll rename n to be nx and m to be ny
% for convenience

%/////////////////////////////Begin Function//////////////////////////////////////////


[nx,ny] = size(f);
d=zeros(size(f));

% loop through all values except right boundary, take the equation for 
% central diff at that point

for i = 1:nx
    for j = 2:(ny-1)
        d(i,j) = (f(i,j+1)-2*f(i,j)+f(i,j-1))/(dy)^2 ;% converted to ngrid from meshgrid then performed forward diff
    end
end

% take the 2nd order forward 2nd diff at the bottom
% and 2nds backwards diff at the top
for  k = 1:nx 
    d(k,1) = (2*f(k,1)-5*f(k,2)+4*f(k,3)-f(k,4))/dy^2;
    d(k,ny) = (2*f(k,ny)-5*f(k,ny-1)+4*f(k,ny-2)-f(k,ny-3))/dy^2;
end


end