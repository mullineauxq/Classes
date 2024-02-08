function d = ddx_bwd(f,dx,periodicity)
% ddx_bwd takes an ngrid matrix f and x-spacing of points 'dx' and computes the matrix of first order backwards 1st differences 
% ddx_bwd takes an n x m input matrix f and returns an n x m matrix d which
% contains the numerical first order 1st derivative of f across the
% columns. At the leftmost boundaries (n,:) a forward difference is taken
% instead. if periodic boundaries are used, pass the char 'p' into the
% periodicity field and at the left bondary, a periodic condition will be
% applied

% matlab works in meshgrid but we use ngrid, therefore n should be #
% columns and m should be # of rows, we'll rename n to be nx and m to be ny
% for convenience

%/////////////////////////////Begin Function//////////////////////////////////////////
arguments
    f
    dx
    periodicity = 'n'
end

[nx,ny] = size(f);
d=zeros(size(f));


% loop through all values except right boundary, take the equation for 
% backward diff at that point

for j = 1:ny
    for i = 2:(nx)
        d(i,j) = (f(i,j)-f(i-1,j))/dx ;% converted to ngrid from meshgrid then performed forward diff
    end
end

% apply a specialized approach to the left boundary, taking a forward diff
% if periodicity is not allowed, and take a periodic backwards diff if a
% periodic condition is passed by the user

if periodicity =='p'
    for  k = 1:ny 
        d(1,k) = (f(1,k)-f(nx,k))/dx; 
    
    end
else
    for  k = 1:ny 
        d(1,k) = (f(2,k)-f(1,k))/dx; 
    
    end
end

end