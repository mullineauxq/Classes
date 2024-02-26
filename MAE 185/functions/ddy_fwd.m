function dfdy = ddy_fwd(f,dy,periodicity)
    arguments 
        f
         dy
         periodicity = 'n'
    end
    % determine field size
    [nx,ny]     = size(f);
    
    % allocate return field
    dfdy        = zeros(size(f));
    
    % forward difference
    for i=1:nx
        for j=1:ny-1
            dfdy(i,j) = (f(i,j+1)-f(i,j))/dy;
        end
    end
    
    if periodicity == 'p'
  
        for i=1:nx
            dfdy(i,ny) = (f(i,1)-f(i,ny))/dy;
        end
    else
        % backward difference for last point
    j = ny;
        for i=1:nx
            dfdy(i,j) = (f(i,j)-f(i,j-1))/dy;
        end
    end

    
end