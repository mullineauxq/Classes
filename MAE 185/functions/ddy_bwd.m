function dfdy = ddy_bwd(f,dy,periodicity)

arguments 
    f
    dy
    periodicity = 'n'
end
    % determine field size
    [nx,ny]     = size(f);

    % allocate return field
    dfdy        = zeros(nx,ny);
    
    % backward difference
    for i=1:nx
        for j=2:ny
            dfdy(i,j) = (f(i,j)-f(i,j-1))/dy;
        end
    end

    if periodicity =='p'
        for  k = 1:nx 
            dfdy(k,1) = (f(k,1)-f(k,ny))/dy;
        end
    else
        for  k = 1:ny 
            dfdy(k,1) = (f(k,2)-f(k,1))/dy; 
        end
    end
end