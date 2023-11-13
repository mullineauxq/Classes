function R = rotate(axis,angle)
%% R = rotate(angle,axis) takes in an angle in degrees and outputs a 3x3 rotation matrix about the axis inputted as a char x y or z by the user

syms t
Rx = [1 0 0; 0 cosd(t) sind(t); 0 -sind(t) cosd(t)];
Ry = [cosd(t) 0 -sind(t); 0 1 0; sind(t) 0 cosd(t)];
Rz = [cosd(t) sind(t) 0; -sind(t) cosd(t) 0; 0 0 1];

switch axis
    case 'x'
    R = subs(Rx,t,angle);
    case 'y'
    R = subs(Ry,t,angle);
    case 'z'
    R = subs(Rz,t,angle);
end

end