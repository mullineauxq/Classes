numer = [.1];
denom = [1 , .1, 0];
d = 6;
h = 2

gs = tf(numer,denom);
gz = c2d(gs,h,['zoh'])









