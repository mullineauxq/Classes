function pixl(image_filename)
%% Written by QBM this code allows you to simply open an image in a pixel mapped format

myfile=image_filename;
myimg=imread(myfile);
h = imshow(myimg);
impixelinfo(h);
end