%% This script is to correlate displacement of the jaws to data collected by the instron tensile tester 
% it imports instron data, then interpolates the force at all displacements
% passed to 


%% CHANGE THESE THESE LINES TO THE NAME OF YOUR INSTRON DATA FILE and photos array
instron_data = readtable('PMMA_dogbone_20240131_134446_1.csv');
photos_data = readtable('Samples Measurments - Custom Design.csv');

%% Now, if you analyzing photos from sample 2, replace this 6 with a 4, if sample 1, use 2
measured_displacements = table2array(photos_data((2:end-1),4));

% should be no need for alterations after this point
photo_count = length(measured_displacements);
estimated_load = zeros(photo_count,1);



time  = squeeze(table2array(instron_data((3:end),2)));
displacement_in = squeeze(table2array(instron_data((3:end),3)));
force_lbf = squeeze(table2array(instron_data((3:end),4)));

for i = 1:photo_count
    estimated_load(i,1) = interp1(displacement_in,force_lbf,measured_displacements(i));
end

force_estimates = array2table(estimated_load);
