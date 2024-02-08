% DICe result reader, N. Boechler (1/2023)
% put this in a folder with your DICe_solution_0*.txt files

clc;
clear all;
close all;

%declare the DICE solution numbers that this trial starts and ends at
start_int = 4472;
final_int  = 4477;

number_of_photos = final_int-start_int+1;

% allocates a vector to store strain at each photo (in the gage)
strain_values = zeros(number_of_photos,1); 


myfiles=dir('*.txt');
myfolder=pwd;

for m = 1:number_of_photos

mydata=readmatrix(myfiles(m).name);

strain_yy=mydata(:,12);
xpos=mydata(:,2);
ypos=mydata(:,3);

figure(01)
scatter(xpos,ypos,[],strain_yy,'filled')
colorbar
xlabel('x')
ylabel('y')
set(gca,'fontsize',20,'linewidth',2)

initial_x  = xpos(1);


%% Custom code

x_length = length(xpos);
x_previous = xpos(1); % store the 1 previous x value in cache

number_of_columns= length(unique(xpos)); % counts how many columns are in the DICe result

newcolumn_index = zeros((number_of_columns+1),1); % will store the indices where a new columns begins 
newcolumn_index(1) = 1; % must preallocate that the first and last instance of a new column is the first eleement of xpos
newcolumn_index(end) = length(xpos); % final instance of a new column is at the end
newcolumn_index_counter = 1; % used for the loop below

% loop thru all x values, record the indices where a new column begins
for i = 1:x_length
    x_current = xpos(i);

    if x_current ~= x_previous
    newcolumn_index(newcolumn_index_counter+1) = i-1; % saves the index of the last x before the column shifts
    newcolumn_index_counter= newcolumn_index_counter +1;
    end

    x_previous = x_current;
end 

dL_column =zeros(1,number_of_columns); % stores the change in length total of the whole section in pixels
overall_strain = zeros(1,number_of_columns);
% loop thru each column calulate the difference between y values in pixels,
% calulate strain vector, multiply and sum ( discrete integrate)
for j = 1:(number_of_columns)

    dy_column = diff(ypos((newcolumn_index(j,1)):(newcolumn_index(j+1,1)),1)); % takes the difference between all y points in the column j
    strain_yy_column  = strain_yy((newcolumn_index(j,1)+1):newcolumn_index(j+1),1); % sacrifice the bottom point, assume dy = 0;
    %now in the column, loop through each value in the column, multiply it
    %by the span of pixels then sum it to the sum of all last products. 
    
    dL_column(j) = sum(dy_column.*(strain_yy_column));
    overall_strain(j) = dL_column(j)/(ypos(newcolumn_index(j+1))-ypos(newcolumn_index(j)));% divide nu
end
%absolute value removes negatives (an artifcact of how i processed data)
overall_strain_all = abs(overall_strain);

% nullify any inf or NaN values 
overall_strain_all(isinf(overall_strain_all)|isnan(overall_strain_all))= [];

% take the mean
overall_strain_all = mean(overall_strain_all);

% print each photo's overall strain
fprintf("Strain of the Gage Section is %f \n", overall_strain_all);
% add the strain from this photo to the global strain values vector
strain_values(m) = overall_strain_all;

end
% convert to table for copy pasting
strain_values = array2table(strain_values);
