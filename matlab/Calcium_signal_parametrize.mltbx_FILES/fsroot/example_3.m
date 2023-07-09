%define path to file containing recordings
%here we use a relative path but an absolute path can be used if the file
%is not located in the same folder as the code
file='example.xls';

%we need to access the global namespace to specify that the algorithm 
%should plot the fitting results
global plt fig_save

%set the plt global to true so we can visaulize the fits
plt=true;
%save figures to png files
fig_save='png';

%execute characterization algorithm 
characterizeDocument(file);