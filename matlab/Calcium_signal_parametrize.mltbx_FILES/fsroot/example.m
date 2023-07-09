%define path to file containing recordings
%here we use a relative path but an absolute path can be used if the file
%is not located in the same folder as the code
file='example.xls'



%execute characterization algorithm 
characterizeDocument(file);