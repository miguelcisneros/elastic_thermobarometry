clc
close all
clear all
					

% CONSTANTS
% % % % % % % % % % % % % % % % % % % % % % % % 
h = 6.62607015*10^-34; % plancks constant J*s
k = 1.380649*10^-23; % boltzmann constant J/K
T = 25 + 273.15; % temperature Kelvin
c = 299792458; % speed of light m/s
% % % % % % % % % % % % % % % % % % % % % % % %


% READING FILES
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
testfiledir = 'C:\Users\cimiguel\Google Drive\primary\raman_data\raman_data_reduction_scripts\bose_einstein_reduction';
matfiles = dir(fullfile(testfiledir, '*.txt'));
nfiles = length(matfiles);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


% INITIALIZE VARIABLES
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
data = readtable(matfiles(1).name);
wavenumber = zeros(length(data{:,{'Var1'}}),1);
intensity_counts = zeros(length(data{:,{'Var2'}}),1);
bose_einstein_reduced_intensity = zeros(length(data{:,{'Var2'}}),1);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


for i = 1 : nfiles
      
    fid = fopen( fullfile(testfiledir, matfiles(i).name) );
    
    data = readtable(matfiles(i).name);
    
    
        for x = 1:length(data{:,{'Var1'}})
        
            raman_data = data{x,{'Var1','Var2'}};
    
            wavenumber = raman_data(1,1) * 100;
            intensity_counts = raman_data(1,2);
     
            bose_einstein_function = 1 - exp((-h * c * wavenumber) / (k * T));
    
            bose_einstein_reduced_intensity(x,1) = intensity_counts * bose_einstein_function;
        
        end
    
    
   reduced_intensity_counts = table(bose_einstein_reduced_intensity);

   data.Properties.VariableNames{'Var1'} = 'wavenumber';

   new_table = [data(:,'wavenumber') reduced_intensity_counts];
   
   new_name = erase(matfiles(i).name, '.txt');
   
   filename = strcat(new_name,'_reduced_intensity.txt');
   
   writetable(new_table,filename);
   
   fclose(fid);
   
end







 

    




 
 
 
 




