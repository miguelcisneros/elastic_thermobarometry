clc
close all
clear all 

% Values from Wang et al., 2015

% C11 (GPa)	(GPa)	C33 (GPa)	(GPa)	C12 (GPa)	(GPa)	C13 (GPa)	(GPa)	C14 (GPa)	(GPa)	C44 (GPa)	(GPa)
% 86.6      0.3     106.4       1.2     6.74        0.09	12.4        1.6     17.8        0.3     58          0.7

% sigma1	c11	c12	c13	c14	c15	c16	?1
% sigma2	c21	c22	c23	c24	c25	c26	?2
% sigma3	c31	c32	c33	c34	c35	c36	?3
% sigma4	c41	c42	c43	c44	c45	c46	?4
% sigma5	c51	c52	c53	c54	c55	c56	?5
% sigma6	c61	c62	c63	c64	c65	c66	?6

% stiffness matrix
%       x1      x2      x3      x4      x5      x6
% x1	86.6	6.74	12.4			
% x2	6.74	86.6	12.4			
% x3	12.4	12.4	106.4			
% x4                            58		
% x5						
% x6						

data = readtable('input_raman_data.xlsx','ReadRowNames',true);

p_incl_128_MPa = zeros(length(data{:,{'shift_464'}}),1);
p_incl_206_MPa = zeros(length(data{:,{'shift_464'}}),1);
p_incl_464_MPa = zeros(length(data{:,{'shift_464'}}),1);
p_incl_128_MPa_error = zeros(length(data{:,{'shift_464'}}),1);
p_incl_206_MPa_error = zeros(length(data{:,{'shift_464'}}),1);
p_incl_464_MPa_error = zeros(length(data{:,{'shift_464'}}),1);
mean_hyrostatic_p_MPa = zeros(length(data{:,{'shift_464'}}),1);

mean_stress_MPa = zeros(length(data{:,{'shift_464'}}),1);
mean_stress_error_MPa = zeros(length(data{:,{'shift_464'}}),1);

difference_MPa = zeros(length(data{:,{'shift_464'}}),1);


for i = 1:length(data{:,{'shift_464'}})
    
%     RAMAN SHIFTS

    raman_shift = data{i,{'shift_128','shift_206','shift_464'}};
    raman_shift_error = data{i,{'shift_128_error','shift_206_error','shift_464_error'}};
 
    raman_shift_128 = raman_shift(1,1);
    raman_shift_206 = raman_shift(1,2);
    raman_shift_464 = raman_shift(1,3);
    
    raman_shift_128_error = raman_shift_error(1,1);
    raman_shift_206_error = raman_shift_error(1,2);
    raman_shift_464_error = raman_shift_error(1,3);
    
    p_incl_128_MPa(i,1) = (1314.3 * raman_shift_128 + 47.5 * (raman_shift_128 ^ 2)) / 10;
    p_incl_206_MPa(i,1) = (307 * raman_shift_206 + 4.76 * (raman_shift_206 ^ 2)) / 10;
    p_incl_464_MPa(i,1) = (1094.5 * raman_shift_464 + 4.204 * (raman_shift_464 ^ 2))/10;
    
    %Thomas and Spear 2018
    p_incl_128_MPa_error(i,1) = (sqrt(((raman_shift_128_error) * ((2 * 47.5) * raman_shift_128 + 1314.3 ))^2 + (3.1561 * (raman_shift_128 ^ 2)) ^ 2 + (29.6 * raman_shift_128) ^ 2))/10;
    %Steele and Ashley
    p_incl_206_MPa_error(i,1) = (sqrt(((raman_shift_206_error) * ((2 * 4.76) * raman_shift_206 + 307))^2 + (0.21 * (raman_shift_206 ^ 2)) ^ 2 + (7 * raman_shift_206) ^ 2))/10;
    %Ashley 2013
    p_incl_464_MPa_error(i,1) = (sqrt(((raman_shift_464_error) * ((2 * 4.204) * raman_shift_464 + 1094.5))^2 + (0.81 * (raman_shift_464 ^ 2)) ^ 2 + (12 * raman_shift_464) ^ 2))/10;
    
    mean_hyrostatic_p_MPa(i,1) = (p_incl_128_MPa(i,1) + p_incl_206_MPa(i,1) + p_incl_464_MPa(i,1))/3;
    
%     STRAINS
    strains = data{i,{'strain1_strain2','strain3'}};

    strain_1_2 = strains(1,1);
    strain_3 = strains(1,2);

        C11 = 86.6;
        C12 = 6.74;
        C13 = 12.4;
        %C14 = 17.8;
        C33 = 106.4;
        %C44 = 58;

        C21 = C12;
        C22 = C11;
        C23 = C13;
        
        C31 = C13;
        C32 = C23;
        
        C11_error = 0.3;
        C12_error = 0.09;
        C13_error = 1.6;
        %C14_error = 0.3;
        C21_error = 0.09;
        C22_error = 0.3;
        C23_error = 1.6;
        C31_error = 1.6;
        C32_error = 1.6;
        C33_error = 1.2;
        %C44_error = 0.7;

    strain_1 = strain_1_2/2;
    strain_2 = strain_1_2/2;

    stress_1 = C11 * strain_1 + C12 * strain_2 + C13 * strain_3;    
    stress_2 = C21 * strain_1 + C22 * strain_2 + C23 * strain_3;
    stress_3 = C31 * strain_1 + C32 * strain_2 + C33 * strain_3;
    
    stress_1_error = sqrt((strain_1 * C11_error)^2 + (strain_2 * C12_error)^2 + (strain_3 * C13_error)^2);
    stress_2_error = sqrt((strain_1 * C21_error)^2 + (strain_2 * C22_error)^2 + (strain_3 * C23_error)^2);
    stress_3_error = sqrt((strain_1 * C31_error)^2 + (strain_2 * C32_error)^2 + (strain_3 * C33_error)^2);

    mean_stress_MPa(i,1) = (-(2 * stress_1 + stress_3)/3) * 1000;
    
    mean_stress_error_MPa(i,1) = (sqrt(((2 ^ (1/3))/(3 * (-stress_1) ^ (2 / 3)) * stress_1_error) ^ 2 +  (1/(3 * (-stress_3) ^ (2 / 3)) * stress_3_error) ^ 2)) * 1000;
    
    difference_MPa(i,1) = mean_stress_MPa(i,1) - mean_hyrostatic_p_MPa(i,1);
   
end

 p_incl_128_MPa = table(p_incl_128_MPa);
 p_incl_206_MPa = table(p_incl_206_MPa);
 p_incl_464_MPa = table(p_incl_464_MPa);
 
 p_incl_128_MPa_error = table(p_incl_128_MPa_error);
 p_incl_206_MPa_error = table(p_incl_206_MPa_error);
 p_incl_464_MPa_error = table(p_incl_464_MPa_error);
 
 mean_hyrostatic_p_MPa = table(mean_hyrostatic_p_MPa);

 mean_stress_MPa = table(mean_stress_MPa);
 mean_stress_error_MPa = table(mean_stress_error_MPa);
 
 difference_MPa = table(difference_MPa);

 new_table = [data p_incl_128_MPa p_incl_206_MPa p_incl_464_MPa mean_hyrostatic_p_MPa mean_stress_MPa difference_MPa p_incl_128_MPa_error p_incl_206_MPa_error p_incl_464_MPa_error mean_stress_error_MPa];
 
 filename = 'inclusion_pressure_calculations_output.xlsx';
 
 writetable(new_table,filename);
 
 
 
 




