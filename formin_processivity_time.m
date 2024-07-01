%Stochastic simulation of formin mediated filament polymerization
%Using experimentally obtained parameter - probability of formin dissociation per time (p_off-time)
%20% filament nuclaeted at t=0; constant nucleation rate for rest 80%

close all;
clear;
clc;

rng('shuffle');         %initialize and shuffle random number generator seed
format long g;

total_filaments = 12000; starting_filaments = 2400;         %number of total and starting number of filaments
colNames = {'0.5uM profilin', '1uM profilin','2.5uM profilin','5uM profilin','10uM profilin','15uM profilin'};     %column names corresponding to experimental profilin concentrations
elongation_rates = [14.93, 16.73, 18.36, 20.76, 16.53, 11.92];      %experimentally obtained elongation rates (subunits/sec) 
timelimit = 1200;       %simulation timelimit 
time_step=0.01;         %time increment value (sec)
off_rates = [0.001988, 0.001844, 0.002369, 0.003522, 0.002863, 0.001913]; %experimentally determined formin dissociation rates (1/sec) 
lengths = [];           %length counter
dwelltimes = [];        %dwelltime counter
stepcount = 0;

for z = 1:6
    steplimit = timelimit*elongation_rates(z);          %calculation of maximum filament length possible in the timelimit
    surviving = total_filaments;
    len = zeros (total_filaments, 1); 
    p_offtime = 1-exp(-off_rates(z)*time_step);        %calculation of probability of formin dissociation per time unit (p_offtime)
    time = zeros (total_filaments, 1);

    for x = 1:starting_filaments            %formin run for 20% filaments nucleated at t=0
        stepcount = 3;                      %resets stepcount and filament nucleates with 3 actin subunits
        while rand > p_offtime && stepcount < steplimit             %random number generated and compared to p_off-time for determining formin dissociation
            stepcount = stepcount + (time_step*elongation_rates(z));
        end
        if stepcount < steplimit            
            len (x,1) = stepcount;
            time (x,1) = stepcount/elongation_rates(z);
        else
            len (x,1) = 0;                  %filaments that have grown to maximum length are removed
            time(x,1) = 0;
        end
        surviving = surviving - 1;          %surviving filament counter
    end

    for t = 1:1200                  %time counter
        for i = 1: 8                %nucleation rate = 8 filaments per second
            stepcount = 3;          %resets stepcount and filament nucleates with 3 actin subunits
            while rand > p_offtime && stepcount < (steplimit - (elongation_rates(z)*(t-1)))     %random number generated and compared to p_off-time for determining formin dissociation
                stepcount = stepcount + (time_step*elongation_rates(z));
            end
            if stepcount < (steplimit - elongation_rates(z)*(t-1))
                len ((starting_filaments+8*(t-1)+i),1) = stepcount;             %recording the final filament length
                time ((starting_filaments+8*(t-1)+i),1) = stepcount/elongation_rates(z);
            else
                len ((starting_filaments+8*(t-1)+i),1) = 0;         %filaments that have grown to maximum length are removed
                time ((starting_filaments+8*(t-1)+i),1) = 0;
            end
            surviving = surviving -1 ;      %surviving filament counter
        end
    end
    len = nonzeros(len);
    lengths = [lengths, {len}];
    time = nonzeros(time);
    dwelltimes = [dwelltimes, {time}];

end

Run_lengths = array2table(lengths,'VariableNames',colNames);
Dwell_times = array2table(dwelltimes,'VariableNames',colNames);
