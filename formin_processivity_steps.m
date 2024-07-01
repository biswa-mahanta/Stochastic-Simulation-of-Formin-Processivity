%Stochastic simulation of formin mediated filament polymerization
%Using experimentally obtained parameter - probability of formin dissociation per step (p_off-step)
%20% filament nuclaeted at t=0; constant nucleation rate for rest 80%

close all;
clear;
clc;

rng('shuffle');             %initialize and shuffle random number generator seed
format long g;

total_filaments = 12000; starting_filaments = 2400;          %number of total and starting number of filaments
colNames = {'0.5uM profilin', '1uM profilin','2.5uM profilin','5uM profilin','10uM profilin','15uM profilin'};     %column names corresponding to experimental profilin concentrations
elongation_rate = [14.93, 16.73, 18.36, 20.76, 16.53, 11.92];          %experimentally obtained elongation rates (subunits/sec) 
timelimit = 1200;       %simulation timelimit
time_step=0.01;         %time increment value (sec)
k_offstep = [0.00021, 0.00018, 0.00017, 0.00018, 0.00021, 0.00024];     %experimentally determined formin dissociation rates per step (1/step) 
lengths = [];           %length counter
dwelltimes = [];         %dwelltime counter
stepcount = 0;         


for z = 1:6
    steplimit = timelimit*elongation_rate(z);               %calculation of maximum filament length possible in the timelimit
    surviving = total_filaments;
    len = zeros (total_filaments, 1); 
    p_offstep = 1 - exp(-k_offstep(z)*1);                              %calculation of probability of formin dissociation per step (p_offstep)
    time = zeros (total_filaments, 1);

    for x = 1:starting_filaments            %formin run for 20% filaments nucleated at t=0
        stepcount = 3;                      %resets stepcount and filament nucleates with 3 actin subunits
        while rand > p_offstep && stepcount < steplimit             %random number generated and compared to p_off-seep for determining formin dissociation
            stepcount = stepcount + 1;
        end
        if stepcount < steplimit
            len (x,1) = stepcount;
            time (x,1) = stepcount/elongation_rate(z);
        else
            len (x,1) = 0;                  %filaments that have grown to maximum length are removed
            time(x,1) = 0;
        end
        surviving = surviving - 1;          %surviving filament counter
    end

    for t = 1:1200              %time counter
        for i = 1: 8            %nucleation rate = 8 filaments per second
            stepcount = 3;      %resets stepcount and filament nucleates with 3 actin subunits
            while rand > p_offstep && stepcount < (steplimit - (elongation_rate(z)*(t-1)))      %random number generated and compared to p_off-time for determining formin dissociation
                stepcount = stepcount + 1;
            end
            if stepcount < (steplimit - elongation_rate(z)*(t-1))
                len ((starting_filaments+8*(t-1)+i),1) = stepcount;             %recording the final filament length
                time ((starting_filaments+8*(t-1)+i),1) = stepcount/elongation_rate(z);
            else
                len ((starting_filaments+8*(t-1)+i),1) = 0;         %filaments that have grown to maximum length are removed
                time ((starting_filaments+8*(t-1)+i),1) = 0;
            end
            surviving = surviving -1 ;          %surviving filament counter
        end
    end
    len = nonzeros(len);
    lengths = [lengths, {len}];
    time = nonzeros(time);
    dwelltimes = [dwelltimes, {time}];

end

Run_lengths = array2table(lengths,'VariableNames',colNames);
Dwell_times = array2table(dwelltimes,'VariableNames',colNames);
