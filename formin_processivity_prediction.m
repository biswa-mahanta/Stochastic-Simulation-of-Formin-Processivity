%Stochastic simulation of formin mediated filament polymerization
%Prediction of formin processivity from user defined FH1-pathway usage percentage
%Using experimentally obtained parameter - probability of formin dissociation per time (p_off-time)
%20% filament nuclaeted at t=0; constant nucleation rate for rest 80%

close all;
clear;
clc;

rng('shuffle');         %initialize and shuffle random number generator seed
format long g;

total_filaments = 12000; starting_filaments = 2400;         %number of total and starting number of filaments
prompt = "enter the FH1 flux (percentage)";                 %user input of the percentage of monomer addition through the FH1-delivery pathway
flux = input(prompt);
if flux>65
    off_rate = 0.0002*flux - 0.012;        %linear increase in formin dissociation rate when FH1 delivery is used for more than two-third monomers
else
    off_rate = 0.00191;                     %basal formin dissociation rate when FH1 delivery is used for less than two-third monomers 
end

colNames = {'5 subunits/s', '15 subunits/s','25 subunits/s','35 subunits/s'};     %column names corresponding to elongation rates used
elongation_rate = 5; %starting rate in subunits/sec
lengths = []; dwelltimes = [];
timelimit = 1200; time_step=0.01;
p_offtime = 1-exp(-off_rate*time_step);

for v = 1:4
    steplimit = timelimit*elongation_rate;
    len = zeros(total_filaments, 1);
    time = zeros (total_filaments, 1);
    surviving = total_filaments;
    stepcount = 0;

    for x = 1:starting_filaments            %formin run for 20% filaments nucleated at t=0
        stepcount = 3;                      %resets stepcount and filament nucleates with 3 actin subunits
        while rand > p_offtime && stepcount < steplimit             %random number generated and compared to p_off-time for determining formin dissociation
            stepcount = stepcount + (time_step*elongation_rate);
        end
        if stepcount < steplimit            
            len (x,1) = stepcount;
            time (x,1) = stepcount/elongation_rate;
        else
            len (x,1) = 0;                  %filaments that have grown to maximum length are removed
            time(x,1) = 0;
        end
        surviving = surviving - 1;          %surviving filament counter
    end    
    for t = 1:1200                  %time counter
        for i = 1: 8                %nucleation rate = 8 filaments per second
            stepcount = 3;          %resets stepcount and filament nucleates with 3 actin subunits
            while rand > p_offtime && stepcount < (steplimit - (elongation_rate*(t-1)))     %random number generated and compared to p_off-time for determining formin dissociation
                stepcount = stepcount + (time_step*elongation_rate);
            end
            if stepcount < (steplimit - elongation_rate*(t-1))
                len ((10*(t-1)+i),1) = stepcount;                    %recording the final filament length
                time ((10*(t-1)+i),1) = stepcount/elongation_rate;
            else
                len (10*(t-1)+i,1) = 0;         %filaments that have grown to maximum length are removed
                time ((10*(t-1)+i),1) = 0;
            end
            surviving = surviving -1 ;      %surviving filament counter
        end
    end
    len = nonzeros(len);
    lengths = [lengths, {len}];
    time = nonzeros(time);
    dwelltimes = [dwelltimes, {time}];
    elongation_rate = 5 + 10*(v);
end

Run_lengths = array2table(lengths,'VariableNames',colNames);
Dwell_times = array2table(dwelltimes,'VariableNames',colNames);
