%%% extractFraction.m
%%% 16/12/2019
%%% Written by Soma Olasz
%%% 
%%% This script is created to calculate the runaway density
%%% of a NORSE calculation to Python languge.


function double = extractDistribution(NORSEobject)

    % calculate the runaway generation rate
    double =  NORSEobject.runawayFraction(end);
    

end