%%% extractDistribution.m
%%% 12/09/2019
%%% Written by Soma Olasz
%%% 
%%% This script is created to extract the p coordinate values
%%% from the last time step of a NORSE calculation to Python
%%% languge.
%%% 
%%% 

function double = extractDistribution(NORSEobject)

    % take the pBig coordinate values of the NORSE object
    double = NORSEobject.grid.pBig;
    

end