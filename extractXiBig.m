%%% extractDistribution.m
%%% 12/09/2019
%%% Written by Soma Olasz
%%% 
%%% This script is created to extract the Xi coordinate values
%%% from the last time step of a NORSE calculation to Python
%%% languge.
%%% 
%%% 

function double = extractDistribution(NORSEobject)

    % take the XiBig coordinate values of the NORSE object
    double = NORSEobject.grid.xiBig;
    

end