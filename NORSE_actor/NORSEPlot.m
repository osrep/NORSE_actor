%%% NORSEPlot.m
%%% 09/01/2019
%%% Written by Soma Olasz
%%% 
%%% This script is created to be able to plot NORSE results from a NORSE
%%% object in Python language.
%%% Object must be a NORSE object. There can be only one arbitrary 
%%% variable, which must be a valid NORSE Plot object method in string
%%% format. 
%%% 
%%% Usage:
%%% 
%%%     NORSEPlot(object, plottingfunction)
%%% 

function NORSEPlot(object, varargin)

    % Check the number of input arguments. Make sure there is only one
    % varargin.
    if nargin ~= 2
        error('Invalid number of input arguments.');
    end
    
    % Plot results of a NORSE calculation
    object.plot.(varargin{1});
    

end