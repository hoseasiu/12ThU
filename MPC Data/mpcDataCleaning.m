% mpcDataCleaning.m
clear all; close all; clc;

load('neoData')
load('cleaned')
for i = 1:max(size(MPCamorscleaned(:,{'Designation'})))
   if isa(MPCamorscleaned.Designation{i},'double')
       MPCamorscleaned.Designation{i} = num2str(MPCamorscleaned.Designation{i});
   end
end
disp('done')

join(MPCamorscleaned,mpcLightcurveParameters,'Keys','Designation','KeepOneCopy','Period','KeepOneCopy','Variation')