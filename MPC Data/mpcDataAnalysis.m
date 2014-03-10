clear all; close all; clc;
% plots of MPC data

load('neoData');

% EMoid vs H
scatter(MPCamors.EMoid,MPCamors.H)
hold on
scatter(MPCapollos.EMoid,MPCapollos.H)
scatter(MPCatens.EMoid, MPCatens.H)
title('EMoid vs H')
legend('Amors', 'Apollos', 'Atens')
xlabel('Earth MOID (AU)')
ylabel('H Magnitude')

xlim([0,0.05])
%%
% Q and q
figure;
scatter(MPCamors.q,MPCamors.Q)
hold on;
scatter(MPCapollos.q,MPCapollos.Q)
scatter(MPCatens.q,MPCatens.Q)