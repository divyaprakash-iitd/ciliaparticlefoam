clear; clc; close all;

% Description: This script checks the proper function of the softcilia
% library written in fortran

CFile = dir(strcat('C_0','*'));


%% Plot cilia
C = load(CFile(1).name)';
hold on
for i = 1:2:size(C,2)
    plot(C(:,i),C(:,i+1),'k.','linewidth',2)
end
daspect([1,1,1])