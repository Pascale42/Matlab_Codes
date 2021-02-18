% function X=GetmV(x, IR, varargin)
%
% Get signal in mVolts according to the Input Range (IR) to enter in Volts
% For +/- 2V >> IR=4
% By default (varargin), Amplification (A) is 1000



function X=GetmV(x, IR, varargin)

[A] = DefaultArgs(varargin,{1000});

% Marie:
X = x * (IR*1000) / A / 2^16;
% Eeg_mV = EegBrut * voltagerange_mV / amplification / 2^nbits
% 
% pour amplipex par exemple:
% - voltage range c'est -5/+5 V, donc total 10 V = 10 000 mV
% - amplification 400
% - nbits = 16



% X=(x./(2^16) /(IR/(A*1000))); % in mV
% % Vm = Vm./2^15*200;
% 
% % X=(x./(2^16))  * (IR*1000); % in mV
% X=x./(2^15)  * IR; % in mV
