% Code by Gunjandeep Singh; HW 10

%% Problem 1

% We know that the tone output is byproduct of two frequency values.
% Assuming our initial setup is going to be similar to how it is presented
% in extra credit; with sampling frequency = 8000 Hz, duration = .25
% seconds, and so on, theoretically all we need to do is to fft of the
% tone, which should yield a graph with 4 small peaks. The first two peaks
% will be ones that correspond to the two frequencies that comprise the
% tone, while the latter two frquencies will be the ones that equate to our
% main two frequencies subtracted by sampling frequency of 8000 Hz. These
% latter two frequencies can be safely ignored, as they are conveying the
% same information as first two. Afterwards, we can basically rework the
% automatic detection code made available in this homework, and apply it to
% this problem in order to output the two frquencies that comprise our
% tone. We can then feed in that output to a series of if/elseif statements
% and output the single number that corresponds to the touchtone.

%The theoretical reason for why this works is because the entire
% spectrum is basically some particular signal that is replicated again and
% again in accordance with our sampling frequency. By doing discrete
% fourier transform, or some variation of it, we are able to extract out
% the signal that is getting repeated again and again.



%% Problem 2

clear
clc

% copy paste from extra credit. For testing purposes.

c1=1209;
c2=1336;
c3=1477;
r1=697;
r2=770;
r3=852;
r4=941;

fs = 8000; %sampling frequency
t=[0:(1/fs):.25];% 1/fs is sampling period. 0 to .25 is the length of tone in seconds
A1=.5;
p1=0;
A2=.5;
p2=0;
tone=A1*cos(2*pi*c2.*t+p1) + A2*cos(2*pi*r3.*t+p2);
%inputting freq values of c2 = 1336 and r3 = 852



s=decoding(tone);
% current problem. Even if I change my frequencies, my output continues to
% be 2. Not sure why. 







