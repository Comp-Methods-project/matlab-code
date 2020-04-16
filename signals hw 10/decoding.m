function output = decoding(x);
fs = 8000; % sampling frequency; given
N=length(x); % number of data points
df = fs/N; % frequency increment
freq = [0:N-1]*df; % frequency mapping vector
tone = abs(fft(x)); % fast fourier transform

tone2=tone; % working variable
tone2(1:10)=0; % zero out baseline peak
tone2(round(N/2):N)=0;%discarding second half
peaks = maxk(tone2,2);
% maxk finds specified number of max values in a vector. Here the specified amount is 2.
% the first maxk value will always be the biggest in the vector
% the second maxk value will always be the second biggest one in the vector

peak1=peaks(:,1); % extracting out index value for finding 1st max frequency
peak2=peaks(:,2); % extracting out index value for finding 2nd max frequency
plot(freq,tone)
idx1 = find(tone2==peak1);% find index location of 1st peak
idx2 = find(tone2==peak2);% find index location of the 2nd peak

format short g
v1 = freq(idx1) % using index location, this outputs the 1st max freq
format short g
v2 = freq(idx2) % using index location, this outputs the 2nd max freq


if 1333<v1<1338 & 694<v2<700
    output = 2
elseif 1207<v1<1211 & 694<v2<700
    output = 1
elseif 1474<v1<1480 & 694<v2<700
    output = 3
elseif 1207<v1<1211 & 768<v2<772
    output = 4
elseif 1333<v1<1338 & 768<v2<772
    output = 5
elseif 1474<v1<1480 & 768<v2<772
    output = 6
elseif 1207<v1<1211 & 850<v2<854
    output = 7
elseif 1333<v1<1338 & 850<v2<854
    output = 8
elseif 1474<v1<1480 & 850<v2<854
    output = 9
elseif 1333<v1<1338 & 939<v2<943
    output = 0
else output = 99
end
