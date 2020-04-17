function output = ecDecoding(x)
fs=8000; % sampling frequency, given
fcol = [1209 1336 1477]; % column frequency of numbers, given
frow = [697 770 852 941]; % row frequency of numbers, given
t = [0:(1/fs):.25]; % .25 is duration provided in assignment
f_total = [frow fcol];
col_container=0; % values from dot product calculation will be placed here
row_container=0; % values from dot product calculation will be placed here

[row,~]=size(x); % this is necessary because test tone has multiple rows. So the entire code needs to be inside a while loop with respect to rows

output = zeros(1,row);
a=1;

while a<=row
    
    for i = 1:length(frow) % loop for calculating dot product for row frequencies
        basis1=cos(2*pi*frow(i)*t); % basis function
        dotProduct1 = sum(x(a,:).*basis1); % dot product with each basis function
        row_container(i)=dotProduct1; % placing values in vector container
    end
    for j = 1:length(fcol) % loop for calculating dot product for column frequencies
        basis2=cos(2*pi*fcol(j)*t); % our basis func.
        dotProduct2 = sum(x(a,:).*basis2); % dot product with basis function
        col_container(j)=dotProduct2; % placing values into vector container
    end
    
freq_container = 0; % both values from row and column vectors will go in this        
        
    maxRow=max(abs(row_container)); % This is the trick to this assignment. Finding max value from row vector
    maxCol=max(abs(col_container)); % Finding max value from this vector
   
   
    cumulDotProduct = [abs(row_container) abs(col_container)]; % row and column vectors in one single container
    
    for k = 1:length(cumulDotProduct); % this loop is designed to extract and place max value in the cumulative container, while placing a value of zero everywhere else
        if cumulDotProduct(k) == maxRow | cumulDotProduct(k)==maxCol;
            freq_container(k)=f_total(k);
        else
            freq_container(k)= 0;
        end
    end

   l = find(freq_container); % this setup allows our freq container to only report values that correspond to initial frequencies used
   if l == [1 5];
       output(a) = '1';
   elseif l == [1 6];
       output(a) = '2';
   elseif l == [1 7];
       output(a) = '3';
   elseif l == [2 5];
       output(a) = '4';
   elseif l == [2 6];
       output(a) = '5';
   elseif l == [2 7];
       output(a) = '6';
   elseif l == [3 5];
       output(a) = '7';
   elseif l == [3 6];
       output(a) = '8';
   elseif l == [3 7];
       output(a) = '9';
   elseif l == [4 5];
       output(a) = '*';
   elseif l == [4 6];
       output(a) = '0';
   elseif l == [4 7]
       output(a) = '#';
   end
   a = a+1;
end

output = char(output);
