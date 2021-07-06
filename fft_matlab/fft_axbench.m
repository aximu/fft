N = 2048

for i= 1:N
    x(i) = i;
end

y=fft(x,N)