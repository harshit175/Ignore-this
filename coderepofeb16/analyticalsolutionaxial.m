%Analytical solution with one end at 200 C and the other end insulated

clear all
diff=.1;%0.00011646413; %thermal diffusivity
time=10e-1;
n=[1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43 45 47 49];  %odd terms of fourier series

L=1;%10e-3; %length along axis

x=0:0.2e-3:L; %points along axis


Term=0;
for c=1:length(n)
    Term=Term+4/(n(c)*pi)*exp(-(n(c)^2*pi^2*diff*time/(2*L)^2))*sin((n(c)*pi*x)/(2*L)); %summation of terms
end

Tanalytical=200+(50-200)*Term;