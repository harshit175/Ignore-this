clear all
diff=0.00011646413; %thermal diffusivity

q=-2e5;

kc3=401;
rho3=8.92e3;
cp3=386;

time=.4-4e-2;
%n=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];  % terms of fourier series

%n=[1 2 3 4 5 6];

for count=1:100
    n(count)=count;
end

L=15e-3; %length along axis

x=0:0.1e-3:L; %points along axis


Term=x./L-0.5*(x./L).^2-1/3-diff*time/L^2;

for c=1:length(n)
    Term=Term+2/(n(c)*pi).^2*exp(-(n(c)^2*pi^2*diff*time/(L)^2))*cos((n(c)*pi*x)/(L)); %summation of terms
end

Tanalytical=25+(q*L/kc3)*Term;

 plot(1e3*x,Tanalytical,'*')