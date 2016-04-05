clear all
diff=0.00011646413; %thermal diffusivity

q=-2e7;

kc3=401;
rho3=8.92e3;
cp3=386;

L=100e-3; %length along axis

%x=0:0.05e-3:L; %points along axis
x=0:0.05e-3:L;

Tignition=200;

for count=1:100
    n(count)=count;  % 100 terms of fourier series taken
end


%n=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];  % terms of fourier series

%n=[1 2 3 4 5 6];

%time=000e-1:1e-1:141e-1;
time=013.7;
for i=1:length(time)
    Term=x./L-0.5*(x./L).^2-1/3-diff*time(i)/L^2;
    
i/length(time)
    
    
for c=1:length(n)
    Term=Term+2/(n(c)*pi).^2*exp(-(n(c)^2*pi^2*diff*time(i)/(L)^2))*cos((n(c)*pi*x)/(L)); %summation of terms
end

Tanalytical=25+(q*L/kc3)*Term;

for count=length(x):-1:1
    if Tanalytical(count)>Tignition
        xT(i)=1e3*x(count); %converted to mm
        break   
        
    end
    if count==1
        disp('Temperature less than');disp(Tignition);
    end

end
disp('blah')
plot(Tanalytical)
%pause(0.2)
end

 %plot(time,xT,'g*')