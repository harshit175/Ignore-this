for count=1:iteration
residuenorm(count)=norm(Residue(j,count));
end
semilogy(residuenorm,'*')

xlabel('Number of iterations','fontsize',12)
ylabel('Norm of the Cure residue (log scale)','fontsize',12)

for count=2:iteration
convergenceratio(count)=residuenorm(count)/residuenorm(count-1)^2;
end