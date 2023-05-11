function [Phi_theta] = der_phi(fattore_scala,lambda,D,positivi)

L=length(positivi);
Phi_theta=zeros(L,2);
aux=phi_nl(fattore_scala,lambda,D,positivi);
Phi_theta(:,1)=aux/fattore_scala;
Phi_theta(:,2)=aux/lambda;
for t=1:L
        if(D<=t)
            for s=1:t-D+1
                Phi_theta(t,2)=Phi_theta(t,2)-positivi(s)*lambda*(t-D-s)*exp(-lambda*(t-D-s));
            end
        end
end

