function [Yv] = phi_nl(fattore_scala,lambda,D,positivi)
L=length(positivi);
Yv=zeros(L,1);
for t=1:L
        if(D<=t)
            for s=1:t-D+1
                Yv(t)=Yv(t)+positivi(s)*fattore_scala*lambda*exp(-lambda*(t-D-s));
            end
        end
end

