function [Phi_theta] = rapp_incr(fattore_scala,lambda,D,positivi,delta1,delta2)
    L=length(positivi);
    Phi_theta=zeros(L,2);
    Phi_theta(:,1)=(phi_nl(fattore_scala+delta1,lambda,D,positivi)-phi_nl(fattore_scala,lambda,D,positivi))/delta1;
    Phi_theta(:,2)=(phi_nl(fattore_scala,lambda+delta2,D,positivi)-phi_nl(fattore_scala,lambda,D,positivi))/delta2;
end

