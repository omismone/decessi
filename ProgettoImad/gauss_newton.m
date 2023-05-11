function [theta_k] = gauss_newton(Y,fattore_scala,lambda,D,positivi)
    syms fs_var l_var k_var u_var;

    phi_matrix=zeros(151,1);
    phi_matrix_lin=zeros(151,2);

    theta_k=[fattore_scala;lambda];

    %creo phi di theta k
    for k=1:151
        ingresso=table2array(positivi((247-D-k),3));
        phi_matrix(k)=(fattore_scala*exp(-lambda*k))*lambda*ingresso;   %calcolo phi coi valori attuali (k)
    end

    %mat lineare
    
    %*u_var(247-D-k_var)
    
    tmp1=diff(fs_var*exp(-l_var*k_var)*l_var,fs_var);                   %derivate parziali rispetto t1,
    tmp2=diff(fs_var*exp(-l_var*k_var)*l_var,l_var);                    %t2
    for k=1:151
        ingresso=table2array(positivi((247-D-k),3));
        phi_matrix_lin(k,1)=subs(tmp1*ingresso,[fs_var,l_var,k_var],[fattore_scala,lambda,k]);  %phi lineare calcolata in
        phi_matrix_lin(k,2)=subs(tmp2*ingresso,[fs_var,l_var,k_var],[fattore_scala,lambda,k]);  % theta k
    end

    theta_k=theta_k+(inv(phi_matrix_lin'*phi_matrix_lin))*phi_matrix_lin'*(Y-phi_matrix); %creo i valori con k+1
end

