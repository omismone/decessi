function [decessi] = convoluz(positivi,D,fattore_scala,lambda)
positivi=positivi(:);
n=length(positivi);
zeri=zeros(D,1);
L=n-D;
T=[0:L-1]';
g=lambda*exp(-lambda*T);
g=g/sum(g);
g=g*fattore_scala;
g=[zeri;g];
decessi = conv(g,positivi);
decessi=decessi(1:n);
end

