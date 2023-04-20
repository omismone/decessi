function [ssr] = calculateSSR(positivi,D,fattore_scala,lambda,media7ggD)
ingressi=table2array(positivi((247-D):(397-D),3));
conv=convoluz(ingressi,D,fattore_scala,lambda);
ssr=(media7ggD - conv)' * (media7ggD - conv);
end

