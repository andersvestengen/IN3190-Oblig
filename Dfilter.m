%Dirac filter funksjon 

function h = Dfilter(K)
    h = zeros(1,K+1);
    for j = 1:K+1
        h(j) = 1;
    end
   h = h*(1/(K+1));
end