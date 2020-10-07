%Konvolusjons funksjon for obligen

function y = konvin3190(x, h, ylen)
% Oppgave 1a)
      M = length(h);
      N = length(x);
      y = zeros(1,(M + N - 1));
      for m = 1:M
          for n = 1:N
              k = n+m-1;
                  y(k) = y(k)+h(m)*x(n);
          end
      end
      if ylen == 0
      % for len(y) = len(x)
      % removing the (len(h)-1)/2 elements to return accurate x-length.
      et1 = ceil((M-1)/2 + 1);
      et2 = floor((M-1)/2);
      y = y(et1:length(y)-et2);
      end
end