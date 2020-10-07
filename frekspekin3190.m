%DTFT funksjon for obligen. mesteparten hentet fra lagoritme i Manolakis
%Diigital signalbehandling bok.

function [X, f]  = frekspekin3190(x, N, fs)
% Oppgave 1b)
    lenx = length(x);
    X = zeros(1,N);
    omega = ((0:N-1)/N)*pi;   %linspace(0, pi, N);
    for k = 1:N
        for n = 1:lenx
        X(k) = X(k) + x(n)*exp(-1j*omega(k)*(n-1));
        end
    end
   f = (omega/pi)*(fs/2);
end