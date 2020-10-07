%Oppgave 1c plott av DTFT for X og Y 
fs = 100;
N = 400;
t = 0:1/fs:5;
f1 = 10;
f2 = 20;

x = sin(2*pi*f1*t) + sin(2*pi*f2*t);
h = Dfilter(4);

y = konvin3190(x, h, 0);

[X, f1] = frekspekin3190(x, N, fs);
%[T, f] = freqz(x,1,N,fs)

[Y, f2] = frekspekin3190(y, N, fs);
%plot(t, abs(X))

%plot(f1, abs(X))
%hold on 
plot(f2, abs(Y))
title("DTFT signal")
xlabel("Frequencies [Hz]")
ylabel("Radians [pi/N]")
%legend("X", "Y")
legend("Y")
%% Oppgave 1c Plot av DTFT for filteret
h = Dfilter(4);
N = 300;
fs  = 100;
[H, hf] = frekspekin3190(h, N, fs);
plot(hf, abs(H))
title("DTFT H")
xlabel("Frequencies [Hz]")
ylabel("Radians [pi/N]")
legend("H")

%% Oppgave 2 a
N = 300;
fs = 250;

[F1, f1] = frekspekin3190(h1, N, fs);
F1log = 20*log10(abs(F1));

[F2, f2] = frekspekin3190(h2, N, fs);
F2log = 20*log10(abs(F2));

plot(f1, F1log)
title("DTFT Filter h1 / h2")
hold on
plot(f2, F2log)
xlabel("Frequencies [Hz]")
ylabel("decibel [dB]")
legend("F1", "F2")

%% Oppgave 2b signaler

N = 300;
fs = 100;

Seisstart = seismogram1(1:400,1);
W = tukeywin(length(Seisstart), 0.25);
figure('Name', 'Nærtrase n=200')

plot(t(1:400), Seisstart)
hold on
plot(t(1:400), W.*Seisstart)
legend("uten vindu", "Med vindu")
xlabel("tid [ms]")
ylabel("utslag")



%% Oppgave 2b DTFT Nærtraser med og uten vinduer i dB.
N = 300;
fs = 250;
%fk = 1/(t(2)-t(1))
Seisstart = seismogram1(1:400,1);
SeisFreqs = (W.*Seisstart);
[Wf, wf9] = frekspekin3190(W, N, fs);
[Fw, fw] = frekspekin3190(SeisFreqs, N, fs);
[Fs, fs] = frekspekin3190(Seisstart, N, fs);
Fwlog = 20*log10(abs(Fw));
Fslog = 20*log10(abs(Fs));

figure('Name', 'DTFT with decibel')
plot(fw, Fwlog)
hold on 
plot(fs, Fslog)
legend("med vindu", "uten vindu")
xlabel("frekvenser [Hz]")
ylabel("Utslag [dB]")

%% Oppgave 2c 

[m,n] = size(seismogram1);

Filter1 = zeros(m,n);
Filter2 = zeros(m,n);

for i=1:n
    Filter1(:,i) = konvin3190(seismogram1(:,i),h1,0);
    Filter2(:,i) = konvin3190(seismogram1(:,i),h2,0);
end



figure('Name', "uten filter")
imagesc(seismogram1(:,1))
xlabel("Tid [ms]")
ylabel("Offset [m]")
colormap gray;
figure('Name', "Filter1")
imagesc(Filter1(:,1))
xlabel("Tid [ms]")
ylabel("Offset [m]")
colormap gray;
figure('Name', "Filter2")
imagesc(Filter2(:,1))
xlabel("Tid [ms]")
ylabel("Offset [m]")
colormap gray;

%% Oppgave 3a
N = 300;
fs = 1/(t(2)-t(1))
Slice1 = Filter1(20:150,30);
figure('Name', "Direkte ankomst")
plot(t(20:150),Slice1)
xlabel("Tid [ms]")
ylabel("Offset [m]")

[Sl, sf] = frekspekin3190(Slice1, N, fs);
figure('Name', "DTFT direkte ankomst")
plot(sf, Sl)
xlabel("frekvenser [Hz]")
ylabel("Utslag")
offset1(1, 30)

%% Oppgave 3b

Slice1 = Filter1(20:150,30);
W = tukeywin(length(Slice1), 0.25);
FSS = (W.*Slice1);
[Dt, ft] = frekspekin3190(FSS, N, fs);
[Ut, ut] = frekspekin3190(Slice1, N, fs);
FSSFT = 20*log10(abs(Dt));
Utlog = 20*log10(abs(Ut));
figure('Name', 'dB FFT av vindu uten vindu med direkte ankomst')
plot(ut, Utlog)
hold on
plot(ft, FSSFT)
xlabel("frekvenser [Hz]")
ylabel("Utslag [dB]")
legend('uten vindu', 'med vindu')




%% Oppgave 3c

MaxF = find(FSSFT == max(FSSFT(:)));
fes = ut(MaxF)
c = 3000;
Lb = c/fes;
h = (1/8)*Lb

%% Oppgave 4a

Ntrase = Filter1(:,1);
plot(t, Ntrase)


%% Oppgave 4b
t_omega = 0.624;
t2 = t_omega*2;
t3 = t_omega*3;
Ntrase = Filter1(:,1);
plot(t, Ntrase)

for i=1:4
    hold on
    plot(t_omega*i, 0, 'redX')
end


%% Oppgave 5

Slice1 = Filter1(20:150,30);
Slice2 = Filter1(50:180,50);
plot(t(20:150),Slice1)
hold on
plot(t(50:180),Slice2)
off1 = offset1(1, 30)
off2 = offset1(1, 50)

% Hentet hastigheter grafisk ved å se på toppene til toppene i plottet.
Hastighet = (off2 - off1)/(0.476 - 0.340)

%% Oppgave 6a
dt = 1/fs;
Seis2 = nmocorrection2(t, dt, offset1, Filter1, Hastighet);
%reflecttimes = sqrt(t.^2 + offset1.^2./Hastighet.^2);
%reflecttimes = reflecttimes.*(reflecttimes >= dt & reflecttimes <= t(end)-dt);
imagesc(Seis2)
colormap gray;

%% Oppgave 6b
%Observer at 'skyggene' er rette her, som er refleksjonene fra et
%sedimentlag
THastighet = 2500
dt = 1/fs;
Seis2 = nmocorrection2(t, dt, offset1, Filter1, THastighet);
%reflecttimes = sqrt(t.^2 + offset1.^2./Hastighet.^2);
%reflecttimes = reflecttimes.*(reflecttimes >= dt & reflecttimes <= t(end)-dt);
imagesc(Seis2)
colormap gray;

%% Oppgave 7



