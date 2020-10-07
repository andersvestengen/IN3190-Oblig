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

%% til filtrert seismogram1
figure('Name', "Filter2")
imagesc(offset1, t*1000, db(Filter1))
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
N = 300;
fs = 250;
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
%Identifiserer t_omega visuelt med plottet i 4a, mellom de to første
%toppene.
t_omega = (1.496 - 0.74); %t_omega2 - t_omega1
t2_omega = (3.16 - 2.432);
Ntrase = Filter1(:,1);
plot(t, Ntrase, 'black')
hold on
plot(t_omega, 0, 'redX')
plot((2.436), 0, 'blueX')
for i=2:4
    plot(t_omega*i, 0, 'redO')
end

for j=1:2
    plot((t2_omega*j + 2.436), 0, 'blueO')
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
%Hastighet = 1470 m/s

%% Oppgave 6a



%
t_d = t;
[u,v] = size(t);

VNMO_array = t./t;
VNMO_array = VNMO_array.*Hastighet;

dt = 1/fs;
Seis2 = nmocorrection2(t, dt, offset1, Filter1, VNMO_array);
%reflecttimes = sqrt(t.^2 + offset1.^2./Hastighet.^2);
%reflecttimes = reflecttimes.*(reflecttimes >= dt & reflecttimes <= t(end)-dt);
imagesc(offset1, t*1000, Seis2)
colormap gray;

%% Oppgave 6b

[u,v] = size(t);
Sp1 = Hastighet;
Sp2 = 2500;
VNMO_array = zeros(1, u);
VNMO_array(1:u*0.5) = Sp1;
VNMO_array(u*0.5:u) = Sp2;
hbedre = 1/50*ones(1, 50);
VNMO_dyn = konvin3190(VNMO_array, hbedre, 0);


Seis2 = nmocorrection2(t, dt, offset1, Filter1, VNMO_array');
Seis3 = nmocorrection2(t, dt, offset1, Filter1, VNMO_dyn');
%reflecttimes = sqrt(t.^2 + offset1.^2./Hastighet.^2);
%reflecttimes = reflecttimes.*(reflecttimes >= dt & reflecttimes <= t(end)-dt);
figure('name', 'uten filter')
imagesc(offset1, t*1000, Seis2)
colormap gray;
figure('name', 'med filter')
imagesc(offset1, t*1000, Seis3)
colormap gray;

%% 7a sediment-hastighet.
%verdier hentet visuelt fra plottet under oppgave 2c. Se figur i
%innlevering til 7a.
y1 = 1524;
y2 = 1652;
x1 = 2330;
x2 = 2680
hastighet2 = (x2-x1)/(y2-y1)

%% 7c
%konstruksjon av Filtrert seismogram
[p,l] = size(seismogram2);

SFilter1 = zeros(p,l);


for i=1:l
    SFilter1(:,i) = konvin3190(seismogram2(:,i),h1,0);
end

figure('Name', "Seismogram2 Filter1")
imagesc(offset2, t*1000, SFilter1)
xlabel("Tid [ms]")
ylabel("Offset [m]")
colormap gray;

figure('Name', "Seismogram2 Filter1 db")
imagesc(offset2, t*1000, db(SFilter1))
xlabel("Tid [ms]")
ylabel("Offset [m]")
colormap gray;

%% Oppgave 7d

%Bruker trigonometri og løser for rettvinklet trekant, hvor den minste
%siden er offset1[1] lengde, med den andre lengden basert på hastighet i
%vannet og hastighet i sedimentærlaget

LS1 = Hastighet*(0.818-0.112) %Refleksjon1 - direkte ankomst

LSediment1 = sqrt((LS1/2)^2 - 50^2)

LS2 = Hastighet*(2.464-0.112) %Refleksjon2 - direkte ankomst

LSediment2 = sqrt((LS2/2)^2 - 50^2)
