% This application works with voice signals. It calculates:
% - The pitch period obtained directly from the signal
% - The periodogram obtained from the DFT of a segment
% - AR spectrum (also known as LPC spectrum) comparison applying different methods (cross-autocorrelation, covariance)
% - Cepstrum FFT signal calculation.
% - Vowels recognition using the cepstrum LPC calculation (works better that the cepstrum FFT)

clear;
%% Cargamos las seÃ±ales
load senales/a.asc;
load senales/e.asc;
load senales/i.asc;
load senales/o.asc;
load senales/u.asc;
load senales/x.asc;


plot(e);
xlabel('N');
ylabel('Señal e');
title('Señal e.');

fs=8000;
N=length(e);
pause;

% sound(e,8*10^3);
% sound(a,8*10^3);
% sound(i,8*10^3);
% sound(o,8*10^3);
% sound(u,8*10^3);
% sound(x,8*10^3);

%% 2.1. PERIODOGRAMA
% k=(abs(fft(e).^2))/length(e);

PSDe=(abs(fft(e).^2))/length(e);
frecuencia=(0:128)*fs/N;
plot(frecuencia,10*log10(PSDe(1:N/2+1)))
xlabel('w');
ylabel('PSD Periodograma');
title('PSD Periodograma Pitch=157Hz Formantes=469, 2250, 3531');
pause;

% - Rizado debida a la excitación. Frecuencia de pitch --> diferencia entre
% dos picos consecutivos del rizado

% - Formates: picos suaves que corresponden a las frecuencias de resonancia
% del tracto vocal. 

%% 2.2. ESTIMACION AR

%% 1.-PERIODOGRAMA(A partir de la autocorrelacion)
rx=xcorr(e,'biased');
plot(rx);
RX=fft(rx); %Periodograma AR
h=10*log10(abs(RX));
frecuencia1=(0:255)*fs/(2*N);
plot(frecuencia1,h(1:256))
xlabel('w');
ylabel('PSD Periodograma Rx');
title('PSD Periodograma Rx Pitch=157Hz Formantes=485,2234,3531');
pause;

%Caracteristicas de la estimación:
% - Estimación asintoticamente sin desplazamiento: E[Pper(w)] =
% F[E[rx(k)(biased)]] = F[E[wB(k)*rx(k)]] = WB(w)*Px(w)--> Px(w)(N-->inf).
% Resolución del periodograma: 0.89*(2*pi/N)(ancho ventana Barlett 3dB)

% - Estimación no consistente. Ejemplo: ruido blanco gaussiano.
% VAR[Pper(w)] = Px(w)^2 = (sigmax^2)^2

%% 2.-ESPECTRO AR(12)POR EL METODO AUTOCORRELACION

frecuencia=(0:128)*fs/N;

% Matriz Rx
for i=0:1:11;
    aux(i+1)=Autoc(i,e);
end

Rx1=toeplitz(aux,aux);

% Matriz C
for i=1:1:12;
    c1(i)=Autoc(i,e);
end

c1=c1';

% Coeficientes
a_xcorr=inv(Rx1)*(-c1);
a_xcorr=a_xcorr';

suma_xcorr=0;

for k=1:1:12;
   suma_xcorr=suma_xcorr+a_xcorr(k)*Autoc(k,e);
end

b0_xcorr=Autoc(0,e)+suma_xcorr;
PSDFormula_xcorr=b0_xcorr./(abs(fft([1,a_xcorr],256))).^2;

% Coeficientes con la función LPC
a_lpc=lpc(e,12); % Ya lleva el coeficiente a0=1 inclusive  calculado
suma_lpc=0;

for k=1:1:12;
   suma_lpc=suma_lpc+a_lpc(k+1)*Autoc(k,e);
end

b0_lpc=Autoc(0,e)+suma_lpc;
PSDFormula_lpc=b0_lpc./(abs(fft(a_lpc,256))).^2;

plot(frecuencia,10*log10(PSDFormula_xcorr(1:N/2+1)),'red','Linewidth',2)
xlabel('w');
ylabel('PSD');
hold on;
plot(frecuencia,10*log10(PSDFormula_lpc(1:N/2+1)))
title('PSD Autocorrelacion vs LPC. Formantes: 500, 2219, 3562');
legend('PSD autocorrelacion','PSD LPC');
pause;

%% 3.-ESPECTRO AR(12)POR EL METODO COVARIANZA

% Matriz Rx2
for k=1:1:12;
    aux1(k)=AutocD(k,1,12,e);
end

Rx2=toeplitz(aux1,aux1);

% Matriz C2
for l=1:1:12;
    c2(l)=AutocD(0,l,12,e);
end

c2=c2';

% Coeficientes 
a2=inv(Rx2)*(-c2);
a2=a2';

suma2=0;
for k=1:1:12;
   suma2=suma2+a2(k)*Autoc(k,e);
end
b0_2=Autoc(0,e)+suma2;

% NO SE SI b0 SE PUEDE CALCULAR ASI PARA EL METODO COVARIANZA

PSDFormula2=b0_2./(abs(fft([1,a2],256))).^2;
figure;
plot(frecuencia,10*log10(PSDFormula2(1:N/2+1)))
xlabel('w');
ylabel('PSD Covarianza');
title('PSD Covarianza');

pause;

plot(frecuencia,10*log10(PSDe(1:N/2+1)),'r',frecuencia1,h(1:256),'y',frecuencia,10*log10(PSDFormula11(1:N/2+1)),'bd',frecuencia,10*log10(PSDFormula11(1:N/2+1)),'g.',frecuencia,10*log10(PSDFormula2(1:N/2+1)),'k')
legend('Periodograma','PeriodogramaRx','Autocorrelacion','LPC','Covarianza');
title('Comparacion de PSDs')
xlabel('w')
ylabel('Comparacion de PSDs')

pause;


%% 3. ANALISIS HOMOMORFICO


Pe=abs(RX);

ce1=ifft(log(Pe));
figure;
plot((0:127),ce1(1:128),'red')
xlabel('frecuencia')
ylabel('Cepstrum periodograma')
title('Cepstrum periodograma & Cepstrum LPC. Pitch=50 (pico despues n=20)')
ce2=ifft(log(PSDFormula_lpc));
hold on;
plot((0:127),ce2(1:128));
legend('Cepstrum periodograma','Cepstrum LPC');

pause;

% - Cepstrum FFT: cx(n) = ch(n)+cu(n). n dominio temporal -> Cuefrecuencia
% (para diferenciar claramente el cepstrum de la señal original). Dos
% componentes: una debida al filtro vocal(correlaciones de retardo corto) y
% a la excitación(correlaciones de retardo largo), respectivamente. Fácil
% de separar dado que la excitación afecta esencialmente a componentes de
% alta frecuencia, mientras que le filtro contribuye a las bajas.

% Pitch: valor de cuefrecuencia al que se produce el primer pico del
% cepstrum para valores superiones a n = 20



%% 4. RECONOCIMIENTO DE VOCALES

% Una buena manera de representar la información relativa al tracto vocal
%(define el tipo de sonido que hemos emitido) es mediante un vector de parámetros que contenga los L primeros
%coeficientes cepstrales. c(0) no suele incluir en el vector ya que está
%relacionado con la energia de la señal, sometido a una alta variabilidad.
%El cepstrum LPC proporciona mejores resultados que el FFT, por lo que es
%el que usaremos

% a
a41=lpc(a,12); % Ya lleva el coeficiente a0=1 inclusive  calculado

suma41=0;
for k=1:1:12;
   suma41=suma41+a41(k+1)*Autoc(k,a);
end

b0_41=Autoc(0,a)+suma41;
PSDFormula41=b0_41./(abs(fft(a41,256))).^2;
c41=ifft(log(PSDFormula41));

% e
a42=lpc(e,12); % Ya lleva el coeficiente a0=1 inclusive  calculado

suma42=0;
for k=1:1:12;
   suma42=suma42+a42(k+1)*Autoc(k,e);
end

b0_42=Autoc(0,e)+suma42;
PSDFormula42=b0_42./(abs(fft(a42,256))).^2;
c42=ifft(log(PSDFormula42));

% % i
a43=lpc(i,12); % Ya lleva el coeficiente a0=1 inclusive  calculado

suma43=0;
for k=1:1:12;
   suma43=suma43+a43(k+1)*Autoc(k,i);
end

b0_43=Autoc(0,i)+suma43;
PSDFormula43=b0_43./(abs(fft(a43,256))).^2;
c43=ifft(log(PSDFormula43));

% o
a44=lpc(o,12); % Ya lleva el coeficiente a0=1 inclusive  calculado
suma44=0;

for k=1:1:12;
   suma44=suma44+a44(k+1)*Autoc(k,o);
end

b0_44=Autoc(0,o)+suma44;
PSDFormula44=b0_44./(abs(fft(a44,256))).^2;
c44=ifft(log(PSDFormula44));

% u
a45=lpc(u,12); % YEl segmento x es una aa lleva el coeficiente a0=1 inclusive  calculado
suma45=0;

for k=1:1:12;
   suma45=suma45+a45(k+1)*Autoc(k,u);
end

b0_45=Autoc(0,u)+suma45;
PSDFormula45=b0_45./(abs(fft(a45,256))).^2;
c45=ifft(log(PSDFormula45));

% x
aRf=lpc(x,12); % Ya lleva el coeficiente a0=1 inclusive  calculado
sumaRf=0;

for k=1:1:12;
   sumaRf=sumaRf+aRf(k+1)*Autoc(k,x);
end

b0_Rf=Autoc(0,x)+sumaRf;
PSDFormulaRf=b0_Rf./(abs(fft(aRf,256))).^2;
cRf=ifft(log(PSDFormulaRf));

%Calculamos la distancia euclidea entre el cepstrum y cada una de las
%vocales.

suma_a=0;

for k=1:1:12
    suma_a=suma_a+(c41(k)-cRf(k))^2;  
end

Resultado(1)=suma_a;

suma_e=0;

for k=1:1:12
    suma_e=suma_e+(c42(k)-cRf(k))^2;  
end

Resultado(2)=suma_e;

suma_i=0;

for k=1:1:12
    suma_i=suma_i+(c43(k)-cRf(k))^2;  
end

Resultado(3)=suma_i;

suma_o=0;

for k=1:1:12
    suma_o=suma_o+(c44(k)-cRf(k))^2;  
end

Resultado(4)=suma_o;

suma_u=0;

for k=1:1:12
    suma_u=suma_u+(c45(k)-cRf(k))^2;  
end

Resultado(5)=suma_u;

minimo=Resultado(1);
indice=0;

for k=1:1:5
    if(Resultado(k)<= minimo)
        minimo=Resultado(k);
        indice=k;
    end
end

if(indice==1)
    disp('El segmento x es una a');
end
if(indice==2)
    disp('El segmento x es una e');
end
if(indice==3)
    disp('El segmento x es una i');
end
if(indice==4)
    disp('El segmento x es una o');
end
if(indice==5)
    disp('El segmento x es una u');
end
