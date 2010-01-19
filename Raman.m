clc
clear all
close all

%% Se cargan los datos:
% 39 vinos con raman. 66 vinos con mediciones químicas
load Data.mat
load Data2.mat
clear Dark
clear Dark2

medQ=xlsread('Datos.xlsm'); %Mediciones químicas
datosImportados=importdata('Datos.xlsm');
nombresCarac=datosImportados.textdata(1,6:end);
cepas=datosImportados.textdata(1:end,4);
%% Promedio de todos los espectros que tengan el mismo ID

% Promedio de los espectros de Data.mat
espectroTemp=[];
aux=0;
for i=1:21
    rango=(30*(i-1)+1):30*i;
    espectroPromediado(i,:)=mean(Data(rango,5:end),1);
    ID(i)=Data(30*(i-1)+1,1);
end

for i=1:18
    rango=(9*(i-1)+4):(9*i);
    espectroPromediado(i+21,:)=mean(Data2(rango,5:end),1);
    ID(i+21)=Data2(9*i,1);    
end

colores='rgbcmykrgbcmykrgbcmykrgbcmykrgbcmykrgbcmykrgbcmykrgbcmykrgbcmykrgbcmykrgbcmykrgbcmykrgbcmyk';
figure,
hold on
for i=1:39
    plot(espectroPromediado(i,:),num2str(colores(i)))
end
xlabel('Nanometers')
ylabel('Intensity')
%% Ploteo fantabuloso
filtro=fspecial('gaussian', [1 500], 4);
colores='rgbcmykrgbcmykrgbcmykrgbcmykrgbcmykrgbcmykrgbcmykrgbcmykrgbcmykrgbcmykrgbcmykrgbcmykrgbcmyk';
figure,
hold on
for i=1:39
    ruidoso=espectroPromediado(i,:); 
    suave=imfilter(ruidoso,filtro,'symmetric','same');
    diferencia(i,:)=(ruidoso-suave)';
    plot(diferencia(i,:),num2str(colores(i)))
end
hold off

%% Se conservan sólo los rangos utilizables del espectro

diferencia=diferencia(:,[15:286 1100:2048]);


%% Normalización de la diferencia respecto al peak ubicado en la posición 1128         CHAO

diferenciaNormalizada=diferencia;

%% Se elimina la medición 1002 (porque no tiene raman UPDATE:¿NO TIENE MEDICIÓN QUÍMICA?)

ID=ID([1 3:39]);
diferenciaNormalizada=diferenciaNormalizada([1 3:39],:);
diferencia=diferencia([1 3:39],:);

%% Se ordenan las mediciones químicas para que todos los índices calcen

aux=1;
for i=1:38    
    lugar=find(medQ(:,1)==ID(i));    
    if(~isempty(lugar))
        medQOrdenado(aux,:)=medQ(lugar,:);
        cepasOrdenadas(aux,:)=cepas(lugar,:);
        aux=aux+1;
    end    
end

%%
% close all
CV=std(diferenciaNormalizada,1)./mean(diferenciaNormalizada,1);
figure, plot(CV)
hold on
plot(diferenciaNormalizada(1,:),'r')
hold off

CVaceptables=find(abs(CV)<1);

%%
diferenciaNormalizadaSegura=diferenciaNormalizada(:,CVaceptables);

%Se busca la posición con menor CV
temp=std(diferenciaNormalizadaSegura,1)./mean(diferenciaNormalizadaSegura,1);
[minimo indiceMinimo]=min(abs(temp))


%%
figure,
hold on
for i=1:38
    plot(diferenciaNormalizadaSegura(i,:),num2str(colores(i)))    
end
hold off
%% Se normaliza el espectro que va quedando respecto a la posición con menor CV
figure,
hold on
for i=1:38
    plot(diferenciaNormalizadaSegura(i,:)/diferenciaNormalizadaSegura(i,indiceMinimo),num2str(colores(i)))
    diferenciaNormalizadaSegura(i,:)=diferenciaNormalizadaSegura(i,:)/diferenciaNormalizadaSegura(i,indiceMinimo);
end
hold off

%% Mientras tanto... clasificación de cepa ES O NO ES CABERNET SAUVIGNON

for i=1:length(cepasOrdenadas)
    cepas2num(i)=1;
    coincidencia=strcmp(cepasOrdenadas(i,:),'Cabernet Sauvignon');
    if(coincidencia==1)
        cepas2num(i)=2;
    end    
end
cepas2num=cepas2num';
X=diferenciaNormalizadaSegura;

%%
[J,selec] = sfs(X,cepas2num,2,1)
% [J,selec] = exsearch(X(:,1:end),cepas2num,2,1)
%%
figure, plotfeatures(X(:,selec),cepas2num,{'Cabernet Sauvignon'; 'Merlot'; 'Carmenere'; 'Syrah'; 'Cepaje'})

%%
[C,err,P,logp,coeff] = classify(X(:,selec),X(:,selec),cepas2num,'linear');

K = coeff(1,2).const;
L = coeff(1,2).linear; 
f = sprintf('0 = %g+%g*x+%g*y', K,L);
h2 = ezplot(f,[-1 1 -1 1]);
set(h2,'Color','m','LineWidth',2)