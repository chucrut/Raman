% J = jfisher(V,c,p)
%
% Toolbox: Balu
%    Fisher objective function J.
%    V features matrix. V(i,j) is the feature j of sample i.
%    c vector that indicates the ideal classification of the samples
%    p a priori probability of each class
%
% D.Mery, PUC-DCC, Apr. 2008
% http://dmery.ing.puc.cl


function J = jfisher(V,c,p)

% M: numero de caracteristicas
% N: numero de clases
% n: numero de muestras por clase

[n,M] = size(V);
cmin = min(c);
cmax = max(c);

N = cmax;

if not(exist('p'))
    p = ones(N,1)/N;
end

% la media de las caracteristicas es Vm
Vm = mean(V)';

L  = zeros(N,1);
Cw = zeros(M,M);
Cb = zeros(M,M);

for k=1:N
    ii   = find(c==k); % filas que pertenecen a la clase k
    L(k) = length(ii); % numero de muestras de la clase k
    Vk   = V(ii,:);    % muestras de la clase k
    Vkm  = mean(Vk)';  % media de la clase k 
    Ck = cov(Vk);      % covarianza de la clase k
    
    % calculo de la covarianza interclase (within-class covariance)
    Cw = Cw + p(k)*Ck; 
    
    % calculo de la covarianza intraclase (between-class covariance)
    Cb = Cb + p(k)*(Vkm-Vm)*(Vkm-Vm)';
end

% calculo del discriminate Fisher
J = trace(inv(Cw)*Cb);

