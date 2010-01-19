% plotfeatures(z,c,featurename)
%
% Toolbox: Balu
%    Plot features z acording classification c. If the feature names are given
%    then they will labeled in each axis.
%
% D.Mery, PUC-DCC, Apr. 2008
% http://dmery.ing.puc.cl
%

function plotfeatures(z,c,nombres)

m = size(z,2);


cmin = min(c);
cmax = max(c);

n = cmax-cmin+1; % number of classes
col = 'gbrcmykbgrcmykbgrcmykbgrcmykgbrcmykbgrcmykbgrcmykbgrcmykgbrcmykbgrcmykbgrcmykbgrcmyk';
mar = 'ox+v^ox+v^ox+v^ox+v^ox+v^ox+v^ox+v^ox+v^ox+v^ox+v^ox+v^ox+v^ox+v^ox+v^ox+v^ox+v^ox+v^ox+v^';
clf
warning off
r = sprintf('%s',39);
warning on
s = 'legend(';
for k = cmin:cmax
    cell2mat(nombres(k,:))
    s = [s r sprintf('%s',cell2mat(nombres(k,:))) r ];
    if k<cmax
        s = [s ','];
    else
        s = [s ');'];
    end
end
if (m<4)
    switch m
        case 1
            for k = cmin:cmax
                ii = find(c==k);
                [h,x] = hist(z(ii));
                dx = x(2)-x(1);
                A = sum(h)*dx;
                x = [x(1)-dx x x(end)+dx];
                h = [0 h 0];
                plot(x,h/A,col(k+1));
                hold on
            end
            if (exist('featurename'))
                xlabel(featurename);
            else
                xlabel('feature value');
            end
            ylabel('PDF');
        case 2
            for k = cmin:cmax
                ii = find(c==k);
                plot(z(ii,1),z(ii,2),[col(k+1) mar(k+1)]);
                hold on
            end
            if (exist('featurename'))
                xlabel(featurename(1,:));
                ylabel(featurename(2,:));
            else
                xlabel('best feature 1');
                ylabel('best feature 2');
            end
            title('feature space');
        case 3
            for k = cmin:cmax
                ii = find(c==k);
                plot3(z(ii,1),z(ii,2),z(ii,3),[col(k+1) mar(k+1)]);
                hold on
            end
            if (exist('featurename'))
                xlabel(featurename(1,:));
                ylabel(featurename(2,:));
                zlabel(featurename(3,:));
            else
                xlabel('feature value 1');
                ylabel('feature value 2');
                zlabel('feature value 3');
            end
    end
    eval(s)
else
    % disp('I cannot plot it for more than three features :(');
    l = 1;
    for j=1:m
        for i=1:m
           zi = z(:,i);
            zj = z(:,j);
            subplot(m,m,l); l = l+1;

            for k = cmin:cmax
                ii = find(c==k);
                plot(zi(ii),zj(ii),[col(k+1) mar(k+1)]);
                hold on
            end
            if (exist('featurename'))
                xl = featurename(i,:);
                yl = featurename(j,:);
            else
                xl = sprintf('z_%d',i);
                yl = sprintf('z_%d',j);
            end
            if i==1
                ylabel(yl)
            end
            if j==m
                xlabel(xl)
            end
        end
    end


end