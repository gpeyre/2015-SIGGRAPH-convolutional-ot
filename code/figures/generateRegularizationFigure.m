%% Set up problem

n = 500;

x = linspace(0,1,n)';

p1 = .9*normpdf(x,.5,.1)+.1*normpdf(x,.1,.05);
p1 = p1/sum(p1)*n;
p2 = .75*normpdf(x,.25,.05)+.25*normpdf(x,.75,.1);
p2 = p2/sum(p2)*n;

%% Solve ground truth

costs = zeros(n,n);
for i=1:n
    for j=1:n
        costs(i,j) = (x(i)-x(j))^2;
    end
end

cvx_begin
    variable T(n,n)
    cvx_solver gurobi

    minimize sum(sum(T.*costs))
    
    subject to
        T >= 0;
        sum(T,1)'/n == p1;
        sum(T,2)/n == p2;
cvx_end

%% Write transportation 

TT = T - min(T(:));
TT = TT / max(TT(:));
TT = 1 - TT;

imwrite(TT,'no_regularization.png','png');

%% Solve regularized problems

gammas = [.0001 .001 .01 .1];

for experiment=1:length(gammas)
    TR = zeros(n,n);
    gamma = gammas(experiment);
    
    for i=1:n
        for j=1:n
            TR(i,j) = exp(-(x(i)-x(j))^2/gamma);
        end
    end
    
    for i=1:5000 % implemented the simple version of Sinkhorn
        rowSums = sum(TR,2);
        TR = bsxfun(@times,TR,p1./rowSums);
        colSums = sum(TR,1);
        TR = bsxfun(@times,TR,p2'./colSums);

        if mod(i,100)==0
            fprintf('%d\n',i);
        end
    end
    
    TR = TR - min(TR(:));
    TR = TR / max(TR(:));
    TR = 1 - TR;
    
    imwrite(TR,sprintf('gamma_%g.png',gamma),'png');
end

%%

plot(x,p1,'k');
axis off
saveas(gcf,'p1.pdf')

plot(x,p2,'k');
axis off
saveas(gcf,'p2.pdf')