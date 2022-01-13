N = 100000;
L = 8;

figure;
g = multipath(2,N);
cdf_plotting(abs(g),'CDF of channel magnitude')
hold on;
g = multipath(8,N);
cdf_plotting(abs(g),'CDF of channel magnitude')
g = sqrt(1/2)*(randn(1,N)+1i*randn(1,N));
cdf_plotting(abs(g),'CDF of channel magnitude')
legend('Two paths CDF', 'Eight paths CDF', 'Rayleigh CDF')

figure;
histogram(log(abs(g)),Normalization="pdf")
% A function that returns random channel when number of multipaths are
% specified
% When L -> Infinite the distribution closes in on rayleigh fading
function channel = multipath(L,N)
    alpha = 1/sqrt(L);
    theta = rand(L,N)*2*pi;
    channel = alpha*(sum(exp(-1i*theta)));
end

function cdf_plotting(X,title_)
    tmp = sort(reshape(X,numel(X),1));
    Xplot = reshape([tmp tmp].',2*length(tmp),1);
    
    tmp = [1:length(X)].'/length(X);
    Yplot = reshape([tmp tmp].',2*length(tmp),1);
    Yplot = [0; Yplot(1:(end-1))];
    
    figure(gcf);
    hp = plot(Xplot, Yplot);
    
    ColOrd = get(gca, 'ColorOrder'); 
    ord = mod(length(get(gca,'Children')), size(ColOrd,1)); 
    set(hp, 'Color', ColOrd((ord==0) + (ord>0)*ord, :));
    if ~ishold
         xlabel('X', 'FontWeight','b','FontSize',12);
         ylabel('F(X)', 'FontWeight','b','FontSize',12);
         title(title_, 'FontWeight','b','FontSize',12);
         grid on;
    end
end