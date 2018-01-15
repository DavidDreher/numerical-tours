options.ti = 1;
n = 256;
Jmin = 2;
y = nm-toolbox_signal.perform_wavelet_transf(zeros(n,n), Jmin, +1, options);
nJ = 3;
nT = 3;
sel = linspace(0,n-1,nT+1); sel = sel(1:nT); sel = sel+round(sel(2)/2);
clf;
k = 0;
for j=1:nJ
    subplot(nJ,1,j);
    x = [];
    for i=1:nT
        k = k+1;
        y = y*0;
        y(end/2,end/2,3*(j-1)+i+1) = 1;
        x = nm-toolbox_signal.perform_wavelet_transf(y, Jmin, -1, options);
        subplot(nJ,nT,k);
        nm-toolbox_signal.imageplot(x);
    end
end
