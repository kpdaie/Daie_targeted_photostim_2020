% initialize matlab figures
marg = .4;
figure(1)
figure_initialize
set(gcf,'position',[4 4 2 2]);
figure(2)
figure_initialize
set(gcf,'position',[4 2 2 1]);
% order = [0 3 2 1];
order = [0 1];
marg = [.2 .2];
clf

ei_ratio = 0;
amp      = .8;
for monolith = [1 0];
    figure(1);
    N = 1000;
    positions = rand(N,1)*1000;
    dist = squareform(pdist(positions));
    NeuronsPerCluster = 50;
    frac = .5;
    R = rand(N,1)>frac;
    L = ~R;
    numClusters = round(N*(1-frac)/NeuronsPerCluster);
    clustCenters = linspace(0,1000,ceil(sqrt(numClusters)));
    [a,b] = meshgrid(clustCenters,clustCenters);
    clustCenters = [a(:) b(:)];
    x = linspace(0,1,100);
    prob = @(x) exp(-x/200)+.0*exp(-x/1000);
    clustID = zeros(N,1);
    for ci = 1:numClusters
        stp = 0;
        while stp == 0
            d = sqrt(sum((repmat(clustCenters(ci,:),N,1) - positions).^2,2));
            stp = length(d)>(NeuronsPerCluster-10);
            ind = find(clustID == 0);
            cells = randsample(ind,NeuronsPerCluster,1,prob(d(ind)));
            cells = unique(cells);
        end
        clustID(cells) = ci;
        len(ci) = length(cells);
    end
    R = clustID ~= 0;
    L = clustID == 0;
    inh = @(x) exp(-x./200).*(1-exp(-x./200));
    ex = @(x) exp(-x/70);
    ex_clust = @(x) exp(-x/500);
    EI = @(x) 3*ex(x) - inh(x);
    g = EI(dist);g = g - eye(size(g));
    g = g.*(rand(N)>.8);
    g = g/eigs(g,1)*.7;
    yes = double(repmat(clustID,1,N)==repmat(clustID',N,1)).*(clustID*clustID'~=0);
    no = ~yes;
    w = zeros(N);
    
    for i = 1:numClusters;
        pop = find(clustID == i);
        %     v = zeros(1,N);
        %     v(pop) = rand(1,length(pop));
        %     v = v/norm(v);
        % %     w = w + v'*v;
        a = zeros(N);
        d = squareform(pdist(positions(pop,:)));
        d = ex_clust(d);
        d = d/max(eig(d));
        a(pop,pop) = d;
        w = w + a;
    end
    
    if monolith == 1
        ex_clust = @(x) exp(-x/200);
        w = zeros(N);
        a = ex_clust(dist);
        w(R==1,R==1) = a(R==1,R==1);
        w = w .* (rand(N)<.1);
        w = w/max(eig(w));
    end
    
    if monolith ~= 1
        %         if monolith == 0
        %             ei_ratio = 0;amp = 1;
        %         elseif monolith == 2
        %             ei_ratio = 1;amp = 1;
        %         elseif monolith == 3;
        %             ei_ratio = 0;amp = 0;
        %         end
        w_inter = (rand(N)-ei_ratio).*(rand(N)<.1).*no;
        w_inter(R==1,L==1) = 0;w_inter(L==1,R==1) = 0;
        w_inter(L==1,L==1) = 0;
        w_inter = w_inter/abs(real(eigs(w_inter,1)))*.3*amp;
        g(R==1,:) = 0;
        g(:,R==1) = 0;
        w_mod = w - w_inter;
    elseif monolith == 1
        w_mod = w;
    end
    w = w_mod+g;
    wRL = (randn(N).*no.*(rand(N)<.05))*.01;
    wRL(R==1,R==1) = 0;
    wRL(L==1,L==1) = 0;
    w = w + wRL;
    w = w - eye(size(w)).*w;
    w = w/max(real(eigs(w,5)));
    if monolith == 0
        [a,b] = sort(clustID,'descend');
    else
        [a,b] = sort(positions);
    end
    if monolith < 2
        KDsubplot(2,2,[1 monolith+1],marg);
        imagesc(w(b,b),[-1 1]/40)
        gap = ones(1,4);
        EIratio = 1;
        short = floor(100/2/EIratio);long = EIratio*short;
        cm = [linspace(0,1,short) gap exp(-linspace(0,1,long)/29);...
            linspace(0,1,short) gap exp(-linspace(0,1,long)/1.3);...
            %     linspace(0,1,20) linspace(1,0,40);
            linspace(0,1,short) gap exp(-linspace(0,1,long)/.2);];
        colormap(cm');
        set(gca,'xtick',[1 1000],'ytick',[1 1000]);
        if monolith == 1
            set(gca,'ytick',[],'xtick',[1 1000]);
        end
    end
    % subplot(212)
    % mean_bin_plot(dist(:),w_mod(:),55,1,1,'b')
    % mean_bin_plot(dist(:),g(:),55,1,1,'k')
    % plot(xlim,xlim*0,'k:');
    % set(gca,'xtick',[0 500 1000]);
    % set(gca,'ytick',[0 4e-3]);
    
    % figure(2)
    % clf
    dt = .01;
    t = 0:dt:2.4;
    T = length(t);
    fun = @(x) x;
    tau = .1;
    sel = ((rand(N,1)-0).*(R==1)) - ((rand(N,1)-0).*(R==0));
    sel = sel.*(rand(size(sel))>-.5);
    clear selR dbr
    iters = 50;
    I = 0;
    for k = [2 1]
        for iter = 1:iters + 2;
            iter
            x = 1:N;            
            sampR = [t>.1 & t<.2]'*double(R'==1);
            sampL = [t>.1]'*double(R'==0)*.1;
            samp = sampR - sampL;
            i_t = 0*t;i_t(40) = 1;
            first = 100;
            p_t = 0*t;p_t(first:first+20) = 1;
            r = zeros(T,N);
            pert = [0 0 1];
            sound = [1 -1 1];
            i_s = zeros(1,N);
            for jj = 1:3
                if jj == 3
                    i_s = zeros(1,N);
                    if k == 1
                        clustIDs = randsample(1:numClusters,1);
                        cells = randsample(find(clustID==clustIDs & sel>=0),10);
                        if iter == iters+1
                            cells = find(clustID == 1);
                        elseif iter == iters + 2
                            cells = find(clustID == max(clustID));
                        end
                    elseif k == 2
                        cells = randsample(find(clustID==0 & sel<=0),10);
                    end
                    i_s(cells) = 100;
                    win = exp(-x.^2/10^2);
                    %         i_s = conv(i_s,win);
                    i_s = i_s(1:N);
                end
                for i = 1:(T-1)
                    r(i+1,:) = r(i,:) + dt/tau*(-r(i,:) + fun(r(i,:))*w + i_t(i)*i_s*pert(jj) + samp(i,:)*sound(jj));
                end
                rrr(:,:,jj) = r;
                if jj == 2
                    sel = squeeze(sum(rrr));
                    sel = sel(:,1) - sel(:,2);
%                     sel = squeeze(mean(rrr(end-150:end,:,1) - rrr(end-150:end,:,2)))';
                end
            end
            dr = rrr(:,:,3) - rrr(:,:,1);                        
            v = sel>0;
            v = sel;
            readout = dr*v;
            readout = readout./((rrr(:,:,1) -  rrr(:,:,2))*v);
            dbr(:,iter,k) = readout;
            if k == 1 & iter>iters
                I = I + 1;
                rr(:,:,I) = dr;
            end
            selR(iter,k) = mean(sel(cells));
            %             selR(iter,k) = mean(sel(cells));
        end
    end
    if size(selR,1)>1
        figure(2)
        KDsubplot(1,2,[1 order(monolith+1)+1],marg);
        [a,B,c,d,e] = mean_bin_plot(100*squeeze(dbr(end,:,:)),selR/100,5,0,1,'k')
        %         a = a/max(abs(a))*10;
        KD_errorbar(B,a,c,'k','o-',3,1);hold on;
        if monolith == 1
            xlim([-1 4]);ylim([-5 25]);
            set(gca,'xtick',[0 4],'ytick',[0 20]);
        else
            xlim([-1 2]);ylim([-11 11]);
            set(gca,'ytick',[-10 0 10],'xtick',[0 1]);
        end
%         a = ylim;
%         if a(1) >= 0;
%             ylim([-5 a(2)]);
%         end
        plot(xlim,xlim*0,'k:');
        plot(ylim*0,ylim,'k:');
    
        %         title(num2str(monolith));
    end
    
    if monolith < 2
        figure(1)
        KDsubplot(2,2,[2 monolith+1],marg);
        a = zeros(N,T,3);
        a(:,:,1) = rr(:,b,1)';
        a(:,:,3) = rr(:,b,2)';
        image(a);
        % imagesc(rr(:,b,1)',[-1 1]/3);
        set(gca,'xtick',[40 140 240],'xticklabel',{'0','1','2'},'ytick',[1 1000]);
        if monolith == 1
%             set(gca,'ytick',[]);
        end        
        colormap(cm')
    end
end
figure(1)
figure_finalize
figure(2)
figure_finalize