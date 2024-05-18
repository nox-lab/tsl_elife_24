% plotting jump sequence 

function plotJS(Seq_store, para, f1)

d = para(1); I = para(2); Hs = para(3); Ls = para(4); Hv = para(5);
Lv = para(6); mu = para(7); sess = para(8); HVl = para(9); LVl = para(10);
jit = para(11);

bound = (I - 1)*0.5;

glbl = {'Small','Large'};
lgtitle = 'Stochasticity';
hp = nan(1,2); 

minmax_outcome = [-(bound+5):bound+5];

ind = 1;
        
for sub = 1:2:3
        figure(f1);        
        subplot(3,1,ind);
        hold on;
        %scatter(1:sess,s);
        % hp(j) = plot(round(outcome(:,ij)),'-','color', col{j},'linewidth',1);
        hp(1) = plot(1:sess,Seq_store{sub}, 'color', 'red','linewidth',1);
        hp(2) = plot(1:sess,Seq_store{sub+1}, 'color', 'blue','linewidth',1);
        xlabel('trial')
        ylabel('outcome');
        ylb = yline(bound, '--', 'Maximum intensity level'); ylb.Color = 'red';
        ylt = yline(-bound, '--', 'Minimum intensity level'); ylt.Color = 'red';
        
        ylim([min(minmax_outcome) max(minmax_outcome)])
        
        %ylim([min(possmu),max(possmu)]);
        if sub==1
            title('High volatility');
            lg = legend(hp,glbl,'orientation','horizontal');    
            title(lg,lgtitle);
        elseif sub==3
            title('Low volatility');
            
        end
        
        hold off 

        ind = ind + 1;
        

end

%line1 = ['I = ', num2str(I), ', LV = ',num2str(Lv),', HV = ',num2str(Hv),', LS = ',num2str(Ls),', HS = ',num2str(Hs)];
%line2 = [', Mu = ', num2str(mu), ', D = ', num2str(d)];
line1 = ['LV length = ',num2str(LVl),', HV length = ',num2str(HVl),', jitter = ',num2str(jit),', LS = ',num2str(Ls),', HS = ',num2str(Hs)];
%sgtitle([line1, line2])  
sgtitle(line1) 

end
