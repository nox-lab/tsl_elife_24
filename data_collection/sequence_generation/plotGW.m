function plotGW(outcome, para, fig)

% Plotting parameters
nr = 1; % no. rows of subplots
nc = 3; % no. columns of subplots
subplots = 1:2; % subplot numbering 
colstrs = {'Small','Large'};
glbl = {'Small','Large'};
lgtitle = 'Stochasticity';
vol = {'Low volatility', 'high volalitily'};

ii = [1 3;2 4];
I = para(5);
bound = (I - 1)*0.5;

col = {'red', 'blue'};
fn = 8;
fsy = 8; %axes fontsize
fsl = fsy+2;      %title fontsize
yl = [-(bound +5) bound+5];  %y axes range

for i=1:2
    figure(fig)
    h(i) = subplot(nc,nr,subplots(i));
    hp = nan(1,2);    
    hx = nan(1,2);
    % plot(state(:,i),'-','linewidth',2,'color','k');

    hold on;
    title(vol(i));
    ylb = yline(bound, '--', 'Minimum intensity level'); ylb.Color = 'red';
    ylt = yline(-bound, '--', 'Maximum intensity level'); ylt.Color = 'red';
    for j=1:2
        ij = ii(i,j);
        hp(j) = plot(round(outcome(:,ij)),'-','color', col{j},'linewidth',1); hold on;
%         hx(j) = plot(outcome(:,ij),'.','color','yellow','markersize',10); hold on;
    end    
    set(gca,'box','off');
      
    if i==1
        lg = legend(hp,glbl,'orientation','horizontal');    
        title(lg,lgtitle);
    end
    ylim(yl);        
    ylabel('Outcome','fontsize',fsy);
    if i==2
        xlabel('Trial','fontsize',fsy);
    end
    
    hold off;
end

line1 = ['I = ', num2str(I), ',LV = ', num2str(para(1)), ', HV = ', num2str(para(2)), ', LS = ', num2str(para(3)), ', HS = ', num2str(para(4))];
sgtitle(line1)

end
