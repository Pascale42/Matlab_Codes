% 
% function GridNeuCM(uClus,uPre, uPost)


function GridNeuCM(uClus,uPre, uPost)


% lines for layers
uLy=unique(uClus(:,3));
    % in pre
    text(1,1,num2str(uLy(1)),'FontSize',12,'FontWeight','bold') % note layer pre
    for n= 2:length(uLy)
    xpr=find(uClus(:,3) == uLy(n),1, 'first'); 
    line([(xpr -0.5) (xpr-0.5)],[0.5, size(uClus,1)+0.5])
    text(xpr,1,num2str(uLy(n)),'FontSize',12,'FontWeight','bold') % note layer pre
    end
    clear  xpr
    % in post
    text(1,1.5,num2str(uLy(1)),'FontSize',12,'FontWeight','bold') % note layer pre
    for n= 2:length(uLy)
    xpo=find(uClus(:,3) == uLy(n),1, 'first');
    line([0.5,size(uClus,1)+0.5],[(xpo -0.5) (xpo-0.5)])
    text(1.5,xpo,num2str(uLy(n)),'FontSize',12,'FontWeight','bold') % note layer post
    end
    clear uLy xpr xpo
% lines for 
        % Pre neurons
   mb=ismember(uClus(:,1), uPre(:,1));
x=1: length(uClus);
for n= x(mb)
    line([n n],[0.5,size(uClus,1)+0.5],'color','k','LineStyle',':')
end
 % Post neurons
   mb=ismember(uClus(:,1), uPost(:,1));
x=1: length(uClus);
for n= x(mb)
    line([0.5,size(uClus,1)+0.5],[n n],'color','k','LineStyle',':')
end
   clear mb x
