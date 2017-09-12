% -------------------- Visualize coarse mesh --------------------
% clear_onoff = 1 if the figure has to be cleared in advance
% hold_onoff  = 1 if 'hold on' has to be issued afterwards
% num_onoff   = 1 if the nodes should be numbered
%             = 2 if the edges should ne numbered             
%             = 3 if the faces should be numbered
% which_gplot = 1 if gplot_prive is to be called
% ----------------------------------------------------------------

function Bvisual_mesh(Node,Edge,Face,...
               which_gplot,clear_onoff,hold_onoff,num_onoff,fontsize)

if clear_onoff==1, clf, end

np    = 3;
nall  = size(Node,2);
nedge = size(Edge,2);
nface = size(Face,2);
AJ    = spalloc(nall,nall,5*nall);
for k=1:nall,  AC(k,:)=[Node(1,k) Node(2,k)]; end
for iedge=1:nedge
    AJ(Edge(1,iedge),Edge(2,iedge))=1;
    AJ(Edge(2,iedge),Edge(1,iedge))=1;
end

if which_gplot==1,
   gplot_prive(AJ,AC)
else
   gplot(AJ,AC)
end

if num_onoff>=1,
   hold on
   for k=1:nall,
%text(AC(k,1),AC(k,2),int2str(k),'Fontsize',fontsize)
text(AC(k,1),AC(k,2),int2str(k),...
'vertical','top','horizontal','right','fontname','times',...
'Fontsize',fontsize,'color','red');
   end
   axis('off')
end

if num_onoff>=2,
   hold on
   for k=1:nedge,
     mid = (Node(:,Edge(1,k))+Node(:,Edge(2,k)))/2;
     text(mid(1),mid(2),int2str(k),...
'vertical','top','horizontal','right','fontname','times',...
'Fontsize',fontsize-2,'color','black');
   end
   axis('off')
end

if num_onoff>=3,
   hold on
   for k=1:nface,
     s = 0;
     for l=1:size(Face(:,1),1)
         s = s + 1;
         cx(s) = Node(1,Edge(1,Face(l,k)));
         cy(s) = Node(2,Edge(2,Face(l,k)));
     end
     mid(1) = (min(cx)+max(cx))/2;
     mid(2) = (min(cy)+max(cy))/2;
     text(mid(1),mid(2),int2str(k),...
'vertical','middle','horizontal','center','fontname','times',...
'Fontsize',fontsize-2,'color','magenta');
   end
   axis('off')
end

if hold_onoff~=1, hold off,end
return
