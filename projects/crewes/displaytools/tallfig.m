pos=get(gcf,'position');
ssz=get(0,'screensize');
if(pos(4)>ssz(3))
    newid=ssz(3)-100;
else
    newwid=pos(4);
end
if(pos(3)>ssz(4))
    newht=ssz(4)-100;
else
    newht=pos(3);
end
set(gcf,'position',[pos(1:2) newwid newht])