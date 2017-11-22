function [err,cc]=waveleterr(w,tw,w1,tw1)

tbeg=tw(1);
tend=tw(end);
tbeg1=tw1(1);
tend1=tw1(end);
dt=tw(2)-tw(1);
if(tbeg<tbeg1)
    %pad w1 in the beginning
    npad=round((tbeg1-tbeg)/dt);
    tw1=[(tbeg:dt:tbeg1-dt)';tw1];
    w1=[zeros(npad,1);w1];
end
if(tbeg>tbeg1)
    %pad w in beginning
    npad=round((tbeg-tbeg1)/dt);
    tw=[(tbeg1:dt:tbeg-dt)';tw];
    w=[zeros(npad,1);w];
end
if(tend>tend1)
    %paw w1 at end
    npad=round((tend-tend1)/dt);
    tw1=[tw1;(tend1+dt:dt:tend)'];
    w1=[w1;zeros(npad,1)];
end
if(tend<tend1)
    %pad w at end
    npad=round((tend1-tend)/dt);
    tw=[tw;(tend+dt:dt:tend)'];
    w=[w;zeros(npad,1)];
end

err=norm(w-w1);
cc=maxcorr(w,w1);