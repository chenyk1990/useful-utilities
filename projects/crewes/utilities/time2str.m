function str=time2str(time)
%
% str=time2str(time)
%
% convert a time to a string with the deicmals to the nearest milliseconds
%
atime=abs(time);
if(atime<1 && atime>.001)
    str=num2str(time,3);
elseif(atime>1 && atime<10)
    str=num2str(time,4);
elseif(atime>10)
    str=num2str(time,5);
elseif(time<.001 && time>.000499999)
    str=num2str(sign(time)*.001);
else
    str=num2str(0);
end
