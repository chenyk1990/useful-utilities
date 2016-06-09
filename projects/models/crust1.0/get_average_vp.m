% code to calculate the average crustal Vp
% Yunfeng Chen, Oct. 15, 2014
function [ave_vp,vp,thick] = get_average_vp(lat,lon)
% function [ave_vp,vp,thick] = get_avarage_vp(lat,lon)
% Input: lat, lon
% Output: ave_vp -> average crustal vp
%         vp -> vp for each layer
%         thick -> thickness for each layer

% get the structure info for certain location
opt = 'rocks';
[out, ~, ~] = get_crust1_model(lat,lon,opt);
vp = out(:,1);
thick = out(:,end);
ave_vp = sum(thick)/sum(thick./vp);
end