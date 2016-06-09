% code to read in the model results from crust 1.0 model
% Yunfeng Chen, Oct. 14th, 2014
function [out, elevation,moho] = get_crust1_model(lat,lon,opt)
% This function will grab the crustal model and format it as needed.
  % This executes a command that will output a model of the crust near the specified point.
  command=sprintf('./getCN1point_modified << EOC\n%f %f\nEOC\n',lat, lon);
  [~,out]=unix(command);
  out=str2num(out); %#ok<ST2NM>
  
  % Notes on the output of getCN1point_modified:
  % first row has the elevation and then just zeros.
  % A table of values is output columns are: Vp(km/s), Vs(km/s), density(Mg/m^3), bottom of layer (km), 
  % The rows are: elevation, water, ice, soft sed, hard sed, upper crust, middle crust, lower crust.
  elevation=out(1,1);
  
%   % Reversing order so that ice layer is on top of water.
%   out(1,:)=out(3,:);
%   out(3,:)=zeros(1,4);
  
  % Check to see what layers should be kept.
  if(strcmp(opt,'rocks'))
      % User has requested only the rock layers.
      out(1:3,:)=zeros(3,4);
  elseif(strcmp(opt,'all'))
      % All layers.  Reversing order so that ice layer is on top of water.
      out(1,:)=out(3,:);
      out(3,:)=zeros(1,4);
  elseif(strcmp(opt,'none'))
      % Just returning the thickness of the crust
      out=sum(out(4:end,end));
  else
      out=NaN;
  end;
      
  %Remove any layers of zero vp.
  out=out(out(:,1)>0,:);
  
  % Change from depth to thickness.
  out(:, end+1) = out(:, end);
  % find the layer above the see level
  keep = out(:,end) > 0;
  out(keep,end) =  -diff([out(keep,end);0]);
  % find the layer below the see level
  keep = out(:,end) < 0;
  out(keep,end) =  -diff([0;out(keep,end)]);
  % Moho depth
  moho = -out(end,end-1);
return;
