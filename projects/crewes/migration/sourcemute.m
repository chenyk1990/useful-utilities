function shotout=sourcemute(shotin,z,x,xshot,zshot,rad)

x=x(:)';
z=z(:);

xx=x(ones(length(z),1),:);
zz=z(:,ones(1,length(x)));
r2=(xx-xshot).^2+(zz-zshot.^2);
mask=1-exp(-(r2./rad^2));

shotout=shotin.*mask;