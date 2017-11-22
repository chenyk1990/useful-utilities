function [mfilt,tm]=matchT(trin,trdsign,t,mlength,icausal,mu)
% MATCHS: A least-squares, smoothness and time constrained, matching filter, good for wavelet estimation. 
%
% [mfilt,tm]=matchT(trin,trdsign,t,mlength,icausal,mu)
%
% MATCHT designs a match filter of specified temporal length which matches trin to trdsign in
% the least squares sense and is constrained to be smooth and time localized. The smoothness
% constraint is imposed by requiring the second derivative of the filter (the model) to be
% small at the same time as we minimize the error of the data fit. The time constraint is that
% wavelet samples far from time zero are penalized. This is done by applying a window=1-g,where
% g is a Gaussian window centered at time zero, to the solution and requiring that the windowed
% answer be small at the same time as the other two terms.  Thus the total objective function
% is A+mu(1)*B+mu(2)*C where A is the sum-squared error, B is the energy (sum of squares) of
% the second derivative, and C is the energy of the windows signal. The objective function is
% minimized and mu(1) and mu(2) provide the weights for the two constraints. The mus values
% should be non-negative. The filter can be causal or non causal. See Constable et al, 1987,
% Geophysics, "Occams's inversion: a practical algorithm for generating smooth models from
% electromagnetic sounding data".
%
% Details on the time constraint: let g be a gaussian window centered at time zero and of the
% same length as the match filter. The standard deviation, sigma, of the gaussian is set such
% that it is 2*sigma down at the ends of the match filter. Thus if the match filter is time
% symmetric (i.e. t=0 is in the middle) then g is symmetric, otherwise, g is not. Then, the
% window W=1-g, when applied to the match filter, emphasizes those samples near the ends while
% the samples near t=0 get nearly zero weight. Thus requiring that W.*mfilt be small forces
% most of the filter energy to be within +/- one sigma of the origin. This works well if the
% misalignment between trin and trdsign is small but otherwise can be erroneous. The choise of
% a gaussian window and 2 sigma down at the ends is arbitrary. Ths can be changed by modifying
% the internal function designtimeweights at the end of this file.
%
% trin= input trace to be matched to trdsign
% trdsign= input trace which is to be matched
% t= time coordinate vector for trin
% ***** note: trin and trdsign must be the same length
% mlength= length of the match filter in seconds
% icausal=0 ... a fully noncausal operator is desired. (t=0 is in the middle)
%     =1 ... a causal operator is desired (t=0 is the first sample).
%     =N ... where N>1 means that the operator will have N-1 samples in negative time. This
%     allows a continuum of variation between symmetric and causal. (t=0 will be sample N).
% NOTE: icausal=0 is the same as icausal=nsamp/2+1 where nsamp is the closest odd number of samples
%   to the requested filter length. This control is about how many samples before a
%   given time and after that time are needed to predict trdesign. Fully causal operators, even
%   if the wavelet is known to be causal, are generally not as good as those where icausal is
%   slightly greater than 1 (say 2 through 10).
% mu ... tradeoff parameters controlling the tradeoff between filter (model)
%       smoothness and time constraint and data fitting. Should be two nonnegative numberw. The
%       first is the smoothness weight and the second is the time constraint weight. Larger
%       values give stronger constraints while [0 0] gives the best fit to the data with no constraints.
% ************* default mu=[1 1] **********
%
% mfilt= output mlength match filter
% tm= time coordinate for the match filter
%
% by G.F. Margrave, 2016
%
% NOTE: This SOFTWARE may be used by any individual or corporation for any purpose
% with the exception of re-selling or re-distributing the SOFTWARE.
% By using this software, you are agreeing to the terms detailed in this software's
% Matlab source file.

% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by
% its 'AUTHOR' (identified above) and the CREWES Project.  The CREWES
% project may be contacted via email at:  crewesinfo@crewes.org
%
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) This SOFTWARE may be used by any individual or corporation for any purpose
%    with the exception of re-selling or re-distributing the SOFTWARE.
%
% 2) The AUTHOR and CREWES must be acknowledged in any resulting publications or
%    presentations
%
% 3) This SOFTWARE is provided "as is" with no warranty of any kind
%    either expressed or implied. CREWES makes no warranties or representation
%    as to its accuracy, completeness, or fitness for any purpose. CREWES
%    is under no obligation to provide support of any kind for this SOFTWARE.
%
% 4) CREWES periodically adds, changes, improves or updates this SOFTWARE without
%    notice. New versions will be made available at www.crewes.org .
%
% 5) Use this SOFTWARE at your own risk.
%
% END TERMS OF USE LICENSE

%Figure out varargin
if(nargin<6)
    mu=[1 1];
end


% preliminaries
 n=round(mlength/(t(2)-t(1)))+1;
 nn=floor(n/2);
 n=2*nn+1;%force n odd
 trin2=trin(:);
 trdsign=trdsign(:);
 dt=t(2)-t(1);
%  if(flag==1)
%      trdsign=trdsign.*mwhalf(length(trdsign));
%  else
%      trdsign=trdsign.*mwindow(length(trdsign));
%  end
 
% generate the Toeplitz matrix for the convolution equations
 TRIN= convmtx(trin2,n);
 
% generate a 2nd derivative operator
D=zeros(n,n);
for k=2:n-1
    D(k,k-1:k+1)=[1 -2 1];
end
D(1,1:2)=[1 -1];
D(n,n-1:n)=[1 -1];
% if(flag==0)
%     for k=2:n-1
%         D(k,k-1:k+1)=[1 -2 1];
%     end
%     D(1,1:2)=[1 -1];
%     D(n,n-1:n)=[1 -1];
% else
%     for k=1:n-2
%         D(k,k:k+2)=[1 -2 1];
%     end
%     D(n-1,n-1:n)=[1 -1];
%     D(n,n)=1;
% end
D2=D'*D;
%method='mwindow',c=30;
method='gaussian';c=2;
% solve the equations with left division
if icausal>=1
    nneg=icausal-1;%number of negative time samples
    RHS=[zeros(nneg,1);trdsign;zeros(n-1-nneg,1)];
    tm=dt*(-nneg:n-nneg-1)';
%     pct=30;
%     mw1=mwhalf(nneg+1,pct);
%     mw2=mwhalf(n-nneg-1,pct);
%     alpha=1-[flipud(mw1);mw2];%rather crude time constraint. Maybe a Gaussian?
    alpha=designtimeweights(n,nneg,method,c);
    %D=diag((tm/dt).^2);
    C=diag(alpha);
    C2=C*C;
    
    TRTR=TRIN'*TRIN;
    %A=max(abs(TRTR(:)));
    B=mu(1)*D2+mu(2)*C2+TRTR;
    %mfilt=pinv(B)*TRIN'*[trdsign;zeros(n-1,1)];
    mfilt=pinv(B)*TRIN'*RHS;

else
    nh=fix(n/2);
    RHS=[zeros(nh,1);trdsign;zeros(n-nh-1,1)];
    tm=(t(2)-t(1))*(-nh:nh)';
%     pct=30;
%     mw1=mwhalf(nh+1,pct);
%     mw2=mwhalf(n-nh-1,pct);
%     alpha=1-[flipud(mw1);mw2];
    %D=diag((tm/dt).^2);
    alpha=designtimeweights(n,nh,method,c);
    C=diag(alpha);
    C2=C*C;
    
    TRTR=TRIN'*TRIN;
    %A=max(abs(TRTR(:)));
    B=mu(1)*D2+mu(2)*C2+TRTR;
    mfilt=pinv(B)*TRIN'*RHS; %#ok<*MINV>
end

j=size(trin,1);
if j==1, mfilt=mfilt.'; tm=tm'; end
end

function wts=designtimeweights(n,nneg,method,c)
% n=total length of match filter
% nneg=number of negative time samples
% method= 'mwindow' or 'gaussian'
% c= for 'mwindow' c is the percent taper
%    for 'gaussian' c is #of sigmas at the ends
% wts = vector of wts that will become a diagonal matrix
switch method
    case 'mwindow'
        pct=c;
        mw1=mwhalf(nneg+1,pct);
        mw2=mwhalf(n-nneg-1,pct);
        wts=1-[flipud(mw1);mw2];
    case 'gaussian'
        tmp=zeros(n,1);
        sig1=nneg/c;
        x=(1:nneg+1)';
        tmp(x)=exp(-(x-nneg-1).^2/sig1^2);
        sig2=(n-nneg-1)/c;
        x=(nneg+1:n)';
        tmp(x)=exp(-(x-nneg-1).^2/sig2^2);
        wts=1-tmp;
end
        
        
end

   
 
 