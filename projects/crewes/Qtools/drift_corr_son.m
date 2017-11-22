function sond=drift_corr_son(son,q,f,f0)
%
% sond=drift_corr_son(son,q,f,f0)
%
% Method: Convert the sonic values to interval velocity, then correct the
% interval velocities for Q dispersion via the Azimi formula (see function
% velf).
%
% son ... input sonic. Should not be check-shot corrected.
% q ... Q. Either a single scalar or a vector the same size as son of local
%       (interval) Q values.
% f ... dominant frequency (Hz) of seismic data to which sonic is to be compared
% *********** default = 50 **********
% f0 ... dominant frequency of the logging tool (the reference frequency)
% *********** default =12500 ***********
%

if(length(q)==1)
    Q=q*ones(size(son));
else
    Q=q;
end

if(nargin<3)
    f=50;
end
if(nargin<4)
    f0=12500;
end

v=10^6 ./son;

vf=velf(v,Q,f,f0);

sond=10^6 ./vf;
