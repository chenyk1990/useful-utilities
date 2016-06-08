 function Mada_matlab(infile,outfile,inpara1,inpara2)
% Author      : Yangkang Chen
%               Texas Consortium of Computational Seismology
%               Jackson School of Geosciences
%               The University of Texas at Austin
%         
% Date        : Oct, 2014

% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
    
%  Copyright (C) 2014 The University of Texas at Austin
%  Copyright (C) 2014 Yangkang Chen
%  
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%  
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  

%% default values 
para1=
para2=

%%from Madagascar to Matlab
% create memory
din=zeros(n1,n2);
rsf_read(din,infile)

%% Main program goes here !


%% from Matlab to Madagascar
rsf_create(outfile,size(dout)');
rsf_write(dout,outfile);


