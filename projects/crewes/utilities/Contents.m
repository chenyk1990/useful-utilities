% CREWES Utilites toolbox
%
%Logical tests
% BETWEEN - logical test, finds sampls in vector between given bounds
% INSIDEPOLY - identify points inside a polygonal region in 2D.
% ISCOMPLEX - logical test for presence of complex numbers
% NEAR - return indices of those samples in a vector nearest given bounds
% NOTBETWEEN - logical test, the negation of BETWEEN
% SURROUND - analyze how a vector surrounds some test points
% WITHIN - test a point to see if it is inside a polygon or not
%
%Graphical
% CLEARLINES - clear (delete) the rays in a figure
% CLICKALINE - simple utility to draw a line by clicking
% COLORVIEW - puts up an interactive widget to manipulate colormaps
% DRAGLINE - simple dragging of a line with constraints. Better than moveline.
% DRAWLINEINIT - initialize line drawing
% DRAWLINEFINI - finalize line drawing
% DRAWLINE - draw a line on a figure window
% DRAWPICK - draw a line on top of an image (e.g. seismic)
% FIGCENT - create a new figure in the center of the screen
% FIGPOSN - shift a figure to somewhere else
% FIGSIZE - set figure size as a fraction of screen size
% HARDZOOM - subset a matrix according to axis limit specs
% MOVELINE - Utility to move a line by clicking and dragging
% SCA - set current axis utility
% SELBOX - draw a selection box on a figure window
% SELBOXFINI - finalize selection box drawing
% SELBOXINIT - initialize selection box drawing
% SIGNATURE - put a signature on a figure window
% SIMPLEDIT - can be used for simple graphical editing of line graphs (see EDITLINES)
% SIMPLEDIT_DEMO - illustrates the use of simpledit
% SIMPLEZOOM - figure zooming utility using SELBOX
% XTICK - adjust the location of tick marks on the x axis
% YTICK - adjust the location of tick marks on the y axis
% YTICKLBLOFF - turn off ytick labels
%
%Interpolation
% PCINT - piecewise constant interpolation
% PWLINT - piecewise linear interpolation (much faster than interp1)
% INTERPEXTRAP - similar to Matlab's interp1 but always linear and extrapolates as needed
%
%Windowing, padding, trace conversions
% BOXKAR - create a boxcar window with optional taper
% FROMDB - convert from (db,phase) to (real,imaginary) (see TODB)
% GAUSSIAN - create a Gaussian window
% GAUS_RADIAL - 2D radially symmetric Gaussian
% HILBM - Hilbert transform (only use if signal toolbox is unavailable)
% GWINDOW - like mwindow except that you get a Gaussian
% MWINDOW - creates an mwindow (boxcar with raised-cosine tapers)
% MWHALF - half an mwindow (boxcar with raised-cosing taper on one end)
% PAD_TRACE - pads (truncates) one trace with zeros to be the length of another
% PADPOW2 - pad a trace with zeros to the next power of 2 in legth
% TODB - converts from (real,imaginary) to (decibels, phase)
% TRIANGLE - create a triangle window
%
%Other stuff
% BOXKAR - returns a boxcar window
% GAUSS - returns a gaussian distribution sampled in frequency
% NUM2STRMAT - convert a vector of numbers to a string matrix
% NUM2STRCELL - convert a vector of numbers to a cell vector of strings.
%     **** note, the opposite of NUM2STRCELL is vec=str2double(str) . No special function needed.
% TIME2STR - works like num2str except always shows to the millisecond
% SLICEMAT - slices a matrix along a trajectory with a fairway
% TIME2STR - convert a time to a string with decimals to the nearest millisecond
% XCOORD - create a coordinate vector given start, increment, and number
% SAMPLE_STARTUP - example of a startup file to use with the CREWES toolbox
% SIGFIG - round a number to a given number of significant figures 
%
%Other stuff: RJ Ferguson:
% BLTIFFT - returns time-domain data given a positive sided, band-limited spectrum
% TI_IMPULSE - plots and returns 2-D coordinates of a migration impulse response 
%            for a homogeneous TI medium
% 
%Other stuff: KWH
% CHECK_INT_RANGE - checks if input data can be stored as an integer
% CHECK_FP_RANGE - checks if input data can be stored as floating-point
% CHECK_IBMFP_RANGE - checks if input data can be stored as IBM floating-point
% CRISGRAPHICS - return logical true if h is a graphics handle or logical
%            false if it is not
% CRVERSION - returns CREWES toolbox version (exists in toolboxes downloaded 
%            from website
% FILE - File object for opening and closing binary files for fread/fwrite
% FINDDIR - recursively finds a directory and returns the full path to it
% TEXTBAR - display/update a command line progress bar if called recursively
%
