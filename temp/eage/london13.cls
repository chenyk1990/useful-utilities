%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is file 'london13.cls' %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Made by: Bobby Hak
%% Email: bhak@fugro-jason.com
%%
%% This class has been tested with:
%%     MikTeX version 2.9
%%
%% It is based on the standard LaTeX2e "article" class and some
%% standard packages:
%%    times    - for the font Times New Roman, 11pt
%%    mathptmx - for the font Times New Romen with math support
%%    fancyhdr - for defining headnote and footnote
%%    graphicx - for including logo
%%    geometry - for putting page dimensions correct
%%    natbib   - for referencing
%%
%% All packages should be availble in any normal distribution of LaTeX,
%% or can be downloaded from the CTAN homepage:
%%    www.ctan.org
%%
\NeedsTeXFormat{LaTeX2e}[1996/12/01]
\ProvidesClass{london13}[2011/07/20 v1.0 EAGE Abstract Class]
%%
%% Pass all options directly to 'article'
%%
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions\relax
%%
%% prefered fontsize, papersize, and type
%%
\LoadClass[oneside,11pt,a4paper]{article}
%%
%% prefered fonts: Times New Roman becomes the standard font (text)
%%
\RequirePackage[T1]{fontenc}
\RequirePackage{times}
\RequirePackage{mathptmx}
%%
%% Page size formats
%%
\RequirePackage[a4paper]{geometry}
\geometry{left=1.0in,right=1.0in,top=1.0in,bottom=1.0in}
\geometry{head=0.68in,headsep=0.09in,foot=0.5in}
\addtolength{\footskip}{-1.0\baselineskip}
%%
%% some others
%%
\setcounter{secnumdepth}{0} % no numbering of sections
\setlength{\parindent}{0pt} % no indent for new paragraphs
\setlength{\parskip}{11pt} % skip one line between paragraphs
\setlength{\partopsep}{0pt} % skip between environments and paragraphs
\setlength{\topsep}{0pt} % skip between environments and paragraphs
%%
%% bibliography
%%
\RequirePackage{natbib} % special citing commands
\gdef\bibfont{\fontsize{10}{10}\selectfont}
%\newcommand{\bibfont}{\fontsize{10}{10}\selectfont} % bibliography fontsize
\setlength{\bibsep}{0pt} % no lines between bibitems
\setlength{\bibhang}{1em} % indentation of second line in bibitem
\bibliographystyle{firstbreak} % stylefile
%%
%% set header and footer
%%
\RequirePackage{fancyhdr}
\RequirePackage{graphicx}
\pagestyle{fancy} % load package
\fancyhf{} % clear styles
\fancyhead[R]{\includegraphics[height=0.47in,width=2.06in]{london13.eps}} % input logo top right, even and odd side
\fancyfoot[C]{\vspace{-1pt}\fontsize{10}{10}\selectfont{%
75\textsuperscript{th} EAGE Conference \& Exhibition incorporating SPE EUROPEC 2013\\
London, UK, 10-13 June 2013}} % footnote
\renewcommand{\headrulewidth}{0pt} % no line in header
\renewcommand{\footrulewidth}{0.4pt} % add line in footnote
\pagenumbering{gobble} % define page numbering
%%
%% section headers and paragraphs
%%
\renewcommand{\section}{\@startsection{section}%
{1}% level
{0pt}% indent
{.1ex plus -.1ex}% beforeskip
{.1ex plus -.1ex}% afterskip
{\fontsize{11}{11}\selectfont\bf}} % font
%
\renewcommand{\subsection}{\@startsection{subsection}%
{2}% level
{0pt}% indent
{.1ex plus -.1ex}% beforeskip
{.1ex plus -.1ex}% afterskip
{\fontsize{11}{11}\selectfont\emph}} % font
%%
%% caption font
%%
\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{{\bf\itshape #1} \itshape #2}%
  \ifdim \wd\@tempboxa >\hsize
    {{\bf\itshape #1} \itshape #2\par}
  \else
    \hbox to\hsize{\hfil\box\@tempboxa\hfil}%
  \fi
  \vskip\belowcaptionskip}
% EOF
