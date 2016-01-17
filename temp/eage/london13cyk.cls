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
\RequirePackage{ifthen}
\RequirePackage{seg}
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






%%
%% This is file `geophysics.cls',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% geophysics.dtx  (with options: `class')
%% 
%% This is an automatically generated file. Do not edit it directly!
%% 
%%    Copyright (C) 2004 by Sergey Fomel
%% 
%%     This work may be distributed and/or modified under the
%%     conditions of the LaTeX Project Public License, either version 1.3
%%     of this license or (at your option) any later version.
%%     The latest version of this license is in
%%       http://www.latex-project.org/lppl.txt
%%     and version 1.3 or later is part of all distributions of LaTeX
%%     version 2003/12/01 or later.
%% 
%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K\L\M\N\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h\i\j\k\l\m\n\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Double quote  \"     Hash (number) \#
%%   Dollar        \$     Percent       \%     Ampersand     \&
%%   Acute accent  \'     Left paren    \(     Right paren   \)
%%   Asterisk      \*     Plus          \+     Comma         \,
%%   Minus         \-     Point         \.     Solidus       \/
%%   Colon         \:     Semicolon     \;     Less than     \<
%%   Equals        \=     Greater than  \>     Question mark \?
%%   Commercial at \@     Left bracket  \[     Backslash     \\
%%   Right bracket \]     Circumflex    \^     Underscore    \_
%%   Grave accent  \`     Left brace    \{     Vertical bar  \|
%%   Right brace   \}     Tilde         \~}

\newboolean{@manu}
\setboolean{@manu}{false}
\DeclareOption{manuscript}{\setboolean{@manu}{true}}
\DeclareOption{paper}{\setboolean{@manu}{false}}
\newboolean{@twoc}
\setboolean{@twoc}{false}
\DeclareOption{twocolumn}{\setboolean{@twoc}{true}}
\DeclareOption{onecolumn}{\setboolean{@twoc}{false}}
\newboolean{@shrt}
\setboolean{@shrt}{false}
\DeclareOption{short}{\setboolean{@shrt}{true}}
\DeclareOption{long}{\setboolean{@shrt}{false}}
\newcommand{\shortpaper}{\setboolean{@shrt}{true}}
\newcommand{\longpaper}{\setboolean{@shrt}{false}}
\newboolean{@eflt}
\setboolean{@eflt}{false}
\DeclareOption{endfloat}{\setboolean{@eflt}{true}}
\DeclareOption{noendfloat}{\setboolean{@eflt}{false}}
\newboolean{@revd}
\setboolean{@revd}{false}
\DeclareOption{revised}{\setboolean{@revd}{true}}
\DeclareOption{new}{\setboolean{@revd}{false}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ExecuteOptions{paper,onecolumn,long,noendfloat,new}
\ProcessOptions*

\ifthenelse{\boolean{@eflt}}{
  \ifthenelse{\boolean{@manu}}{\RequirePackage{endfloat}}{}
}{
\ifthenelse{\boolean{@manu}}{
  \renewcommand{\inputdir}[1]{\AtEndDocument{\renewcommand{\@path}{#1}}}}{}
\providecommand{\figdir}{.}
\providecommand{\setfigdir}[1]{\renewcommand{\figdir}{#1}}
\providecommand{\fig@file}[1]{\@path/Fig/#1}
\newcounter{@plots}
\setcounter{@plots}{0}
\ifthenelse{\boolean{@manu}}
{\newcommand{\fullplot}[4][X]{
    \stepcounter{@plots}
    \AtEndDocument{
      \begin{figure}[p]
        \centering
        \includegraphics[#3]{\fig@file{#2}}
        \caption [#4]{#4 \\
  \textbf{\seg@lhead} --
  \ifthenelse{\equal{\ms@number}{Manuscript}}{}{\ms@number}}
\label{fig:\@path-#2}
        \label{fig:#2}
      \end{figure}\clearpage
    }}
  \def\plot{\@ifstar{\fullplot}{\fullplot}}
  \setlength\abovecaptionskip{50\p@}
}
{\newcommand{\fullplot}[4][htbp]{%
    \begin{figure*}[#1]
      \centering
      \includegraphics[#3]{\fig@file{#2}}
      \caption{#4}
      \label{fig:\@path-#2}
      \label{fig:#2}
    \end{figure*}
  }}
\ifthenelse{\boolean{@manu}}
{\newcommand{\sideplot}[4][X]{\plot[#1]{#2}{#3}{#4}}}
{\newcommand{\sideplot}[4][htbp]{%
    \ifthenelse{\boolean{@twoc}}{%
      \begin{figure}[#1]
        \centering
        \includegraphics[#3]{\fig@file{#2}}
        \caption{#4}
\label{fig:\@path-#2}
        \label{fig:#2}
      \end{figure}
    }{%
      \begin{figure*}[#1]
        \makebox{%
          \begin{minipage}{0.4\textwidth}
            \caption{#4}
    \label{fig:\@path-#2}
            \label{fig:#2}
          \end{minipage}
          \hfill
          \begin{minipage}{0.6\textwidth}
            \centering
            \includegraphics[#3]{\fig@file{#2}}
          \end{minipage}}
      \end{figure*}
    }
  }
  \ifthenelse{\boolean{@twoc}}{%
    \def\plot{\@ifstar{\fullplot}{\sideplot}}}{%
    \def\plot{\@ifstar{\fullplot}{\fullplot}}}
}
\RequirePackage{subfigure}
\def\next@item#1,#2?{#1}
\def\rest@item#1,#2?{#2}
\newcounter{sub@fig}
\ifthenelse{\boolean{@manu}}
{\newcommand{\seg@multiplot}[5][X]{
    \stepcounter{@plots}
    \AtEndDocument{
      \begin{figure}[p]
        \centering
        \setcounter{sub@fig}{0}
        \edef\list@i{#3}
        \loop
        \edef\item@i{\expandafter\next@item\list@i,\empty?}
        \edef\list@i{\expandafter\rest@item\list@i,\empty?}
        \ifx\item@i\empty
        \else
        \stepcounter{sub@fig}
        \subfigure[]{\includegraphics[#4]{\fig@file{\item@i}}%
  \label{fig:\@path-\item@i}
  \label{fig:\item@i}}
        \ifthenelse{\value{sub@fig} = #2}{\\ \setcounter{sub@fig}{0}}{}
        \repeat
        \caption[#5]{#5 \\
  \textbf{\seg@lhead} --
  \ifthenelse{\equal{\ms@number}{Manuscript}}{}{\ms@number}}
\label{fig:\@path-#3}
        \label{fig:#3}
      \end{figure}\clearpage
}}
\def\multiplot{\@ifstar{\seg@multiplot}{\seg@multiplot}}
}{%
  \newcommand{\seg@fullmultiplot}[5][htbp]{
    \begin{figure*}[#1]
      \centering
      \setcounter{sub@fig}{0}
      \edef\list@i{#3}
      \loop
      \edef\item@i{\expandafter\next@item\list@i,\empty?}
      \edef\list@i{\expandafter\rest@item\list@i,\empty?}
      \ifx\item@i\empty
      \else
      \stepcounter{sub@fig}
      \subfigure[]{\includegraphics[#4]{\fig@file{\item@i}}%
\label{fig:\@path-\item@i}
        \label{fig:\item@i}}
      \ifthenelse{\value{sub@fig} = #2}{\\ \setcounter{sub@fig}{0}}{}
      \repeat
      \caption{#5}
      \label{fig:\@path-#3}
      \label{fig:#3}
    \end{figure*}
  }
  \newcommand{\seg@sidemultiplot}[5][htbp]{
    \begin{figure}[#1]
      \centering
      \setcounter{sub@fig}{0}
      \edef\list@i{#3}
      \loop
      \edef\item@i{\expandafter\next@item\list@i,\empty?}
      \edef\list@i{\expandafter\rest@item\list@i,\empty?}
      \ifx\item@i\empty
      \else
      \stepcounter{sub@fig}
      \subfigure[]{\includegraphics[#4]{\fig@file{\item@i}}%
\label{fig:\@path-\item@i}
        \label{fig:\item@i}}
      \ifthenelse{\value{sub@fig} = #2}{\\ \setcounter{sub@fig}{0}}{}
      \repeat
      \caption{#5}
      \label{fig:\@path-#3}
      \label{fig:#3}
    \end{figure}
  }
  \ifthenelse{\boolean{@twoc}}{%
    \def\multiplot{\@ifstar{\seg@fullmultiplot}{\seg@sidemultiplot}}}{%
    \def\multiplot{\@ifstar{\seg@fullmultiplot}{\seg@fullmultiplot}}}
}
\newcounter{@tabls}
\setcounter{@tabls}{0}
\ifthenelse{\boolean{@manu}}{
  \providecommand{\tabl}[4][X]{
    \stepcounter{@tabls}
    \AtEndDocument{
      \begin{table}[p]
        #4
        \caption{#3}
        \label{tbl:#2}
      \end{table}\clearpage
    }
  }
}{
  \providecommand{\seg@sidetabl}[4][htbp]{
    \begin{table}[#1]
      #4
      \caption{#3}
      \label{tbl:#2}
    \end{table}
  }
  \providecommand{\seg@fulltabl}[4][htbp]{
    \begin{table*}[#1]
      #4
      \caption{#3}
      \label{tbl:#2}
    \end{table*}
  }
  \def\tabl{\@ifstar{\seg@fulltabl}{\seg@sidetabl}}
}
\renewcommand*{\l@figure}[2]{%
  \setlength{\@tempdima}{2.3em}%
  \noindent\hspace*{1.5em}#1\hfil\newline}
\let\l@table\l@figure
\renewcommand\listoffigures{%
  \section{\listfigurename
    \@mkboth{\MakeUppercase\listfigurename}%
    {\MakeUppercase\listfigurename}}%
  \@starttoc{lof}%
}
\renewcommand\listoftables{%
  \section{\listtablename
    \@mkboth{\MakeUppercase\listtablename}%
    {\MakeUppercase\listtablename}}%
  \@starttoc{lot}%
}
\AtEndDocument{
  \ifthenelse{\value{@tabls} > 0}{\newpage\listoftables}{}
  \ifthenelse{\value{@plots} > 0}{\newpage\listoffigures}{}
}}
\ifthenelse{\boolean{@revd}}{
  \RequirePackage{color}%
  \RequirePackage{ulem}%
  \normalem
}{}

\endinput
%%
%% End of file `geophysics.cls'.
