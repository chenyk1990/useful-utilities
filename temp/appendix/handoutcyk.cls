\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{handout}[2003/09/01 v0.1 handouts for imaging class]

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions*

\LoadClass{article}

\RequirePackage{times}
\RequirePackage{mathptm}
\RequirePackage{ifthen}
\RequirePackage{seg}

\setlength{\textheight}{9.5in}
\setlength{\textwidth}{6.5in}
\setlength{\parskip}{1.5ex plus0.1ex minus0.1ex} % should be a *rubber* length
\setlength{\oddsidemargin}{.01in}
\setlength{\evensidemargin}{.01in}
\setlength{\topmargin}{-0.51in}
\setlength{\textfloatsep}{30pt plus 3pt minus 6pt}
\setlength{\headsep}{0in}

\newcommand{\@name}{Yangkang Chen}
\renewcommand{\author}[1]{\renewcommand{\@name}{#1}}

\renewcommand{\title}[1]{%
  \begin{center}
{\small
  #1 \hfill \@name} \\ 
\hrulefill
\end{center}}

% to make section headings smaller
\renewcommand{\Large}{\protect\normalsize\centering}%
\renewcommand{\large}{\protect\normalsize}

\let\sep@section\section
\let\sep@subsection\subsection
\let\sep@subsubsection\subsubsection
\newcommand{\sepsection}[1]{\sep@section*{#1}}
\renewcommand{\section}[2][]{\sepsection{\MakeUppercase{#2}}}%
\def\subsection{\sep@subsection*}%
\def\subsubsection{\sep@subsubsection*}

\RequirePackage{natbib}
\DeclareRobustCommand{\refname}{REFERENCES}
\DeclareRobustCommand{\bibname}{\protect\large~REFERENCES}
\renewcommand{\bibsection}{\section{\refname}}

\RequirePackage{graphicx}

\newcommand{\figdir}{.}


\newcommand{\plot}[4][htbp]{%
\begin{figure*}[#1]
\centering
\includegraphics[#3]{\figdir/Fig/#2}
\caption{#4}
\label{fig:#2}
\end{figure*}}




\newlength{\box@width}
\setlength{\box@width}{\textwidth}
\addtolength{\box@width}{-20pt} 
\setlength{\fboxsep}{10pt}

\newcommand{\boxit}[1]{\noindent\fbox{\parbox{\box@width}{#1}}}



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
  
  
  
  
  
  
