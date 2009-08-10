%% From Ken Parker (kparker@facility.cs.utah.edu) Fri May  5 16:34:08 1995
%% uuthesis2e.cls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\ProvidesClass{uuthesis2e}
              [1998/04/01 v1.3 University of Utah Thesis Class]
\newcommand\@ptsize{}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\typeout{Document Style `uuthesis2e'. Version 1.8d - 6 Apr 1996}
\def\@warning#1{\typeout{LaTeX Warning [l.\the\inputlineno]: #1.}}

%%%%%%%%%%%%%%%%%%%%
% 11pt and 12pt are the only possibilities for type size - default is 12pt.
% These defns relate document style options to style file names loaded later.

\DeclareOption{11pt}{\renewcommand\@ptsize{1}}
\DeclareOption{12pt}{\renewcommand\@ptsize{2}}

%%% Used to use \input size\@ptsize.clo, not now. See below.

%%%%%%%%%%%%%%%%%%%%
% Choice a style guild.

\DeclareOption{Chicago}{
  \@ifundefined{@StyleGuild}
    {\gdef\@StyleGuild{Chicago}}
    {\@latex@error{Can't specify more than one style guild}\@eha}
}
\DeclareOption{MLA}{
  \@ifundefined{@StyleGuild}
    {\gdef\@StyleGuild{MLA}}
    {\@latex@error{Can't specify more than one style guild}\@eha}
}
\DeclareOption{Turabian}{
  \@ifundefined{@StyleGuild}
    {\gdef\@StyleGuild{Turabian}}
    {\@latex@error{Can't specify more than one style guild}\@eha}
}

%%%%%%%%%%%%%%%%%%%%
% "twoside" printing option not for theses, but maybe for drafts?

\DeclareOption{oneside}{\@twosidefalse \@mparswitchfalse}
\DeclareOption{twoside}{\@twosidetrue  \@mparswitchtrue}

%%%%%%%%%%%%%%%%%%%%
% "draft" option relaxes pickiness about overfull hboxes.

\DeclareOption{draft}{\setlength\overfullrule{5pt}}
\DeclareOption{final}{\setlength\overfullrule{0pt}}

%%%%%%%%%%%%%%%%%%%%
% "report" option does single-spacing and other stuff suitable for a
% tech report version of the thesis.  The option just sets a flag that
% everybody else will look at.

\newif\if@report \@reportfalse
\DeclareOption{reoprt}{\@reporttrue}

%%%%%%%%%%%%%%%%%%%%
% Other stuff, not necessarly working

\newif\if@titlepage
\DeclareOption{titlepage}{\@titlepagetrue}
\DeclareOption{notitlepage}{\@titlepagefalse}

\newif\if@openright
\DeclareOption{openright}{\@openrighttrue}
\DeclareOption{openany}{\@openrightfalse}

\DeclareOption{onecolumn}{\@twocolumnfalse}
\DeclareOption{twocolumn}{\@twocolumntrue}
\DeclareOption{leqno}{\input{leqno.clo}}
\DeclareOption{fleqn}{\input{fleqn.clo}}

\newif\if@honors\@honorsfalse
\DeclareOption{honors}{\@honorstrue}

\newif\if@csreport\@csreportfalse
\DeclareOption{csreport}{\@csreporttrue}

\newif\if@openbib \@openbibfalse
\DeclareOption{openbib}{\@openbibtrue}

\def\reportitle{Define control sequence reporttitle to produce titlepage}

\ExecuteOptions{letterpaper,12pt,oneside,onecolumn,final,openany}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\ProcessOptions

%% Check that the user selected a style guild.
\@ifundefined{@StyleGuild}
  {\@latex@error{You must specify a style guild.}\@ehc}
  {}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\input{\@StyleGuild.clo}
\input{size1\@ptsize.clo}

\if@csreport\global\def\reporttitle{\CSREPORTTITLE}\fi
\if@honors\global\def\reporttitle{\HONORSTITLE}\fi



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thesis 12-point type is being used.
% Each size-changing command \SIZE executes the command
%        \@setsize\SIZE{BASELINESKIP}\FONTSIZE\@FONTSIZE
% where:
%   BASELINESKIP = Normal value of \baselineskip for that size.  (Actual
%                  value will be \baselinestretch * BASELINESKIP.)
%
%  \FONTSIZE     = Name of font-size command.  The currently available
%                  (preloaded) font sizes are: \vpt (5pt), \vipt (6pt),
%                  \viipt (etc.), \viiipt, \ixpt, \xpt, \xipt, \xiipt,
%                  \xivpt, \xviipt, \xxpt, \xxvpt.
%  \@FONTSIZE    = The same as the font-size command except with an
%                  '@' in front---e.g., if \FONTSIZE = \xivpt then
%                  \@FONTSIZE = \@xivpt.
%
% For reasons of efficiency that needn't concern the designer,
% the document style defines \@normalsize instead of \normalsize .  This is
% done only for \normalsize, not for any other size-changing commands.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Set all the basic size parameters of the UU thesis for 12pt
%% 1.2em= 14.09981pt, 2em= 23.49976pt in 12pt font
\lineskip 1pt           % \lineskip is 1pt for all font sizes.
\normallineskip 1pt

\def\@optionONE{2}
\if\@optionONE\@ptsize
\typeout{Using 12pt option}

\def\singlespace{\baselineskip=13.5pt}
\def\doublespace{\baselineskip=21pt}

% Define sectional baseline skip for each sectional font
%
%% Sometimes the space after a main heading is too large or too small.
%% This is created by a chapter command followed by a section command,
%% or in some cases, a table insertion where the caption interferes with
%% the chapter heading. Use \fixchapterheading to remove extra space.
%
\def\doublespacedheadings{
\def\HFmainheadHT{24pt}
\def\HFpartHT{24pt}
\def\HFchapterHT{24pt}
\def\HFsectionHT{20pt}
\def\HFsubsectionHT{20pt}
\def\HFsubsubsectionHT{20pt}
\def\HFsubsubsubsectionHT{20pt}
\def\HFparagraphHT{14.5pt}
% Skips
\def\HFmainheadSKIP{24pt}
\def\HFpartSKIP{14.5pt}
\def\HFchapterSKIP{14.5pt}
\def\HFsectionSKIP{14.5pt}
\def\HFsubsectionSKIP{14.5pt}
\def\HFsubsubsectionSKIP{15.5pt}
\def\HFsubsubsubsectionSKIP{15.5pt}
\def\fixchapterheading{\vspace{-18pt}}
\def\fixmainheadingSKIP{\vskip 12pt}
}

\def\singlespacedheadings{
\def\HFmainheadHT{24pt}
\def\HFpartHT{24pt}
\def\HFchapterHT{24pt}
\def\HFsectionHT{15.5pt}
\def\HFsubsectionHT{13.5pt}
\def\HFsubsubsectionHT{13.5pt}
\def\HFsubsubsubsectionHT{13.5pt}
\def\HFparagraphHT{14.5pt}
% skips
\def\HFmainheadSKIP{24pt}
\def\HFpartSKIP{14.5pt}
\def\HFchapterSKIP{14.5pt}
\def\HFsectionSKIP{21pt}
\def\HFsubsectionSKIP{21pt}
\def\HFsubsubsectionSKIP{21pt}
\def\HFsubsubsubsectionSKIP{21pt}
\def\fixchapterheading{\vspace{-24pt}}
\def\fixmainheadingSKIP{\vskip 12pt}
}

% Define special fonts for sectional headings
%
\def\HFmainhead{\@setsize\HFmainhead{\HFmainheadHT}\xivpt\@xivpt}
\def\HFpart{\@setsize\HFpart{\HFpartHT}\xivpt\@xivpt}
\def\HFchapter{\@setsize\HFchapter{\HFchapterHT}\xivpt\@xivpt}
\def\HFsection{\@setsize\HFsection{\HFsectionHT}\xivpt\@xivpt}
\def\HFsubsection{\@setsize\HFsubsection{\HFsubsectionHT}\xiipt\@xiipt}
\def\HFsubsubsection{\@setsize\HFsubsubsection{\HFsubsubsectionHT}\xiipt\@xiipt}
\def\HFsubsubsubsection{\@setsize\HFsubsubsubsection{\HFsubsubsubsectionHT}\xiipt\@xiipt}
\def\HFparagraph{\@setsize\HFparagraph{\HFparagraphHT}\xiipt\@xiipt}
\def\HFapproval{\@setsize\HFapproval{13pt}\xipt\@xipt}
\def\HFapprovalsmall{\@setsize\HFapprovalsmall{12pt}\xpt\@xpt}
\fi
%% End of all 12pt specific stuff

%% Begin 11pt specific settings

\def\@optionONE{1}
\if\@optionONE\@ptsize
\typeout{Using 11pt option}
\def\singlespace{\baselineskip=13pt}
\def\doublespace{\baselineskip=19pt}

\def\doublespacedheadings{
\def\HFmainheadHT{24pt}
\def\HFpartHT{24pt}
\def\HFchapterHT{24pt}
\def\HFsectionHT{20pt}
\def\HFsubsectionHT{20pt}
\def\HFsubsubsectionHT{20pt}
\def\HFsubsubsubsectionHT{20pt}
\def\HFparagraphHT{13.5pt}
% Skips
\def\HFmainheadSKIP{24pt}
\def\HFpartSKIP{14.5pt}
\def\HFchapterSKIP{14.5pt}
\def\HFsectionSKIP{14.5pt}
\def\HFsubsectionSKIP{14.5pt}
\def\HFsubsubsectionSKIP{15.5pt}
\def\HFsubsubsubsectionSKIP{15.5pt}
\def\fixchapterheading{\vspace{-18pt}}
\def\fixmainheadingSKIP{\vskip 12pt}
}

\def\singlespacedheadings{
  \def\HFmainheadHT{24pt}
  \def\HFpartHT{24pt}
  \def\HFchapterHT{24pt}
  \def\HFsectionHT{15.5pt}
  \def\HFsubsectionHT{13.5pt}
  \def\HFsubsubsectionHT{13.5pt}
  \def\HFsubsubsubsectionHT{13.5pt}
  \def\HFparagraphHT{13.5pt}
% skips
  \def\HFmainheadSKIP{24pt}
  \def\HFpartSKIP{14.5pt}
  \def\HFchapterSKIP{14.5pt}
  \def\HFsectionSKIP{21pt}
  \def\HFsubsectionSKIP{21pt}
  \def\HFsubsubsectionSKIP{21pt}
  \def\HFsubsubsubsectionSKIP{21pt}
  \def\fixchapterheading{\vspace{-24pt}}
  \def\fixmainheadingSKIP{\vskip 12pt}
}

% Define special fonts for sectional headings
%
\def\HFmainhead{\@setsize\HFmainhead{\HFmainheadHT}\xivpt\@xivpt}
\def\HFpart{\@setsize\HFpart{\HFpartHT}\xivpt\@xivpt}
\def\HFchapter{\@setsize\HFchapter{\HFchapterHT}\xivpt\@xivpt}
\def\HFsection{\@setsize\HFsection{\HFsectionHT}\xivpt\@xivpt}
\def\HFsubsection{\@setsize\HFsubsection{\HFsubsectionHT}\xiipt\@xiipt}
\def\HFsubsubsection{\@setsize\HFsubsubsection{\HFsubsubsectionHT}\xiipt\@xiipt}
\def\HFsubsubsubsection{\@setsize\HFsubsubsubsection{\HFsubsubsubsectionHT}\xiipt\@xiipt}
\def\HFparagraph{\@setsize\HFparagraph{\HFparagraphHT}\xipt\@xipt}
\def\HFapproval{\@setsize\HFapproval{13pt}\xipt\@xipt}
\def\HFapprovalsmall{\@setsize\HFapprovalsmall{12pt}\xpt\@xpt}

\typeout{Option 11pt uses 5pt for 'tiny': microfilm needs 2mm.}
\fi
%% End of 11pt settings

% Default is single-spaced headings
\singlespacedheadings

%% Document Style Parameters.

% Many specific numbers are in the type size files, while others are
% absolute distance mandated by the thesis format itself.
% The names of type-size-specific parameters are given here, but commented out.

%% Page Style Parameters.

\oddsidemargin 0.5in
\evensidemargin 0in     % does this work right with \textwidth below?

\marginparwidth = 0pt   % theses don't have marginal notes
%\marginparsep

\topmargin -0.25in
\headheight 1em
\headsep 2em
\textheight 8.7in       % 9in - 1em more or less
\textwidth 5.85in       % 6in - 1em more or less
\topskip 10pt           % default from rep12.sty (dubious)

\newdimen\footheight
\footheight 2em         % page foot is used by front matter page numbers
\footskip 10pt

%% Paragraph Style Parameters.

\parindent 1.5em        % Width of paragraph indentation.
\def\baselinestretch{1} % \baselinestretch * \baselineskip = \normalbaselineskip

%% Displayed Paragraph Style Parameters.

% Fix by Nelson Beebe Nov 1992
%% NB: None of these may have glue; otherwise, spacing around section
%% headings violates thesis guidelines!
\topsep=10pt
\partopsep=3pt
\itemsep=5pt
\parsep=5pt
\parskip=0pt

% The following commands are used to set the default values for the list
% environment's parameters. See the LaTeX manual for an explanation of
% the meanings of the parameters.  Defaults for the list environment are
% set as follows.  First, \rightmargin, \listparindent and \itemindent
% are set to 0pt.  Then, for a Kth level list, the command \@listK is
% called, where 'K' denotes 'i', 'ii', ... , 'vi'.  (I.e., \@listiii is
% called for a third-level list.)  By convention, \@listK should set
% \leftmargin to \leftmarginK.
%
% For efficiency, level-one list's values are defined at top level, and
% \@listi is defined to set only \leftmargin.

\leftmargini 2.5em
\leftmarginii 2.2em     % > \labelsep + width of '(m)'
\leftmarginiii 1.87em   % > \labelsep + width of 'vii.'
\leftmarginiv 1.7em     % > \labelsep + width of 'M.'
\leftmarginv 1em        % is 1em a random value?
\leftmarginvi 1em       % random value?

\def\@listi{\leftmargin\leftmargini}

\def\@listii{\leftmargin \leftmarginii
 \labelwidth \leftmarginii \advance\labelwidth-\labelsep
 \topsep 5pt plus 2.5pt minus 1pt
 \parsep 2.5pt plus 1pt minus 1pt
 \itemsep \parsep
}

\def\@listiii{\leftmargin \leftmarginiii
  \labelwidth \leftmarginiii \advance\labelwidth-\labelsep
  \topsep 2.5pt plus 1pt minus 1pt
  \parsep \z@ \partopsep 1pt plus 0pt minus 1pt
  \itemsep \topsep
}

\def\@listiv{\leftmargin \leftmarginiv
  \labelwidth \leftmarginiv \advance\labelwidth-\labelsep
}

\def\@listv{\leftmargin \leftmarginv
  \labelwidth \leftmarginv \advance\labelwidth-\labelsep
}

\def\@listvi{\leftmargin \leftmarginvi
  \labelwidth \leftmarginvi \advance\labelwidth-\labelsep
}

\leftmargin \leftmargini
%\rightmargin
%\listparindent
%\itemindent
\labelsep .5em
\labelwidth \leftmargini \advance\labelwidth-\labelsep

%% Math Mode Style Parameters.

%\jot
%\mathindent
%...

%% Tabbing Style Parameters.

\tabbingsep \labelsep

%% Array Style Parameters.

% The numbers are conventional, since thesis format does not say much.

\arraycolsep 5pt
\tabcolsep 6pt
\arrayrulewidth .4pt
\doublerulesep 2pt
\def\arraystretch{1}

%% Box Style Parameters.

% Again, numbers are conventional.

\fboxrule .4pt
\fboxsep 3pt

%% Page-Breaking Penalties.  (Mysterious numbers, fitted by experiment.)

\@lowpenalty   51      % Produced by \nopagebreak[1] or \nolinebreak[1]
\@medpenalty  151      % Produced by \nopagebreak[2] or \nolinebreak[2]
\@highpenalty 301      % Produced by \nopagebreak[3] or \nolinebreak[3]

\@beginparpenalty -\@lowpenalty    % Before a list or paragraph environment.
\@endparpenalty   -\@lowpenalty    % After a list or paragraph environment.
\@itempenalty     -\@lowpenalty    % Between list items.

\hyphenpenalty=400
\exhyphenpenalty=400

\clubpenalty=4000       % Try to eliminate club lines.
\widowpenalty=4000      % Ditto for widows.
\displaywidowpenalty=3000       % Again, for math display widow line.
\predisplaypenalty=3000         % Breaking before a math display.
\postdisplaypenalty=2000        % Breaking after a math display.

% Committee approval page format depends on committee size (flags set
% automatically by use of \fourthreader etc commands).

\newif\iffourth
\newif\iffifth

\fourthfalse
\fifthfalse

% Only do a dedication page if a dedication has been supplied (auto set by
% \dedication command).

\newif\ifdedication

\dedicationfalse

% Correct front matter requires knowing if lists of figures and tables present.
% Unfortunately these have to be under manual control.

\newif\iflistoffigures
\newif\iflistoftables

\listoffigurestrue
\listoftablestrue

% For controlling the various chit-chat that might be issued.

\newif\ifnoisy \noisytrue

%%% Define the macros that the user will define everything with.

% \title already defined.
% \author already defined.
\def\title#1{\gdef\@title{\uppercase{#1}}}      % 17-Apr-1992 to uppercase
\def\thesistype#1{\gdef\@thesistype{#1}}        % thesis or dissertation
\def\graduatedean#1{\gdef\@graduatedean{#1}}    % dean of the graduate school
\def\department#1{\gdef\@department{#1}}        % academic department
\def\degree#1{\gdef\@degree{#1}}                % exact degree being gotten
\def\committeechair#1{\gdef\@committeechair{#1}} % Committee chair-critter.
\def\firstreader#1{\gdef\@firstreader{#1}}      % The other committee members.
\def\secondreader#1{\gdef\@secondreader{#1}}
\def\thirdreader#1{\gdef\@thirdreader{#1}\fourthtrue}
\def\fourthreader#1{\gdef\@fourthreader{#1}\fifthtrue}
\def\departmentchair#1{\gdef\@departmentchair{#1}} % Department chair-person.
\def\chairtitle#1{\gdef\@chairtitle{#1}}        % Title of chair-person.
\def\submitdate#1{\gdef\@submitdate{#1}}        % Submission date.
\def\copyrightyear#1{\gdef\@copyrightyear{#1}}  % Copyright year.
\def\dedication#1{\gdef\@dedication{#1}\dedicationtrue} % The dedication.
\def\reportnumber#1{\gdef\@reportnumber{#1}}    % Number/id of tech report.
\def\honorsdepartment#1{\gdef\@honorsdepartment{#1}}
\def\honorsadvisor#1{\gdef\@honorsadvisor{#1}}
\def\honorssupervisor#1{\gdef\@honorssupervisor{#1}}
\def\honorsdirector#1{\gdef\@honorsdirector{#1}}

% Number of appendices could be computed I suppose, but probably hard and
% (more importantly) unreliable.  Correct TOC format requires this number.

\newcount\numberofappendices
\numberofappendices=0

% Constrain the submission dates to be only the four possible months for
% a valid thesis submission.

\def\@submitdate{\ifcase\the\month\or March\or June\or August\or December\fi
  \space\number\the\year}

%%% Define all our counters.  (Which style guide approves of these?)

\def\labelenumi{\arabic{enumi}.}
\def\theenumi{\arabic{enumi}}
\def\labelenumii{(\alph{enumii})}
\def\theenumii{\alph{enumii}}
\def\p@enumii{\theenumi}
\def\labelenumiii{\roman{enumiii}.}
\def\theenumiii{\roman{enumiii}}
\def\p@enumiii{\theenumi(\theenumii)}
\def\labelenumiv{\Alph{enumiv}.}
\def\theenumiv{\Alph{enumiv}}
\def\p@enumiv{\p@enumiii\theenumiii}

% Bullets for itemization. (Who says these are the way to go?)

\def\labelitemi{$\bullet$}      % filled-in disk
\def\labelitemii{\bfseries --}        % en-dash
\def\labelitemiii{$\ast$}       % asterisk
\def\labelitemiv{$\cdot$}       % centered dot

% Parameters for leader dot spacing in toc.
%  They should be defined with \def's.
%    \@pnumwidth : Width of box in which page number is set.
%    \@tocrmarg  : Right margin indentation for all but last line of
%                  multiple-line entries.
%    \@dotsep    : Separation between dots, in mu units.  Should be \def'd to
%                  a number like 2 or 1.7
%

\def\@pnumwidth{1.55em}
\def\@tocrmarg{2.55em}%was 2.55em.
\def\@dotsep{2.5}


%%% Section Headings.
%
% Old pre-1993 \chapter and \part controls sometimes made mistakes in
% the toc because too-long lines were broken in a bad way. Indents
% failed for the second and succeeding lines (James Rieger, Feb 1993).
%  -Problem verified and fixed by GBG Oct 1993
%
% A "groupheader" is one of those phrases that goes all by itself in the TOC.

\def\l@groupheader#1#2{\par\noindent \protect\\ {\bfseries #1}\par}

% A part gets big Roman numerals and a centered title on its own page.

\newcounter{part}

\def\thepart{\Roman{part}}

\def\part{\cleardoublepage
  \thispagestyle{empty}
  \hbox{}\vfil
  \bgroup  \centering \secdef \@part \@spart}

% Parts can be numbered or unnumbered - this handles both.
% Section depth for a part is 1.
%
% Adjust interline spacing and title width
\def\mainheadingwidth{4.5in}% Increased from 4in to 4.8in for Tina Ma's thesis
%
\def\@part[#1]#2{%
  \ifnoisy\testsize{#2}{\HFpart\bfseries}\fi
  \ifnum \c@secnumdepth >-2\relax
      \refstepcounter{part}%
      \addcontentsline{toc}{part}
      {\protect\numberline{PART \thepart}\thesisTOC{#1}}%
      \ifnoisy \typeout{PART \arabic{part}.} \fi
    \else
      \addcontentsline{toc}{part}{\thesisTOC{#1}}
  \fi
  \markboth{}{}
  \parbox{\mainheadingwidth}{%
   \begin{center}
    \ifnum \c@secnumdepth >-2\relax
      \HFpart\bfseries PART \thepart
      \par \vspace{\HFpartSKIP}%
    \fi
    \HFpart\bfseries \uppercase{\nohyphenation #2}
   \end{center}
  }\@endpart
}
% Define \part*{#1}
\def\@spart#1{%
  \ifnoisy\testsize{#1}{\HFpart\bfseries}\fi
  \parbox{\mainheadingwidth}{%
   \begin{center}
    \HFpart\bfseries \uppercase{\nohyphenation #1}
   \end{center}
  }\@endpart
}

% Finish the part page and maybe zip out a blank backside page.

\def\@endpart{\par \egroup \vfil
   \newpage
   \if@twoside
      \hbox{}
      \thispagestyle{empty}
      \newpage
   \fi
}

% Composition of the part's mention in the TOC.
%
\def\l@part{\@PARTdottedtocline{1}{0.0em}{4.6em}}
%
% \@dottedtocline{LEVEL}{INDENT}{NUMWIDTH}{TITLE}{PAGE} :
%   Macro to produce a table of contents line with the following
%   parameters:
%     LEVEL    : If LEVEL > \c@tocdepth, then no line produced.
%     INDENT   : Total indentation from the left margin.
%     NUMWIDTH : Width of box for number if the TITLE has a
%                \numberline command.
%                As of 25 Jan 88, this is also the amount of extra indentation
%                added to second and later lines of a multiple line entry.
%     TITLE    : Contents of entry.
%     PAGE     : Page number.
%
%

\def\@PARTdottedtocline#1#2#3#4#5{\ifnum #1>\c@tocdepth\relax \else
  \addpenalty{-\@highpenalty}%
%  \vskip \z@ plus .2pt
  \vskip 1.25em plus 1pt
  {%
   \leftskip #2\relax \rightskip \@tocrmarg \parfillskip -\rightskip
    \parindent #2\relax\@afterindenttrue
   \interlinepenalty\@M
   \normalsize\bfseries
   \leavevmode
   \@tempdima #3\relax \advance\leftskip \@tempdima \hbox{}\hskip -\leftskip
    #4\nobreak\normalsize\rmfamily
    \leaders\hbox{$\m@th \mkern \@dotsep mu.\mkern \@dotsep
       mu$}\hfill \nobreak \hbox to\@pnumwidth{\hfil\rmfamily #5}\par}\fi}




% "Main headings" are a common generalization throughout the thesis format.
% They are invariably two inches down, centered, and all uppercase boldface.
% (Arithmetic is such that only need 1/2 inch of blank space - invisible
% hrules are kludge to keep the space from vanishing.)
% 9-Feb-1990 D.Eyre. Oct 1993 GBG
%
\def\mainheading#1{
  \hrule width0in \vspace{0.4in} \hrule width0in
  \mainheadingtext{#1}{\HFmainhead\bfseries}%
}
\def\mainheadingtext#1#2{% #1=text, #2=font size and style
  \begin{center}%
    \parbox{\mainheadingwidth}{%
     \begin{center}%
        #2%
        \uppercase{{\nohyphenation #1}}%
      \end{center}%
      \removelastskip
    }% end parbox
   \end{center}%
}

% Chapters are an important case of main headings.

\newcounter{chapter}

\def\thechapter{\arabic{chapter}}

% Turns into ``APPENDIX'' when it's appendix time.

% We need a shorthand to suppress hyphenations in various places
\def\nohyphenation{\hyphenpenalty=10000\exhyphenpenalty=10000\relax}

% Modify uuthesis \chapter macro to take an optional argument

\def\chapter{\@ifnextchar[{\@chapter}{\@@chapter}}

\def\@@chapter#1{\@chapter[#1]{{\nohyphenation #1}}}

\def\@chapapp{CHAPTER}

% Chapter title pages include text but no page number.
% Repair table of contents spacing problems
% Replaced baselineskip commands by special font command. GBG Oct 1993
% Fixed missing \global\@nobreakfalse at end of routine.
% Define \mainheadingwidth=4.0in so adjustments can be made. GBG Oct 1993
%  Affects appendix and chapter headings.
% Fixed hyphen trouble Nov 1992 Nelson Beebe
% Nelson Beebe Nov 1992: Fixes to \uppercase{#1}
% G.B.G 17-April-1992, fixed \vskip to \vspace.
% Try to get 4-inch Thesis Editor rule to work;D.Eyre 9-Feb-1990
%   The actual rule is up to 4.5in heading.
%  Appendix problem: Fix applied to \refstepcounter{chapter}-1990
%
%\chapter command \def's
\def\@optionA{A}
\def\thesisTOC#1{\nohyphenation\linepenalty=0\uppercase{#1}}
%
\def\@chapter[#1]#2{\clearpage
  \thispagestyle{empty}%
  \global \@topnum \z@
  \if@oneappendix% Different treatment for one appendix
           \ifnoisy \typeout{\@chapapp:} \fi
           \refstepcounter{chapter}
           \addcontentsline{toc}{chapter}{\thesisTOC{\@chapapp: #1}}%
           \mainheading{\@chapapp}%
  \else
    \ifnum \c@secnumdepth >\m@ne\relax
      \refstepcounter{chapter}
% If chapter==1 or chapter==A, then TOC gets header
      \edef\@tmp{\thechapter}
      \ifx\@tmp\@optionA
        \addcontentsline{toc}{groupheader}{APPENDICES}%
      \fi
      \ifx\@tmp\@optionONE
        \addcontentsline{toc}{groupheader}{CHAPTERS}
      \fi
        \addcontentsline{toc}{chapter}%
        {\protect\numberline{\thechapter.}\thesisTOC{#1}}%
        \addcontentsline{toc}{groupheader}{\protect\vspace{-22pt}}
    \else
       \addcontentsline{toc}{chapter}{\thesisTOC{#1}}%
    \fi
    \ifnoisy \typeout{\@chapapp\space\thechapter.} \fi
    \mainheading{\@chapapp\space\thechapter}
  \fi
  \mainheadingtext{#2}{\HFchapter\bfseries}% Text=boldface chapter header font
  \vskip \HFchapterSKIP
  \nobreak\@afterheading\global\@nobreakfalse
}
%%
%% Use \minusline, \plusline to remove or insert whole lines.
%% Use \minushalfline, \plushalfline to remove or insert half lines.
%%
\def\minusline{\vspace{-\normalbaselineskip}}
\def\plusline{\vspace{\normalbaselineskip}}
\def\minushalfline{\vspace{-0.5\normalbaselineskip}}
\def\plushalfline{\vspace{0.5\normalbaselineskip}}
\def\minusfourthline{\vspace{-0.25\normalbaselineskip}}
\def\plusfourthline{\vspace{0.25\normalbaselineskip}}
%
% 1988: In the TOC, the chapter number and title appear in bf, followed
%       by leaders then right-justified page number.
%       Original 1.0em changed to .75em and changed #1\hfil to #1\mydotfill.
%
\def\l@chapter{\@CHAPTERdottedtocline{2}{0.0em}{2em}}
%
\def\@CHAPTERdottedtocline#1#2#3#4#5{\pagebreak[3]
  \ifnum #1>\c@tocdepth\relax \else
%  \vskip \z@ plus .2pt
  \vskip .75em plus 1pt
  {%
    \leftskip #2\relax \rightskip \@tocrmarg \parfillskip -\rightskip
    \parindent #2\relax\@afterindenttrue
   \interlinepenalty\@M
   \leavevmode
   \@tempdima #3\relax \advance\leftskip \@tempdima \hbox{}\hskip -\leftskip
   \normalsize\bfseries
    #4\nobreak\normalsize\rmfamily
    \leaders\hbox{$\m@th \mkern \@dotsep mu.\mkern \@dotsep
       mu$}\hfill\nobreak \hbox to\@pnumwidth{\hfil\rmfamily #5}\par}\fi}


% Section numbering/titling is a problem for UU theses, since the format is
% different, depending on how many subsections one has.  Therefore, we have
% some commands to set the number of levels.

% \@Ustartsection {NAME}{LEVEL}{INDENT}{BEFORESKIP}{AFTERSKIP}{STYLE}
%            optional * [ALTHEADING]{HEADING}
%    Generic command to start a section.
%    NAME       : e.g., 'subsection'
%    LEVEL      : a number, denoting depth of section -- e.g., chapter=1,
%                 section = 2, etc.  A section number will be printed if
%                 and only if LEVEL < or = the value of the secnumdepth
%                 counter.
%    INDENT     : Indentation of heading from left margin
%    BEFORESKIP : Absolute value = skip to leave above the heading.
%                 If negative, then paragraph indent of text following
%                 heading is suppressed.
%    AFTERSKIP  : if positive, then skip to leave below heading,
%                       else - skip to leave to right of run-in heading.
%    STYLE      : commands to set style
%  If '*' missing, then increments the counter.  If it is present, then
%  there should be no [ALTHEADING] argument.  A sectioning command
%  is normally defined to \@Ustartsection + its first six arguments.
%
\def\@Ustartsection#1#2#3#4#5#6{%
   \@testsizetrue\ifnum\c@secnumdepth=#2\relax \@testsizefalse \fi
   \if@noskipsec \leavevmode \fi
   \par\@tempskipa #4\relax
   \@afterindenttrue
   \ifdim \@tempskipa <\z@ \@tempskipa -\@tempskipa \@afterindentfalse\fi
   \if@nobreak \everypar{}
     \addpenalty{\@secpenalty}% Added by GBG
     \advance\@tempskipa by -14.5pt% Added by GBG
     \addvspace{\@tempskipa}% Added by GBG
   \else
     \addpenalty{\@secpenalty}\addvspace{\@tempskipa}\fi
     \@ifstar
     {\@Ussect{#3}{#4}{#5}{#6}}{\@dblarg{\@Usect{#1}{#2}{#3}{#4}{#5}{#6}}}}
%
% \testsize checks for title width and issues error message if too wide
%
% Build a tabular array from the title string. Measure the width of the
% table and test it against the maximum allowed width of 4.0 inches (about
% 1/2 inch is used for the section number). Issue error message if the
% width is exceeded but otherwise do nothing. -GBG OCT 1993
%
\newdimen\@titlesize
\newdimen\@maxsize
\newif\if@testsize \@testsizetrue
\def\testsize#1#2{%
\if@testsize
\let\@CENTER\center
\let\center\relax
\@maxsize=\mainheadingwidth
\advance\@maxsize by -0.5in
\setbox\@tempboxa=\hbox{{#2\begin{tabular}{@{}l@{}}#1\end{tabular}}}%
\@titlesize\wd\@tempboxa\relax
\ifdim\@titlesize>\@maxsize
\ifnoisy
\advance\@titlesize by -\@maxsize
\typeout{Title Error: (\the\@titlesize\space too wide)}
\typeout{Page \thepage, Title="#1"}\typeout{}
\fi
\fi
\let\center\@CENTER
\fi
}

\def\@Usect#1#2#3#4#5#6[#7]#8{%
\testsize{#8}{#6}\@testsizetrue
\ifnum #2>\c@secnumdepth\relax
     \def\@svsec{}\else
     \refstepcounter{#1}\edef\@svsec{\csname the#1\endcsname\hskip 1em }\fi
     \@tempskipa #5\relax
     \ifdim \@tempskipa>\z@
        \begingroup #6\relax
          \@hangfrom{\hskip #3\relax\@svsec}{\interlinepenalty \@M #8\par}%
        \endgroup
        \csname #1mark\endcsname{#7}%
        \ifnum\c@tocdepth>\z@\relax
          \addcontentsline{toc}{#1}%
           {\ifnum #2>\c@secnumdepth \else
            \protect\numberline{\csname the#1\endcsname}\fi #7}%
        \fi
     \else
        \def\@svsechd{#6\hskip #3\relax  %% \relax added 2 May 90
            \@svsec #8\csname #1mark\endcsname {#7}%
            \ifnum\c@tocdepth>\z@\relax
               \addcontentsline{toc}{#1}%
                {\ifnum #2>\c@secnumdepth \else
                 \protect\numberline{\csname the#1\endcsname}\fi #7}%
            \fi
        }
     \fi
     \@xsect{#5}%
}
%

\def\@Ussect#1#2#3#4#5{\testsize{#5}{#4}\@testsizetrue
\@tempskipa #3\relax
  \ifdim \@tempskipa>\z@
     \begingroup #4%
     \@hangfrom{\hskip #1}{\interlinepenalty \@M%
     \nopagebreak[4]%
     #5\par}%
     \endgroup
   \else \def\@svsechd{#4\hskip #1\relax #5}\fi
    \@xsect{#3}}


% Three levels of heading is most typical.

% Put in fixed spacing above and below headings, rather than keying it
% to \normalbaselineskip, which depends on the size of the current
% font.
%
% UofU thesis guidelines require two blank lines before some headings,
% and one blank line after after.  Paragraph headings get indent equal
% to the usual indent for a paragraph. The heading ends with a period.
% Then the text of the paragraph continues on the same line as the
% heading. If two displayed subheadings follow one another without any
% text, then they are separated by only one blank line.
%
% I translate these requirements into typography as normal font height =
% 14.5pt, large font height = 18pt. These are the only two font sizes
% used for headings in this style, hence everything can be set from
% them.
%
% Fixes by Nelson Beebe 22 June 1993 and GBG Oct 1993
%
\def\threelevels{
% The value of the counter secnumdepth gives the depth of the
% highest-level sectioning command that is to produce section numbers.
  \setcounter{secnumdepth}{5}
  \setcounter{tocdepth}{5}
  \newcounter{section}[chapter]
  \newcounter{subsection}[section]
  \newcounter{paragraph}[subsection]
  \def\thesection{\thechapter.\arabic{section}}
  \def\thesubsection{\thesection.\arabic{subsection}}
  \def\theparagraph{\thesubsection.\arabic{paragraph}}
  \def\section{\@Ustartsection{section}{3}
                        {\z@}   %indent. \z@ == 0pt
             {\HFsectionSKIP}   % before=\baselineskip for \large
                        {1sp}   %after
        {\center                %style
         \HFsection\bfseries}}
  \def\subsection{\@Ustartsection{subsection}{4}
                        {\z@}   %indent
          {\HFsubsectionSKIP}   % before=\baselineskip for \normalsize
                        {1sp}   %after
          {\HFsubsection\bfseries}}   %style
  \def\subsubsection{\typeout{No subsubsections allowed with three levels.}}
  \def\paragraph{\@Ustartsection{paragraph}{5}
                 {\parindent}   %indent
                        {1sp}   %before
                       {-1em}   %after
           {\HFparagraph\bfseries}}   %style
  \def\subsubsubsection{\typeout{No subsubsubsections allowed in threelevels.}}
  \def\subparagraph{\typeout{No subparagraphs allowed in this style.}}
  \def\l@section{\@dottedtocline{3}{2em}{2em}}
  \def\l@subsection{\@dottedtocline{4}{3.5em}{2.8em}}
  \def\l@paragraph{\@dottedtocline{5}{5em}{3.6em}}
}


% Four levels of sectioning
% Fixes by Nelson Beebe Nov 1992. Fixes by GBG Oct 1993
\def\fourlevels{
  \threelevels
  \setcounter{secnumdepth}{6}
  \setcounter{tocdepth}{6}
  \newcounter{subsubsection}[subsection]
  \def\thesubsubsection{\thesubsection.\arabic{subsubsection}}
  \def\theparagraph{\thesubsubsection.\arabic{paragraph}}
  \def\subsection{\@Ustartsection{subsection}{4}% Fix Nov 1992
                        {\z@}   %indent. \z@ == 0pt
          {\HFsubsectionSKIP}   %before
                        {1sp}   %after
   {\center                     %style
    \HFsubsection\bfseries}}
  \def\subsubsection{\@Ustartsection{subsubsection}{5}% Fix Nov 1992
                        {\z@}   %indent
       {\HFsubsubsectionSKIP}   %before
                        {1sp}   %after
  {\HFsubsubsection\bfseries          %style
   \hsize=\mainheadingwidth}}
  \def\paragraph{\@Ustartsection{paragraph}{6}
                 {\parindent}   %indent
                        {1sp}   %before
                       {-1em}   %after
           {\HFparagraph\bfseries}}   %style
  \def\subsubsubsection{\typeout{No subsubsubsections allowed in fourlevels.}}
  \def\subparagraph{\typeout{No subparagraphs allowed in this style.}}
  \def\l@subsubsection{\@dottedtocline{5}{5em}{3.6em}}
  \def\l@paragraph{\@dottedtocline{6}{6.5em}{4.2em}}
}

% Five levels of sectioning
%  by GBG Oct 1993
\def\fivelevels{
  \fourlevels
  \setcounter{secnumdepth}{7}
  \setcounter{tocdepth}{7}
  \newcounter{subsubsubsection}[subsubsection]
  \def\thesubsubsubsection{\thesubsubsection.\arabic{subsubsubsection}}
  \def\theparagraph{\thesubsubsubsection.\arabic{paragraph}}
  \def\subsubsubsection{\@Ustartsection{subsubsubsection}{6}
                        {\z@}   %indent
       {\HFsubsubsectionSKIP}   %before
                        {1sp}   %after
  {\HFsubsubsection\bfseries          %style
   \hsize=\mainheadingwidth}}
  \def\paragraph{\@Ustartsection{paragraph}{7}
                 {\parindent}   %indent
                        {1sp}   %before
                       {-1em}   %after
           {\HFparagraph\bfseries}}   %style
  \def\l@subsubsubsection{\@dottedtocline{6}{6.5em}{4.2em}}
  \def\l@paragraph{\@dottedtocline{7}{8em}{5.0em}}
}


%%% Footnote parameters
%
%   \skip\footins   : Space between main text and footnotes.  The rule
%                     separating footnotes from text occurs in this space.
%                     This space lies above the strut of height \footnotesep
%                     which is at the beginning of the first footnote.
\skip\footins 1em plus 4pt minus 2pt
\skip\@mpfootins = \skip\footins

% 1.8 inches is more or less the 20 spaces that the footnote rule should be.
\def\footnoterule{\kern-3\p@ \hrule width 1.8in \kern 2.6\p@}

% Footnotes numbered on a per-chapter basis (Chicago Manual, 2.22).
\@addtoreset{footnote}{chapter}

% Wrapper for footnote text, does paragraph indent and adds footnote number.
\long\def\@makefntext#1{\parindent 1em \noindent
  \hbox to 1.8em {\hss$^{\@thefnmark}$}#1
}

% Height of strut placed at the beginning of every footnote is enough to
% make something like doublespacing, but will apply to all footnotes, not
% just the first (sigh).
\footnotesep 20pt

%%
%%  End of footnote options
%%

%    ****************************************
%    *         FIGURES AND TABLES           *
%    ****************************************
%
% Fix by Nelson Beebe Nov 1992 and Sep 1993
% More fixes by GBG Oct 1993
%
% Figures may appear physically after a reference at the bottom of
% the same page or the top of the next page. Figures deferred to the
% next page must all appear at once, producing a page of floats if necessary.
%
% To avoid a single figure centered on a float page the [t] option can be
% applied to the last figure in a list. This is done experimentally,
% because the figure placement depends on text also.
%
%% FIGURE CAPTION
% The actual creator of table/figure captions takes a number and the caption
% text.  This macro should be used with specific line breaks
% "\protect\\" because the title is to be in inverted pyramid format if
% the title extends beyond the figure (or table). A parbox can be used
% to force the text to a certain width (the figure width, for example).
% Often the text is already set in an inverted pyramid and nothing else
% needs to be done. Suggestion: try it first, then fix as needed.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No longer used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% \newif\ifoldstylecaption \oldstylecaptiontrue
% \def\captionstyleparagraph{\oldstylecaptiontrue}
% \def\captionstylecentered{\oldstylecaptionfalse}
% \def\captionlineskip{13.5pt}% Single spaced captions
%
% \long\def\@makecaption#1#2{
% \ifoldstylecaption
% The old one from 1988:
%   \vskip 10pt
%   \baselineskip=\captionlineskip
%   \setbox\@tempboxa\hbox{{#1}. #2}%
%   \ifdim \wd\@tempboxa >\hsize
%       \unhbox\@tempboxa\par
%     \else
%       \hbox to\hsize{\hfil\box\@tempboxa\hfil}%
%   \fi
% \else
% % Inverted pyramid centered in width of figure is desired.
%  \begin{center}
%  \baselineskip=\captionlineskip
%  {#1}. #2 \end{center}%
% \fi
% }
%
% \newif\ifcaptionON \captionONtrue
%
% \def\caption{
%   \ifcaptionON\refstepcounter\@captype\fi \@dblarg{\@caption\@captype}}
%
% \long\def\@caption#1[#2]#3{\par
%   \ifcaptionON\addcontentsline{\csname
%   ext@#1\endcsname}{#1}{\protect\numberline{\csname
%   the#1\endcsname}{\ignorespaces #2}}\fi\begingroup
%     \@parboxrestore
%     \normalsize
%     \@makecaption{\csname fnum@#1\endcsname}{\ignorespaces #3}\par
%   \endgroup}
%

%% Figure and Table Style Parameters.

% The thesis format dictates that figures "should be placed so they look
% balanced on the page".
% Defined limits: No more than 6.5 inches (about 3/4) of page.  Pages of
% floats are not constrained too much, they can have many many figures
% on them.  Everything is separated by triple spacing (2em), plus just
% enough stretchiness to make things work.
%
% The total number of floats on a page is rather large because the thesis
% editor permits any number of floats if they're on their own page.
%
% Float placement parameters.  See LaTeX manual for their definition.
% This fix by GBG Oct-1993 to set defaults for figure placement.
%
% Instructions: Use only [b], [t] or [p] in a figure or table
%               placement. The default is [p], a page of floats.
%               The spec [h] will work, if that is the only way
%               that the figure can be placed where it is wanted.
%
%               The figure environment uses the entire \textwidth
%               to place a figure. However, the environment may contain
%               multiple figures and possibly text or tables.
%
% Police work for figures:
%
%       1. A figure must follow its reference, at the bottom of the
%          current page or immediately on the top of the next page.
%          Start with [b] and change to [t] if it doesn't fit.
%
%       2. An option of [t], which causes proper figure placement after
%          the reference on the preceding page, can be followed by a
%          figure with [b] option. Both figures can be used on the same
%          page. Warning:  usually three won't fit.
%
%       3. Three figures in a row should generate one figure at the
%          bottom and two figures on the next page. Often the text looks
%          best if the two deferred figures are on a page of floats by
%          themselves. Use [p] to make it happen. The decision is based
%          on 6 lines of text. If six lines won't fit on the same page
%          as the two figures (at the bottom), then make a page of
%          floats for the two figures.
%
%        4. If pages of floats have a widow figure at the end, then it is
%         to appear as a [t] option on the next page of text.
%
\setcounter{topnumber}{10}
\def\topfraction{1}
\setcounter{bottomnumber}{10}
\def\bottomfraction{1}
\setcounter{totalnumber}{20}
\def\textfraction{0}
\def\floatpagefraction{.5}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The double-column float parameters aren't used, but things fail without them.
\setcounter{dbltopnumber}{10}
\def\dbltopfraction{1}
\def\dblfloatpagefraction{.5}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\floatsep 2em plus 2pt minus 4pt        % vertical space between floats
\textfloatsep 2em plus 2pt minus 4pt    % vertical between text and floats
\intextsep 2em plus 4pt minus 4pt       % vertical around mid-page float (?)
\dblfloatsep
\dbltextfloatsep

% Stretch at top of float page/column. (Must be 0pt plus ...)
\@fptop 0pt plus 1fil
\@fpsep 10pt plus 2fil
\@fpbot 0pt plus 1fil

% Figures are numbered by chapters.
\newcounter{figure}[chapter]

%
% Default float placement option should be a page of floats
% according to the thesis office. This means that deferred figures
% that add up to \floatpagefraction * \pageheight or more will generate
% a full page of floats containing no text.
%
\def\fps@figure{p}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\thefigure{\thechapter.\@arabic\c@figure}
\def\fnum@figure{\bfseries Figure \thefigure}

% Define figures as a type of floating object (see report.doc for details).
% Figure environment definition proper.
%
\def\figure{\@float{figure}}
\let\endfigure\end@float

% Table definitions nearly identical to those for figures, since
% they are treated the same internally, except for list of tables
% and list of figures generation. In particular, tables and figures
% are output in the same interleaved order as received.
%
\newcounter{table}[chapter]
\def\fps@table{p}
\def\ftype@table{2}
\def\ext@table{lot}
\def\thetable{\thechapter.\@arabic\c@table}
\def\fnum@table{\bfseries Table \thetable}
\def\table{\@float{table}}
\let\endtable\end@float

%%% Page style commands.

% Initialize all of TeX's marks.
%
\mark{{}{}}

% Page style is pretty simple, only need a number in the corner.
%
\if@twoside
  \def\ps@headings{
    \let\@mkboth\markboth
    \def\@oddfoot{}
    \def\@evenfoot{}
    \def\@evenhead{\rmfamily \thepage\hfil \slshape \leftmark}
    \def\@oddhead{\hbox{}\slshape \rightmark \hfil \rmfamily \thepage}
  }
\else
  \def\ps@headings{
    \let\@mkboth\markboth
    \def\@oddfoot{}
    \def\@evenfoot{}
    \def\@oddhead{\hbox {}\slshape \rightmark \hfil \rmfamily \thepage}
  }
\fi

% The 'myheadings' page style is used here and there.
%
\def\ps@myheadings{
  \let\@mkboth\@gobbletwo
  \def\@oddhead{\hbox{}\slshape\rightmark \hfil \rmfamily \thepage}
  \def\@oddfoot{}
  \def\@evenhead{\rmfamily \thepage \hfil \slshape \leftmark \hbox{}}
  \def\@evenfoot{}
}

% Definitions of standard environments.
%
\def\verse{
  \let\\=\@centercr
  \list{}{\itemsep\z@ \itemindent -1.5em\listparindent \itemindent
          \rightmargin\leftmargin\advance\leftmargin 1.5em}
  \item[]
}
\let\endverse=\endlist

\def\quotation{
  \list{}{\listparindent 1.5em
          \itemindent\listparindent
          \rightmargin\leftmargin\parsep 0pt plus 1pt}
  \item[]
}
\let\endquotation=\endlist

\def\quote{
  \list{}{\rightmargin\leftmargin}
  \item[]
}
\let\endquote=\endlist

\def\descriptionlabel#1{\hspace\labelsep \bfseries #1}

\def\description{
  \list{}{\labelwidth\z@ \itemindent-\leftmargin
          \let\makelabel\descriptionlabel}
}
\let\enddescription=\endlist

\newdimen\descriptionmargin

\descriptionmargin=3em

\def\normalspace{\if@report \singlespace \else \doublespace \fi}

\def\endsinglespace{\normalspace}

\def\enddoublespace{\normalspace}

\def\epigraph{
  \singlespace
  \list{}{\listparindent 1.5em
          \itemindent \listparindent
          \rightmargin \leftmargin \parsep 0pt plus 1pt}
  \item[]
}

\def\endepigraph{\endlist\endsinglespace}

%%% The following commands are used in defining the pieces of a thesis.
%%% They are defined in the same order that they will appear in the work.

\def\frontmatter#1#2#3{
  \requiredfrontmatter{#1}
  \listoftables
  \listoffigures
  \preface{#2}{#3}
}

% The first part of the front matter always looks the same.

\def\requiredfrontmatter#1{
  \frontmatterformat
  \if@honors \@reporttrue\fi
  \if@csreport \@reportrue\fi
  \if@report \reporttitlepage \else \titlepage \fi
  \copyrightpage
  \if@report \else \committeeapproval \fi
  \if@report \else \readingapproval \fi
% by GBG 17-April-1992
\ifnoisy
\typeout{Abstract: Under 350 words, no citations.}
\fi
  \preface{#1}{ABSTRACT}
  \ifdedication \dedicationpage \fi
  \tableofcontents
}

% Do whatever is needed to set up the basic page format of the front matter.
%
\def\frontmatterformat{%
  \pagestyle{plain}%
  \pagenumbering{roman}%
  \global\c@page = -1
  \parskip 0pt plus 1pt
}

% The format of the title page is laid out in elaborate detail.
% Fixes by Nelson Beebe Nov 1992:
% The title is centered in a vbox of fixed height so that it can
% have up to 9 lines.  Version 1.3 and earlier used explicit blank lines,
% which caused spill of the entire titlepage onto a second
% numbered page if a multiline title was specified.
%
\newlength{\minilength}
\setlength{\minilength}{.95\textwidth}
%
\def\titlepage{
  \thispagestyle{empty}%
  \noindent\hspace{1em}%
  \begin{minipage}[c]{\minilength}%
    \begin{center}%
      \vbox to 9\baselineskip {%
          \HFmainhead\bfseries
          \vss
          \begin{center}
          \@title
          \end{center}
          \vss
          }%
      \mbox{} \\
    by\\
      \mbox{} \\
      \@author \\
      \vbox to 15\baselineskip
      {
        \vss
        A \@thesistype ~submitted to the faculty of\\
        The University of Utah\\
        in partial fulfillment of the requirements for the degree of
        \vss
      }
      \@degree \\
      \vskip 7\baselineskip
      \expandafter{\@department}\\
      \mbox{} \\
      The University of Utah \\
      \mbox{} \\
      \@submitdate%
    \end{center}%
  \end{minipage}%
  \ifnoisy \typeout{Title.} \fi
  \newpage% Added 1995 by GBG; see Eric Eide notes.
}

% The copyright page is fairly simple, but we have to make sure it is not
% counted in the page numbering.

\def\copyrightpage{\newpage
  \thispagestyle{empty}%
  \rule{0in}{4.25in}%
  \noindent\begin{minipage}{.95\textwidth}%
    \begin{center}%
      Copyright \copyright\ \@author\ \@copyrightyear \\
      \mbox{}  \\
      All Rights Reserved
    \end{center}%
  \end{minipage}%
  \rule{0in}{4.2in}%
  \if@report \else \addtocounter{page}{-1}\fi
  \ifnoisy \typeout{Copyright.} \fi
  \newpage% Added 1995 by GBG; see Eric Eide notes.
}

% The committee approval page is highly constrained.  Its format is not
% exactly identical to the preprinted ones.
%

\long\def\ReaderPerson#1{%
    \mbox{}\\  \mbox{}\\ \mbox{} \\
    \mbox{\rule{1.5in}{0.4pt}\rule{0.75in}{0in}\rule{3.375in}{0.4pt}}\\
    \rule{2.375in}{0in}%
    \phantom{{\HFapprovalsmall Chair:\hspace{1.5em}}}#1\\}

\long\def\committeeapproval{\newpage\HFapproval
  \thispagestyle{empty}
  \begin{center}
    \mbox{} \\ \mbox{} \\
    {\scshape THE UNIVERSITY OF UTAH GRADUATE SCHOOL} \\
    \mbox{} \\ \mbox{} \\ \mbox{} \\
    {\Large\bfseries SUPERVISORY COMMITTEE APPROVAL} \\
    \mbox{} \\ \mbox{} \\ \mbox{} \\
     of a \@thesistype{} submitted by \\
    \mbox{} \\ \mbox{} \\ \mbox{} \\
    \@author \\
    \mbox{} \\ \mbox{} \\ \mbox{} \\
  \end{center}
  This \@thesistype{} has been read
  by each member of the following supervisory
  committee and by majority vote has been found to be satisfactory.
  \begin{flushleft}
    \mbox{}\\ \mbox{}\\  \mbox{}\\ \mbox{} \\
    \mbox{\rule{1.5in}{0.4pt}\rule{0.75in}{0in}%
    \rule{3.375in}{0.4pt}}\\
    \rule{2.375in}{0in}{\HFapprovalsmall Chair:\hspace{1em}}
    \@committeechair\\
    \ReaderPerson{\@firstreader}
    \ReaderPerson{\@secondreader}
    \iffourth \ReaderPerson{\@thirdreader} \fi
    \iffifth \ReaderPerson{\@fourthreader} \fi
  \end{flushleft}
  \ifnoisy \typeout{Committee Approval.} \fi
  \normalsize\rmfamily\normalspace
}

% The reading approval is page iii, and its format is also highly constrained.

\long\def\readingapproval{\newpage
  \thispagestyle{empty}\HFapproval
  \begin{center}
    \mbox{} \\ \mbox{} \\
    {\scshape THE UNIVERSITY OF UTAH GRADUATE SCHOOL} \\
    \mbox{} \\ \mbox{} \\ \mbox{} \\
    {\Large\bfseries FINAL READING APPROVAL}
    \mbox{} \\ \mbox{} \\
    \mbox{} \\ \mbox{} \\ \mbox{} \\
    \begin{minipage}{\textwidth}
      To the Graduate Council of the University of Utah:\\
      \mbox{} \\ \mbox{} \\
      I have read the \@thesistype{} of
      \underline{\makebox[2.62in]{\HFapproval \@author}}
      in its final form and have found that (1) its format, citations,
      and bibliographic style are consistent and acceptable; (2) its
      illustrative materials including figures, tables, and charts are in
      place; and (3) the final manuscript is satisfactory to the Supervisory
      Committee and is ready for submission to The Graduate School.\\[0.25in]
      \mbox{\rule{0.33\linewidth}{0.4pt}\rule{0.066\linewidth}{0in}
      \rule{0.58\linewidth}{0.4pt}} \\
      {\HFapprovalsmall Date}\rule{0.35\linewidth}{0in}
      \@committeechair
      \\ \rule{0.41\linewidth}{0in}
      {\HFapprovalsmall Chair, Supervisory Committee} \\
      \mbox{}
    \end{minipage} \\[0.9in]
    Approved for the Major Department \\[0.5in]
    \rule{.45\textwidth}{0.4pt} \\
    \@departmentchair \\ {\HFapprovalsmall Chair/Dean} \\[0.45in]
     Approved for the Graduate Council \\[0.5in]
    \rule{.45\textwidth}{0.4pt} \\
    \@graduatedean \\ {\HFapprovalsmall Dean of The Graduate School} \\
  \end{center}
  \ifnoisy \typeout{Reading Approval.} \fi
  \normalsize\rmfamily\normalspace
}

% Optional dedication page (thanx Eric!)

\def\dedicationpage{\newpage
  \thispagestyle{empty}
  \vspace*{3.0in}
  \begin{center}{\@dedication}\end{center}
  \ifnoisy \typeout{Dedication.} \fi
}

% The first table of contents page is headed with the centered word "CONTENTS".
% No page number on it.

\def\tableofcontents{\newpage
  \thispagestyle{empty}
  \mainheading{CONTENTS}
  \par\removelastskip\singlespace\par\removelastskip% GBG Oct 1993
  \vskip 5pt
  \@starttoc{toc}
  \ifnoisy \typeout{Contents.} \fi
}

% Make lists of figures
% Entries are single-spaced but between entries there is a double-space.
%
\def\listoffigures{
  \iflistoffigures
    \newpage
    \thispagestyle{empty}
    \addcontentsline{toc}{chapter}{LIST OF FIGURES}
    \mainheading{LIST OF FIGURES}\vskip 3pt
    \ifnoisy \typeout{List of Figures.} \fi
    \@starttoc{lof}
  \else
    \ifnoisy \typeout{No list of figures.} \fi
  \fi
}

% Doublespace between each entry in the list and use some leader dots.
% This sometimes breaks for too many figures: number runs into text.
% Fixed by using a larger width "3.3em" and different justification.

% Fill number box and right justify (old one is left justify)
\def\uunumberline#1{\hbox to\@tempdima{\hfil #1\hspace{6pt}}}

\def\l@figure#1#2{\let\numberline\uunumberline\vskip 6pt \@dottedtocline{1}{0.0em}{3.3em}{#1}{#2}}

% Tables are just like figures, but the list is on a different page.

\def\listoftables{
  \iflistoftables
    \newpage
    \thispagestyle{empty}
    \addcontentsline{toc}{chapter}{LIST OF TABLES}
    \mainheading{LIST OF TABLES}\vskip 3pt
    \ifnoisy \typeout{List of Tables.} \fi
    \@starttoc{lot}
  \else
    \ifnoisy\typeout{No list of tables.}\fi
  \fi
}

% Entries in list handled exactly as for figures.

\let\l@table\l@figure

% For use in making special front matter sections.

\def\optionalfront#1{\newpage
  \thispagestyle{empty}
  \addcontentsline{toc}{chapter}{\uppercase{#1}}
  \mainheading{#1}\fixmainheadingSKIP
  \ifnoisy \typeout{#1.} \fi
}

% Acknowledgment page is also treated like a chapter.

\def\preface#1#2{
  \prefacesection{#2}
  \par \input{#1}
  \ifnoisy \typeout{#2.} \fi
}

% Optional preface is treated almost like a chapter - goes in toc, but no
% number assigned.

\def\prefacesection#1{\newpage
  \thispagestyle{empty}
  \normalspace
  \addcontentsline{toc}{chapter}{\uppercase{#1}}
  \mainheading{#1}\fixmainheadingSKIP
}

% Once all the front matter has been processed, we need to switch into
% arabic page numbering for the remainder of the thesis.  This command
% should also be usable by itself, which saves 10+ pages on drafts.

\def\maintext{
  \newpage
  \normalspace
  \pagestyle{myheadings}
  \markright{}
  \pagenumbering{arabic}
  \parskip 0pt plus 1pt minus 1pt
}

% The \appendix command must reset counters, redefine \@chapapp, and
% redefine \thechapter so as to get different titles and headings.
% It is a switch rather than a heading command; \chapter marks out appendices
% still.
%
% Put colon instead of period after Appendix in TOC; the thesis
% office requested this in Nov 1992. Fix by Nelson Beebe.
% Various past fixes by David Eyre, Grant Gustafson, Nelson Beebe 1988-1992.

\newcounter{oldchapter}
\newcounter{oldtocdepth}
\newif\if@oneappendix \@oneappendixfalse

\def\appendix{\newpage
  \ifnum\numberofappendices=0\relax
    \typeout{Error: 'numberofappendices' was 0, now set to 2.}%
    \typeout{Change it to 1 or higher so TOC is correct.}
    \numberofappendices=2
  \fi
  \ifnum \numberofappendices=1\relax
    \@oneappendixtrue
  \else
    \@oneappendixfalse
  \fi
  \setcounter{oldchapter}{\value{chapter}}
  \setcounter{chapter}{0}%
  \setcounter{section}{0}%
  \setcounter{oldtocdepth}{\value{tocdepth}}
  \setcounter{tocdepth}{0}% Only chapters go into appendix TOC
  \def\thechapter{\Alph{chapter}}%
  \def\@chapapp{APPENDIX}%
}

\def\noappendix{
\@oneappendixfalse
\def\@chapapp{CHAPTER}
\def\thechapter{\arabic{chapter}}
\setcounter{chapter}{\value{oldchapter}}
\setcounter{section}{0}
\setcounter{tocdepth}{\value{oldtocdepth}}
}

% Whether or not bibtex is used, \thebibliography eventually shows up
% (in .bbl file if bibtex used).  This command makes the heading,
% defines a blocking command, sets formatting to be loose, and doesn't
% let periods cause extra space to be inserted.

\newif\ifrawbibliography \rawbibliographyfalse
\def\thebibliography#1{%
 \ifrawbibliography
 \else
  \newpage
  \thispagestyle{empty}%
  \addcontentsline{toc}{chapter}{REFERENCES}%
% Switch singlespace to after the heading gets printed.
  \mainheading{REFERENCES}%
  \par\removelastskip\singlespace\par\removelastskip% GBG Oct 1993
  \fixmainheadingSKIP
 \fi
  \list{[\arabic{enumi}]}%
  {\settowidth\labelwidth{[#1]}\leftmargin%
     \labelwidth\advance\leftmargin\labelsep\usecounter{enumi}}%
  \def\newblock{\hskip .11em plus .33em minus -.07em}%
  \sloppy\clubpenalty4000\widowpenalty4000%
  \sfcode`\.=1000\relax
}

\let\endthebibliography=\endlist% why not \endsinglespace?

% The thesis format *does* allow for an index, believe it or not.

\def\theindex{\newpage
  \columnseprule \z@
  \columnsep 35pt\relax
  \twocolumn[\@makeschapterhead{Index}]
  \@mkboth{INDEX}{INDEX}
  \thispagestyle{plain}
  \parindent \z@
  \parskip\z@ plus .3pt\relax
  \let\item\@idxitem
}

\def\@idxitem{\par \hangindent 40pt}

\def\subitem{\par \hangindent 40pt \hspace*{20pt}}

\def\subsubitem{\par \hangindent 40pt \hspace*{30pt}}

\def\endtheindex{\if@restonecol\onecolumn\else\clearpage\fi}

\def\indexspace{\par \vskip 10pt plus 5pt minus 3pt\relax}

% The thesis format allows a vita at the end, but we leave all the details
% up to the author.

\def\vita{\newpage
  \thispagestyle{empty}
  \addcontentsline{toc}{chapter}{VITA}
  \mainheading{VITA}\fixmainheadingSKIP
}

% The abstract title page is a variation on the regular title page.

\def\abstracttitlepage{
  \thispagestyle{empty}
  \begin{minipage}[c]{\minilength}
    \begin{center}
      \mbox{} \\ \mbox{} \\ \mbox{}\\ \mbox{} \\ \mbox{} \\
      {\HFmainhead\bfseries \@title}
      \mbox{} \\ \mbox{} \\ \mbox{} \\
      by\\
      \mbox{} \\
      \@author \\
      \mbox{} \\ \mbox{} \\ \mbox{} \\ \mbox{} \\ \mbox{} \\
      An abstract of a  \@thesistype ~submitted to the faculty of\\
      The University of Utah\\
      in partial fulfillment of the requirements for the degree of\\
      \mbox{}\\  \mbox{} \\ \mbox{}\\
      \@degree \\
      \mbox{}\\ \mbox{} \\ \mbox{}\\ \mbox{} \\ \mbox{} \\ \mbox{}\\
      \@committeechair ~~~~~~~~~Chair, Supervisory Committee \\
      \mbox{}\\
      \@chairtitle
      \mbox{}\\ \mbox{} \\ \mbox{} \\
      \expandafter{\@department}\\
      \mbox{} \\
      The University of Utah \\
      \mbox{} \\
      \@submitdate
    \end{center}
  \end{minipage}
  \ifnoisy \typeout{Abstract Title.} \fi
  \newpage% Added 1995 by GBG; see Eric Eide notes.
}

% Titlepage for a tech report.  This is what the CS dept likes - redefine
% this command to get something appropriate for other departments.
% Should use option [csreport] at top of thesis source file.
%
\long\def\CSREPORTTITLE{
  \unitlength = 1in
  \begin{picture}(5.80,8.6)(0,0)
    \put(.2,6.0){\framebox(4.2,2.5){%
      \begin{minipage}[c]{4.19in}%
        \begin{center}
          \sf \@title \\*[.5cm]
          \sf \@author \\*[.5cm]
          \tt Tech. Report. \@reportnumber
        \end{center}%
      \end{minipage}}}%
  \end{picture}%
  \unitlength=1pt
\ifnoisy \typeout{CS Technical Report Title.} \fi
}

%
%
% Honors program title page. Don't print the two approval pages.
% Should use option [honors] at top of thesis source file.
%
\long\def\HONORSTITLE{
\vbox to \textheight{
  \vspace{-0.75in}% Top margin is 1in from paper edge
  \begin{minipage}[c]{\minilength}%
    \begin{center}%
      \vbox to 9\baselineskip {%
          \HFmainhead\bf
          \vss
          \begin{center}
          \@title
          \end{center}
          \vss
          }%
      \mbox{} \\
    by\\
      \mbox{} \\
      \@author \\
      \vbox to 7\baselineskip
      {
        \vss
        A Senior Honors Project submitted to the faculty of\\
        The University of Utah\\
        in partial fulfillment of the requirements for the
        \vss
      }
      Honors Degree of \expandafter{\@degree} \\
      \mbox{} \\
      in\\
      \mbox{} \\
      \expandafter{\@honorsdepartment}
    \end{center}%
  \end{minipage}%
\vfill
{
\vskip 2\baselineskip
\noindent{\normalsize\rm APPROVED:}
\vskip 2\baselineskip
\begin{center}
\begin{tabular}{@{}p{7.2cm}@{\hspace{1cm}}p{7.2cm}}
\underline{\hspace{7cm}}\par
\expandafter{\@honorssupervisor}\par
 Supervisor
&
\underline{\hspace{7cm}}\par
\@departmentchair\par
Chairperson,  \expandafter{\@honorsdepartment}\\[2em]
\underline{\hspace{7cm}}\par
\expandafter{\@honorsadvisor}\par
 Departmental Honors Advisor
&
\underline{\hspace{7cm}}\par
\expandafter{\@honorsdirector}\par
 Director, Honors Program \\
\end{tabular}
\end{center}
}
  \unitlength = 1in
  \noindent
  \begin{picture}(5.85,0)(0,0)
    \put(0,-0.25){\mbox{\centerline{\@submitdate}}}%
  \end{picture}%
  \unitlength=1pt
}
\eject
\ifnoisy \typeout{Ignore Last Overfull hbox} \fi
\ifnoisy \typeout{Honors Title Page.} \fi
}


\long\def\reporttitlepage{
\thispagestyle{empty}
\reporttitle
\newpage
}

% Final defaults (not too important, since usually overridden).

\pagestyle{plain}

\pagenumbering{arabic}

\onecolumn

\if@twoside \else \raggedbottom \fi

\normalsize

% New label macro, GBG and Greg Conner, 17 April 1992:
%usage:
%\ulabel{MyLabelForRef}{RefMessage}
%Example:
%         \begin{theorem}% Theorem 2.4
%         \ulabel{hottheorem}{Theorem}
%         \end{theorem}
%         Later on, "\ref{hottheorem}" produces string "Theorem 2.4"
%         because "\thetheorem" equals "2.4".
%
\def\ulabel#1#2{\@bsphack\if@filesw {\let\thepage\relax
   \def\protect{\noexpand\noexpand\noexpand}%
   \xdef\@gtempa{\write\@auxout{\string
   \newlabel{#1}{{#2 \@currentlabel}{\thepage}}}}}\@gtempa
   \if@nobreak \ifvmode\nobreak\fi\fi\fi\@esphack}

%%% ============================================================= %%%
%%%   Local macro additions to UofU thesis style made by persons
%%%   who wrote theses. Use or change these according to taste.
%%% ============================================================= %%%

% This is a special two-page environment for when you paste a photograph
% on the second page. (thanx to JWP)
% Example: \twopagefigure{Gauss at Age Five}{gaussphoto}
%          Puts "Gauss at Age Five" in center of first page, makes
%          a blank second page for photo paste-up. Label is "gaussphoto".
%
\newcommand{\twopagefigure}[2]{
\begin{figure}[p]
\vspace{4in}
\caption{#1}
\label{#2}
\vspace{4in}
\end{figure}
\begin{figure}[p]
\vspace{8in}
\end{figure}
}

% Box for ends of proofs.

%\def\qed{\nobreak \hfil \rule{1ex}{1em}}
\def\qed{\protect\par\noindent{}$\square$\protect\par}

% New proof macros, Greg Conner, 17 April 1992
\def\proofline#1{\protect\par\noindent{\bfseries #1}:~}
\def\proof#1{\proofline{Proof #1}}
\def\pf{\protect\proofline{Proof}}

% Proof environment by Nelson Beebe Nov 1992
% A Proof environment is in normal Roman text, with a trailing
% black square at the end of the text, supplied automatically
% by the \end{Proof}.
\newcommand{\boxx}{\unskip \nopagebreak \hfill \rule{0.75em}{0.75em}}
\newenvironment{Proof}%
    {\noindent{\bfseries Proof.}\begin{rmfamily}}%
    {\boxx\end{rmfamily}}

\@addtoreset{equation}{chapter}

\def\theequation{\thechapter.\arabic{equation}}

% Turn @ back into a normal character and define \@ .

% Definitions for math environments.
% Turn on feature with control \theoremsetup in preamble

\def\theoremsetup{
\newtheorem{theorem}{Theorem}[chapter]
\newtheorem{proposition}{Proposition}[chapter]
\newtheorem{corollary}[proposition]{Corollary}
\ifnoisy
\message{Using default theorem, proposition, corollary.}
\fi
}

% \begin{topics}{template} ... \end{topics}
% used like enumerate, itemize, description environments
% The item[] is never optional!
%
\newcommand{\topicslabel}[1]{\mbox{#1}\hfil}
\newenvironment{topics}[1]{
    \begin{list}{}{%
      \let\makelabel\topicslabel
      \settowidth{\labelwidth}{\topicslabel{#1}}
      \setlength{\leftmargin}{\labelwidth}
      \advance\leftmargin by 1.5em}
}{\end{list}}

%% Used for 1997 Spiro thesis
%% Affects the start command of the environment only.
%% Usage: \begin{equation}\EMX .... \end{equation}
%% It is not desirable to insert "\EMX" into every equation
%% manually. For uniformity, do it with a macro definition.
%%
\def\EMX{%
\ifx\@optionONE\@ptsize
  \abovedisplayskip 24pt plus 3pt minus 7pt
\else
  \abovedisplayskip 22pt plus 2pt minus 5pt
\fi
\belowdisplayskip \abovedisplayskip
\abovedisplayshortskip  \abovedisplayskip
\belowdisplayshortskip  \belowdisplayskip
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Allow for negative page numbering.
%%% Removed by GBG March 5 1998 to restore LaTeX basic performance.

%% \def\@outputpage{%
%%   \let \protect \@unexpandable@protect
%%   \shipout \vbox{%
%%     \let \protect \relax
%%     \aftergroup\let \aftergroup\protect \aftergroup\relax
%%     \@shipoutsetup
%%     \vskip \topmargin
%%     \moveright\@themargin \vbox {%
%%       \setbox\@tempboxa \vbox to\headheight{%
%%         \vfil
%%         \color@begingroup
%%           \normalcolor
%%           \hbox to\textwidth {%
%%             \let \label \@gobble
%%             \let \index \@gobble
%%             \let \glossary \@gobble %% 21 Jun 91
%%             \@thehead
%%             }%
%%         \color@endgroup
%%         }%                        %% 22 Feb 87
%%       \dp\@tempboxa \z@
%%       \box\@tempboxa
%%       \vskip \headsep
%%       \box\@outputbox
%%       \baselineskip \footskip
%%       \color@begingroup
%%         \normalcolor
%%         \hbox to\textwidth{%
%%           \let \label \@gobble
%%           \let \index \@gobble      %% 22 Feb 87
%%           \let \glossary \@gobble   %% 21 Jun 91
%%           \@thefoot
%%           }%
%%       \color@endgroup
%%       }%
%%     }%
%%   \global \@colht \textheight
%% %
%%   \ifnum \c@page < \z@
%%     \addtocounter{page}\m@ne
%%   \else
%%     \addtocounter{page}\@ne
%%   \fi
%%   \begingroup\let\@elt\@stpelt \csname cl@page\endcsname\endgroup
%% %
%%   \let\firstmark\botmark
%% }

\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
\DeclareRobustCommand*\cal{\@fontswitch\relax\mathcal}
\DeclareRobustCommand*\mit{\@fontswitch\relax\mathnormal}

\makeatletter

%% end of uuthesis2eBETA.cls

%% uuthesis.history
%% ======================================================================
% Fixes by Grant Gustafson and Nelson Beebe, Oct, Nov 1993
% 1. TOC problem for \chapter and \part fixed
% 2. Added control sequences for mainheading width
% 3. Fixed chapter heading spacing.
% 4. Added control sequences for vertical spacing adjustment by line.
%    Example: \plusline adds \normalbaselineskip of vertical space
% 5. Fixed figure troubles and added uuguide.tex documentation.
% 6. Removed trash from 6 years of fixes.
% 7. Introduced font height parameters \largeHT and \normalHT.
% 8. New font controls for sectional commands.
% 9. Converted baselineskips to fixed spacing.
% 10. Found TOC error which made CONTENTS appear too high on page.
% 11. New controls for single spaced headings (where applicable).
% 12. Fixed TOC, LOF, LOT to maximum entry width 5.0in
% 13. Part starts the numbering not the chapter. Fixed TOC format.
% 14. CHAPTERS appears only once in the TOC.
% 15. Appendix fix for correct TOC content (no sections, etc).
% 16. Fixed "The" 3 times in final reading approval form.
% 17. Added \begin{topics}..\end{topics}, a list env like description.
% 18. Added 11pt support. All options are in the one style file, making
%     uut11.sty and uut12.sty obsolete. Other sizes are not possible.
% 19. Repaired \part, which had argument error (21-Nov-93).
% 20. Changed "\parskip 0pt plus 3pt minus 3pt" to 1pt tolerance
% 21. Fixed the APPROVAL pages to look like the Handbook samples.
% 22. Added \@warning hack to print line numbers for warnings.
% 23. Added \rawbibliography boolean to control Chemistry-style
%     references at the end of each chapter. In this case, references
%     are coded as a separate section. Jan 1994
% 24. Fixed involuntary change to 11pt font caused by frontmatter pages.
%     May 1994. Reported by Physics.
% 25. Removed \@warning hack, incompatible with latex2e.
% 26. Added [csreport] and [honors] options for making \reporttitlepage
%     into something special. Fixed automatic formatting to use the
%     correct option or else issue an error message. Sept 1994.
% 27. Added \nopagebreak[4] into section headings because of bug report
%     sent by Eric Eide.
% 28. Added \newpage into \titlepage, \copyrightpage and
%     \reporttitlepage because of a bug report sent by Eric Eide in
%     April 1995.
% 29. Added macros \EX and \EXM for math mode, to control uniformity of
%     spacing before and after a math display. Control \EX can be used
%     on inline equations to keep the display from encroaching on nearby
%     lines. Fix caused by demands on Jose Burillo math thesis, April
%     1996. Fix by GBG.
% 30. Added macro \EMX for math mode. Disabled \EX and \EMX. For Spiro
%     thesis May 1997. Creates "triple spaced equations".
%
% Fixes by Nelson H. F. Beebe <beebe@math.utah.edu> 22 June 1993
% 1. Remove duplicate obsolete code.
% 2. Insert critical \global\@nobreakfalse in \@chapter, thereby repairing
%    a long-standing problem of improper spacing before first section
%    in each chapter.
% 3. Correct before and after skip dimensions in all headings to bring
%    into line with UofU thesis requirements.
% Fixes by Nelson H. F. Beebe 01 June 1993
% 1. If \listoffiguresfalse or \listoftablesfalse has been specified,
%    do not typeset any exist .lof or .lot files
% Fixes 17-Apr-1992 by GBG and Greg Conner
% 1. Version printout changed to 1.4 27-April-1992
% 2. PART not printed for \part{} macro into TOC.
% 3. \title should be uppercase according to thesis editor.
% 4. Chapter headers now have control for spacing set by control
%    sequence \mainheadingwidth
% 5. PART uses ..... in TOC entry to match rest of TOC.
% 6. Control sequences \proofline, \proof, \pf, \ulabel added.
% 7. Control sequence \qed changed to reflect better ideas.
% 8. Fixed \vskip -22pt to \protect\vspace{-22pt} in TOC.
%    Caused underfull hboxes.
% 9. Fixed REFERENCES spacing, ABSTRACT spacing. Added \vbox
%    documentation to uuguide.tex to explain how to get chapter headings
%    spaced correctly for the thesis editor.
%
% Fix by Nelson H. F. Beebe 30-Apr-1990
% 1. Change + xxpt to plus xxpt in two places
% 2. Add double spacing for multiline chapter titles, as required by
%    UofU thesis office.
%   Fixes by David Eyre and Grant Gustafson, Feb-1990.
%
% 1. Strange \ps@arabic error found.
%    Fix: \markright{}, not \markright{\pagestyle{arabic}}
% 2. Uppercase error in two places.
% 3. 4-inch Thesis Format rule implemented with \parbox{4in}
% 4. Checked JDW 5/27/88 changes to see that all were implemented in
%    the July 1988 fixes by Stan Shebs.
% 5. Fixed three occurrences of subsection and subsubsection spacing errors.
%
% This file sets many of the formatting variables so that the output
% is acceptable for graduate school theses and dissertations.
% Original sources: REPORT.STY and REP12.STY.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixes by Kenneth W. Parker, May 1995.
% These changes where part of convertion from LaTeX2 to LaTeX2e
%
% 1. Declared options define \@ptsize with the \long atribute so \ifx
%    comparisions must be replaceed with \if comparisions.
% 2. The new size11.clo and size12.clo files declare \normalsize
%    \samll, and the like, so we shouldn't do so too.
% 3. Change \bf to \bfseries, \sc to \scshape, \rm to \rmfamily, and
%    \sl to \slshape within the body of the class.
% 4. Changed the numbering on the front-matter pages to be negitive.
%    This makes it posable to use dvips to select a range of page
%    numbers without conflicts.
% 5. Added the style-guild as a class option.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Things to work on:
%
% 1. old_fashion and modern options to turn on and off varable glue
%    and raggid bottems.
