#!/usr/bin/env python3

from sympy import *
from mpmath import *
from matplotlib.pyplot import *
#init_printing()     # make things prettier when we print stuff for debugging.


# ************************************************************************** #
# Error Analysis for Conductivity sigma                                      #
# ************************************************************************** #

# All values are in standard SI units unless otherwise noted.

# ---------------------------------------------------------#
# Constants                                                #
# ---------------------------------------------------------#
rho_kuchling_alu   = 0.027e-6  # resistivity Kuchling 17th edition, p.649, tab. 45
sigma_kuchling_alu = 1/rho_kuchling_alu
rho_kuchling_cu    = 0.0172e-6 # resistivity Kuchling 17th edition, p.649, tab. 45
sigma_kuchling_cu  = 1/rho_kuchling_cu

# http://www.aksteel.com/pdf/markets_products/stainless/austenitic/304_304L_Data_Sheet.pdf
# Converting from microOhm / inch to standard SI
rho_aksteel = 28.4 * 25.4 * 1e-3 * 1e-6
sigma_aksteel = 1/rho_aksteel
# http://hypertextbook.com/facts/2006/UmranUgur.shtml
sigma_glenEbert304 = 1.450e6
sigma_glenEbert347 = 1.392e6
sigma_glenEbert316 = 1.334e6
sigma_glenEbert = (sigma_glenEbert304+sigma_glenEbert347+sigma_glenEbert316) / 3
# http://www.dew-stahl.com/fileadmin/files/dew-stahl.com/documents/Publikationen/Werkstoffdatenblaetter/RSH/1.4301_de.pdf
rho_stahlwerke = 0.73e-6
sigma_stahlwerke = 1 / rho_stahlwerke

sigma_ref_st = ( sigma_aksteel + sigma_glenEbert + sigma_stahlwerke) / 3


# ---------------------------------------------------------#
# Load exported Values from Python Scripts                 #
# ---------------------------------------------------------#

files = [
    'numpy-txt/hollow--cu--freq--approx2.txt',
    'numpy-txt/hollow--cu--freq--approx.txt',
    'numpy-txt/hollow--cu--freq--exact.txt',
    'numpy-txt/hollow--st--freq--approx2.txt',
    'numpy-txt/hollow--st--freq--approx.txt',
    'numpy-txt/hollow--st--freq--exact.txt',
    'numpy-txt/massive--alu--freq.txt',
    'numpy-txt/massive--alu--high-freq-approx--high.txt',
    'numpy-txt/massive--alu--high-freq--exact.txt',
    'numpy-txt/massive--alu--low-freq--exact.txt'
    ]

valueDict = {}
for file in files:
    valueDict[file] = np.loadtxt(file)


# ---------------------------------------------------------#
# Functions                                                #
# ---------------------------------------------------------#
def getAbsError(array,avg):
    # array: numpy array
    # avg: number
    # Calculates absolute error for an array of values

    enum = np.array(array[0])
    for value in array:
        enum = np.append(enum,(value - avg)**2)

    denom = array.size * (array.size - 1)

    return sqrt(sum(enum)/denom)

def getStandardDeviation(array,avg):
    # array: numpy array
    # avg: number
    # Calculates standard deviation for an array of values

    enum = np.array(array[0])
    for value in array:
        enum = np.append(enum,(value - avg)**2)

    denom = array.size - 1

    return sqrt(sum(enum)/denom)


# ---------------------------------------------------------#
# Do Error Analysis                                        #
# ---------------------------------------------------------#

# copper, hollow:
sigmaCu = np.array(valueDict['numpy-txt/hollow--cu--freq--exact.txt'])
sigmaCu = np.append(sigmaCu,valueDict['numpy-txt/hollow--cu--freq--approx2.txt'])
sigmaCu = np.append(sigmaCu,valueDict['numpy-txt/hollow--cu--freq--approx.txt'])

sigmaAvgCu = sigmaCu.mean()

sigmaAbsErrCu = getAbsError(sigmaCu,sigmaAvgCu)
sigmaRelErrCu = sigmaAbsErrCu / sigmaAvgCu

sigmaStandDevCu = getStandardDeviation(sigmaCu,sigmaAvgCu)


# stainless steel, hollow:
sigmaSt = np.array(valueDict['numpy-txt/hollow--st--freq--exact.txt'])
sigmaSt = np.append(sigmaSt,valueDict['numpy-txt/hollow--st--freq--approx2.txt'])
sigmaSt = np.append(sigmaSt,valueDict['numpy-txt/hollow--st--freq--approx.txt'])

sigmaAvgSt = sigmaSt.mean()

sigmaAbsErrSt = getAbsError(sigmaSt,sigmaAvgSt)
sigmaRelErrSt = sigmaAbsErrSt / sigmaAvgSt

sigmaStandDevSt = getStandardDeviation(sigmaSt,sigmaAvgSt)


# aluminium, massive:
sigmaAlu = np.array(valueDict['numpy-txt/massive--alu--freq.txt'])
sigmaAlu = np.append(sigmaAlu,valueDict['numpy-txt/massive--alu--high-freq-approx--high.txt'])
sigmaAlu = np.append(sigmaAlu,valueDict['numpy-txt/massive--alu--high-freq--exact.txt'])
sigmaAlu = np.append(sigmaAlu,valueDict['numpy-txt/massive--alu--low-freq--exact.txt'])

sigmaAvgAlu = sigmaAlu.mean()

sigmaAbsErrAlu = getAbsError(sigmaAlu,sigmaAvgAlu)
sigmaRelErrAlu = sigmaAbsErrAlu / sigmaAvgAlu

sigmaStandDevAlu = getStandardDeviation(sigmaAlu,sigmaAvgAlu)


# LaTeX table
lines_table = [
        '        ' + r'$\overline{\sigma}_{Alu}    '    + r'$ & $\SI{'  +  str(sigmaAvgAlu)      + r'}{\ampere\per\volt\per\meter}$ \\' + "\n",
        '        ' + r'$s_{\overline{\sigma}_{Alu}}'    + r'$ & $\SI{'  +  str(sigmaAbsErrAlu)   + r'}{\ampere\per\volt\per\meter}$ \\' + "\n",
        '        ' + r'$r_{\overline{\sigma}_{Alu}}'    + r'$ & $\num{' +  str(sigmaRelErrAlu)   + r'}$ \\'                             + "\n",
        '        ' + r'$sdev_{\overline{\sigma}_{Alu}}' + r'$ & $\SI{'  +  str(sigmaStandDevAlu )+ r'}{\ampere\per\volt\per\meter}$ \\' + "\n",
        '        ' + r'$\overline{\sigma}_{Cu}     '    + r'$ & $\SI{'  +  str(sigmaAvgCu)       + r'}{\ampere\per\volt\per\meter}$ \\' + "\n",
        '        ' + r'$s_{\overline{\sigma}_{Cu}} '    + r'$ & $\SI{'  +  str(sigmaAbsErrCu)    + r'}{\ampere\per\volt\per\meter}$ \\' + "\n",
        '        ' + r'$r_{\overline{\sigma}_{Cu}} '    + r'$ & $\num{' +  str(sigmaRelErrCu)    + r'}$ \\'                             + "\n",
        '        ' + r'$sdev_{\overline{\sigma}_{Cu}} ' + r'$ & $\SI{'  +  str(sigmaStandDevCu ) + r'}{\ampere\per\volt\per\meter}$ \\' + "\n",
        '        ' + r'$\overline{\sigma}_{St}     '    + r'$ & $\SI{'  +  str(sigmaAvgSt)       + r'}{\ampere\per\volt\per\meter}$ \\' + "\n",
        '        ' + r'$s_{\overline{\sigma}_{St}} '    + r'$ & $\SI{'  +  str(sigmaAbsErrSt)    + r'}{\ampere\per\volt\per\meter}$ \\' + "\n",
        '        ' + r'$r_{\overline{\sigma}_{St}} '    + r'$ & $\num{' +  str(sigmaRelErrSt)    + r'}$ \\'                             + "\n",
        '        ' + r'$sdev_{\overline{\sigma}_{St}} ' + r'$ & $\SI{'  +  str(sigmaStandDevSt ) + r'}{\ampere\per\volt\per\meter}$ \\' + "\n",
        ]


# ---------------------------------------------------------#
# Save Results to LaTeX Table                              #
# ---------------------------------------------------------#

resultsOpening = r"""

"""

table_opening = r"""
{%
    \begin{center}
    \captionof{table}{%
        Mittelwerte, Fehler und Standardabweichung f\"ur die Leitf\"ahigkeiten
    }
    \label{tab:errorAnalysis}
    \sisetup{%
        %math-rm=\mathtt,
        scientific-notation=engineering,
        table-format = +3.4e+2,
        round-precision = 4,
        round-mode = figures,
    }
    \begin{tabular}{lr}
    \toprule
"""
table_closing = r"""
    \bottomrule
    \end{tabular}
    \end{center}
}
"""


dumpfileTable = open('listings/error-analysis.tex', 'w')
dumpfileTable.writelines(table_opening)
for line in lines_table:
    dumpfileTable.writelines(line)
dumpfileTable.writelines(table_closing)
dumpfileTable.close()


# -------------------------------------------------------- #
# Create Math-style Results                                #
#                                                          #
# NOTE: These values require editing  by hand in the *.tex #
# files  because the  siunitx  package in  LaTeX does  not #
# round values  with uncertainties  and the values  we get #
# from Python are far too precise for our needs.           #
# -------------------------------------------------------- #

equationCu_opening = r"""
\begin{equation}
    \sigma_{Cu}
        = \overline{\sigma}_{Cu} \pm s_{\overline{\sigma}_{Cu}}
        = \SI[separate-uncertainty = true]{"""
equationCu_closing = r" }{\ampere\per\volt\per\meter}" + "\n" + r"\end{equation}"

equationCu = equationCu_opening + str(sigmaAvgCu) + r" \pm " + str(sigmaAbsErrCu) + equationCu_closing

file = open('equations/resultCu.tex', 'w')
file.writelines(equationCu)
file.close()

equationSt_opening = r"""
\begin{equation}
    \sigma_{St}
        = \overline{\sigma}_{St} \pm s_{\overline{\sigma}_{St}}
        = \SI[separate-uncertainty = true]{"""
equationSt_closing = r" }{\ampere\per\volt\per\meter}" + "\n" + r"\end{equation}"

equationSt = equationSt_opening + str(sigmaAvgSt) + r" \pm " + str(sigmaAbsErrSt) + equationSt_closing

file = open('equations/resultSt.tex', 'w')
file.writelines(equationSt)
file.close()


equationAlu_opening = r"""
\begin{equation}
    \sigma_{Alu}
        = \overline{\sigma}_{Alu} \pm s_{\overline{\sigma}_{Alu}}
        = \SI[separate-uncertainty = true]{"""
equationAlu_closing = r" }{\ampere\per\volt\per\meter}" + "\n" + r"\end{equation}"

equationAlu = equationAlu_opening + str(sigmaAvgAlu) + r" \pm " + str(sigmaAbsErrAlu) + equationAlu_closing

file = open('equations/resultAlu.tex', 'w')
file.writelines(equationAlu)
file.close()


# ---------------------------------------------------------#
# Create Error Bar Plots                                   #
# ---------------------------------------------------------#

pictureStartStringCu = r"""
\begin{tikzpicture}
    \begin{axis}[
        try min ticks=2,
        width=.67\textwidth,
        height=.24\textwidth,
        title = {Leitf\"ahigkeit Kupfer, Vergleich experimenteller Wert mit Literaturwert},
        xlabel = {Leitf\"ahigkeit ($\si{\ampere\per\volt\per\meter}$)},
        symbolic y coords = {Literatur,Experiment}
    ]
    \addplot+[
        only marks,error bars/.cd,
        x dir=both,x explicit,
        error bar style={line width=0.5pt},
        ]
    coordinates {%
"""
pictureCloseStringCu = r"""
    };
    \end{axis}
\end{tikzpicture}
\captionof{figure}{%
    Vergleich  der  experimentell   bestimmten  Leitf\"ahigkeit  f\"ur  Kupfer
    mit  dem   Literaturwert  aus   Kuchlings  \emph{Taschebuch   der  Physik}
    \cite{ref:kuchling:resistivityTable}%
    }
"""
pictureCoordStringCu = '        (' + str(sigma_kuchling_cu) + ',Literatur)' + "\n" + '        (' + str(sigmaAvgCu) + ',Experiment) +- (' + str(sigmaAbsErrCu) + ',0)'
pictureStringCu = pictureStartStringCu + pictureCoordStringCu + pictureCloseStringCu

pictureStartStringSt = r"""
\begin{tikzpicture}
    \begin{axis}[
        try min ticks=2,
        width=.67\textwidth,
        height=.24\textwidth,
        title = {Leitf\"ahigkeit rostfreier Stahl, Vergleich experimenteller Wert mit Literaturwerten},
        xlabel = {Leitf\"ahigkeit ($\si{\ampere\per\volt\per\meter}$)},
        symbolic y coords = {Literatur,Experiment}
    ]
    \addplot+[
        only marks,error bars/.cd,
        x dir=both,x explicit,
        error bar style={line width=0.5pt},
        ]
    coordinates {%
"""
pictureCloseStringSt = r"""
    };
    \end{axis}
\end{tikzpicture}
\captionof{figure}{%
    Vergleich  der experimentell  bestimmten Leitf\"ahigkeit  f\"ur rostfreien
    Stahl  mit  einem  aus   Literaturwerten  bestimmten  Wert  (siehe  Anhang
    \ref{app:steel})%
    }
"""
pictureCoordStringSt = '        (' + str(sigma_ref_st) + ',Literatur)' + "\n" + '        (' + str(sigmaAvgSt) + ',Experiment) +- (' + str(sigmaAbsErrSt) + ',0)'
pictureStringSt = pictureStartStringSt + pictureCoordStringSt + pictureCloseStringSt


pictureStartStringAlu = r"""
\begin{tikzpicture}
    \begin{axis}[
        try min ticks=2,
        width=.67\textwidth,
        height=.24\textwidth,
        title = {Leitf\"ahigkeit Aluminium, Vergleich experimenteller Wert mit Literaturwert},
        xlabel = {Leitf\"ahigkeit ($\si{\ampere\per\volt\per\meter}$)},
        symbolic y coords = {Literatur,Experiment}
    ]
    \addplot+[
        only marks,error bars/.cd,
        x dir=both,x explicit,
        error bar style={line width=0.5pt},
        ]
    coordinates {%
"""
pictureCloseStringAlu = r"""
    };
    \end{axis}
\end{tikzpicture}
\captionof{figure}{%
    Vergleich  der experimentell  bestimmten  Leitf\"ahigkeit f\"ur  Aluminium
    mit  dem   Literaturwert  aus   Kuchlings  \emph{Taschebuch   der  Physik}
    \cite{ref:kuchling:resistivityTable}%
    }
"""
pictureCoordStringAlu = '        (' + str(sigma_kuchling_alu) + ',Literatur)' + "\n" + '        (' + str(sigmaAvgAlu) + ',Experiment) +- (' + str(sigmaAbsErrAlu) + ',0)'
pictureStringAlu = pictureStartStringAlu + pictureCoordStringAlu + pictureCloseStringAlu

file = open('plots-pgf/finalResultsCu.tex', 'w')
file.writelines(pictureStringCu)
file.close()

file = open('plots-pgf/finalResultsSt.tex', 'w')
file.writelines(pictureStringSt)
file.close()

file = open('plots-pgf/finalResultsAlu.tex', 'w')
file.writelines(pictureStringAlu)
file.close()
