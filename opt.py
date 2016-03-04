#!/usr/bin/env python
"""
Main script for the algorithm.

Handles parameters, and export of results to pdf and text files.
"""
from optlib import Optimize, capture_command, augment
import argparse
from inspect import cleandoc
import numpy as np

parser = argparse.ArgumentParser(
    usage="./opt.py [optional number of grid points] [optional parameters]")
parser.add_argument(
    "grid", nargs="?", default=64, type=int,
    help="Integer N giving NxN rectangular grid. The triangle is defined " +
         "by the main diagonal of that rectangle. Default value 64.")
parser.add_argument(
    "-d", "--divide", action="store", type=int, default=0,
    help="How many times to refine the grid for eigenvalue calculation. " +
         "Default value 0 only subdivides each rectangle into two congruent " +
         "triangles. Higher values force further spliting of each triangle " +
         "into four congruent triangles.")
parser.add_argument(
    "-t", "--threshold", action="store", default=12.25, type=float,
    help="The eigenvalue threshold leading tocontradiction with " +
         "domain monotonicity.")
parser.add_argument(
    "-l", "--left", action="store", type=int, default=93,
    help="The length of the left side.")
parser.add_argument(
    "-b", "--bottom", action="store", type=int, default=128,
    help="The length of the bottom side.")
args = parser.parse_args()

# script parameters
n = args.grid
print "Grid points: ", n
divide = args.divide
print "Mesh refinements: ", divide
threshold = args.threshold
print "Eigenvalue threshold: ", threshold

left_side_num = args.left
left_side_denom = args.bottom
left_side_len = left_side_num * 1.0 / left_side_denom

opt = Optimize(n, divide, threshold, left_side_num, left_side_denom)


def init_export():
    """ Write preamble to tex file. """
    f = open(opt.filename, "w")
    f.write(cleandoc("""
    \\documentclass{{article}}
    \\usepackage{{tikz}}
    \\begin{{document}}
    \\begin{{tikzpicture}}[scale=10]
        \\draw[clip] (0,{0}) |- (1,0) --cycle;
        % domains
        \\draw[help lines,yscale={0}] (0,0) grid[step={1}] (1,1);
    \\end{{tikzpicture}}
    \\begin{{enumerate}}
    \\end{{enumerate}}
    \\end{{document}}
    """.format(left_side_len, 1.0 / n)))
    f.close()
    capture_command(opt.latex)


def export(domain, eig):
    """ Extend the tex file and run latex. """
    f = open(opt.filename, "r")
    flist = f.readlines()
    f.close()
    # insert domain at line 6
    flist.insert(6, domain)
    # insert figure at line -2
    flist.insert(-2, eig)
    f = open(opt.filename, "w")
    f.write("".join(flist))
    f.close()
    capture_command(opt.latex)


def export_lower(lower, upper, raw, excluded, old_excluded, lower_eig, rect):
    """ Export figure and explanation for lower domains. """
    lower = zip(*[iter(lower[2:-2])] * 2)
    lower.append((0, left_side_len))
    domain = "\\fill[red,fill opacity=0.2] " + \
        "--".join([str(x) for x in lower]) + " --cycle;\n"
    upper = zip(*[iter(upper[:-2])] * 2)
    test_domain = "\\draw[very thick, red]" + \
                  "--".join([str(x) for x in upper]) + \
                  "coordinate (a);\n\draw[very thick, blue] (a) -- " + \
                  "(0,{}) -- (0,0);\n".format(left_side_len)

    f = open(opt.genericname + '_upper_plots.tex', 'a')
    f.write("\n" + cleandoc(r"""
    \begin{tikzpicture}[scale=5]
    \clip
    """) + "--".join([str(x) for x in upper]) +
            " -- (0,{});\n".format(left_side_len) +
            test_domain.replace('very ', '') + cleandoc(r"""
    \draw[help lines,yscale={},opacity=0.5] (0,0) grid[step={}] (1,1);
    \end{{tikzpicture}}
    """.format(left_side_len, 1.0/n))
            )
    f.close()

    f = open(opt.textname, 'a')
    raw = raw[:-1] + [[(0, n)]]
    print raw
    f.write('Upper domain:\n' + str(raw) + '\n')
    f.close()

    excluded = zip(*[iter(excluded[:-2])] * 2)
    excluded.extend([(1, 0), (0, 0)])
    excluded_domain = "\\fill[red,fill opacity=0.3] " + \
        "--".join([str(x) for x in excluded]) + " -- cycle;\n"
    vertex = (float(rect[1]) / n, float(rect[0]) / n * left_side_len)
    old_excluded = zip(*[iter(old_excluded[:-2])] * 2)
    old_excluded.extend([(0, left_side_len)])
    old_excluded_domain = "\\fill[magenta,fill opacity=0.3] " + \
        "--".join([str(x) for x in old_excluded]) + " -- cycle;\n"
    eig_string = cleandoc("""
    \\newpage
    \\item New vertex $p$ in upper excluded region $U^{{{}}}$:
    {} (in grid coordinates).

    Vertices of the upper domain $D_U(p)$:

    {}

    Unprocessed eigenvalues for the upper domain:

    {}

    Residuals:

    {}

    Rescaled eigenvalues (by the bottom side length), but not postprocessed:

    {}

    Postprocesses eigenvalue (lower bound):

    {}

    Eigenvalue for lower test domain $D_L^{{test}}$: {} (does not need to
    be calculated until it reaches the threshold {}).
    """.format((counter+1)/2, rect[::-1], raw, opt.raweigs, opt.residual,
               np.array(opt.raweigs)*left_side_denom**2/opt.solver.scale,
               opt.lowerbound, lower_eig, threshold)) + """

    \\begin{{tikzpicture}}[scale=10]
    \\begin{{scope}}
    \\draw[clip] (0,{0}) |- (1,0) --cycle;
    \\draw[help lines,yscale={0}] (0,0) grid[step={1}] (1,1);
    """.format(left_side_len, 1.0/n) + excluded_domain + \
        old_excluded_domain + """
    \\fill[green!50!black,opacity=0.3] (0,{}) rectangle {} node
    [shape=circle,fill=black,fill opacity=1,minimum size=1.5mm,inner sep=0pt,
    outer sep=0pt] {{}} node [above left=-3pt,black,opacity=1]
    """.format(left_side_len, vertex) + " { $V$};\\end{scope}" + \
        test_domain + "\\end{tikzpicture}\n\n"

    export(domain, eig_string)
    return lower_eig


def export_upper(upper, lower, raw, excluded, old_excluded, upper_eig, rect):
    """ Export figure and explanation for upper excluded set. """
    upper = zip(*[iter(upper[:-2])] * 2)  # put points together
    upper.extend([(1, 0), (0, 0)])  # add missing points
    domain = "\\fill[red,fill opacity=0.2] " + \
        "--".join([str(x) for x in upper]) + " -- cycle;\n"
    lower = zip(*[iter(lower[2:-2])] * 2)
    test_domain = "\\draw[very thick, red]" + \
                  "--".join([str(x) for x in lower]) + \
                  "coordinate (a);\n\draw[very thick, blue] (a) -- (1,0) " + \
                  "-- (0,0) -- {};\n".format(str(lower[0]))

    f = open(opt.genericname + '_lower_plots.tex', 'a')
    f.write("\n" + cleandoc(r"""
    \begin{tikzpicture}[scale=5]
    \clip (1,0) -- (0,0) --
    """) + "--".join([str(x) for x in lower]) + ";\n" +
            test_domain.replace('very ', '') + cleandoc(r"""
    \draw[help lines,yscale={},opacity=0.5] (0,0) grid[step={}] (1,1);
    \end{{tikzpicture}}
    """.format(left_side_len, 1.0/n))
            )
    f.close()

    f = open(opt.textname, 'a')
    raw = [(y, x) for x, y in raw[1:-1]] + [[(n, 0), (0, 0)]]
    print raw
    f.write('Lower domain:\n' + str(raw) + '\n')
    f.close()

    vertex = (float(rect[0]) / n, float(rect[1]) / n * left_side_len)
    excluded = zip(*[iter(excluded[:-2])] * 2)
    excluded.extend([(0, left_side_len)])
    excluded_domain = "\\fill[red,fill opacity=0.3] " + \
        "--".join([str(x) for x in excluded]) + " -- cycle;\n"
    old_excluded = zip(*[iter(old_excluded[:-2])] * 2)
    old_excluded.extend([(1, 0), (0, 0)])
    old_excluded_domain = "\\fill[magenta,fill opacity=0.3] " + \
        "--".join([str(x) for x in old_excluded]) + " -- cycle;\n"
    eig_string = cleandoc("""
    \\newpage
    \\item New vertex $p$ in lower excluded region $L^{{{}}}$:
    {} (in grid coordinates).

    Vertices for the lower domain $D_L(p)$:

    {}

    Unprocessed eigenvalues for the lower domain:

    {}

    Residuals:

    {}

    Rescaled eigenvalues (by bottom side length), but not postprocessed:

    {}

    Postprocesses eigenvalue (lower bound):

    {}

    Eigenvalue for upper test domain $D_U^{{test}}$: {} (does not need to
    be calculated until it reaches the threshold {}).

    """.format((counter+1)/2, rect, raw, opt.raweigs, opt.residual,
               np.array(opt.raweigs)*left_side_denom**2/opt.solver.scale,
               opt.lowerbound, upper_eig, threshold)) + """

    \\begin{{tikzpicture}}[scale=10]
    \\begin{{scope}}
    \\draw[clip] (0,{0}) |- (1,0) --cycle;
    \\draw[help lines,yscale={0}] (0,0) grid[step={1}] (1,1);
    """.format(left_side_len, 1.0/n) + excluded_domain + \
        old_excluded_domain + """
    \\fill[green!50!black,opacity=0.3] (1,0) rectangle {} node
    [shape=circle,fill=black,fill opacity=1,minimum size=1.5mm,
    inner sep=0pt,outer sep=0pt] {{}} node [below right=-3pt,black,opacity=1]
    """.format(vertex) + " { $V$};\\end{scope}" + \
        test_domain + "\\end{tikzpicture}\n\n"
    export(domain, eig_string)
    return upper_eig


def final(which=None):
    """
    Finish tex file when algorithm succeeded.

    Add the plot of the domain giving contradiction.
    """
    opt.stiff = opt.teststiff
    opt.mass = opt.testmass
    opt.save_matrices(counter)
    f = open(opt.filename, "r")
    flist = f.readlines()
    f.close()
    ind = [i for i, elem in enumerate(flist) if 'enumerate' in elem][0]
    if which is "lower":
        lower = zip(*[iter(opt.lower_domain()[2:-2])] * 2)
        txt = open(opt.textname, 'a')
        raw = [(y, x) for x, y in opt.lower_domain(raw=True)[1:-1]] + \
            [[(n, 0), (0, 0)]]
        print raw
        txt.write('Lower domain:\n' + str(raw) + '\n')
        txt.close()

        flist.insert(ind, "\n\n" + cleandoc("""
    We found the excluded regions shown above.

    \\newpage
    The complement of the upper
    excluded region (picture below) has the smallest eigenvalue above the
    threshold {}. Yet it must contain the nodal domain which has eigenvalue
    below that threshold. This gives a contradiction with the domain
    monotonicity.

    Vertices:

    {}

    Unprocessed eigenvalues:

    {}

    Residuals:

    {}

    Rescaled eigenvalues (but not postprocessed):

    {}

    Postprocessed eigenvalue (lower bound):

    {}

    """.format(threshold, raw, opt.testraweigs, opt.testresidual,
               np.array(opt.testraweigs), new_lower_eig)) + """

    \\begin{{tikzpicture}}[scale=10]
    \\begin{{scope}}
    \\draw[clip] (0,{0}) |- (1,0) --cycle;
    \\draw[help lines,yscale={0}] (0,0) grid[step={1}] (1,1);
    \\end{{scope}}
    \\draw[very thick, red]""".format(left_side_len, 1.0 / n) +
            "--".join([str(x) for x in lower]) +
            "coordinate (a); \draw[very thick,blue] (a) -- (1,0) -- " +
            "(0,0) -- {};\n".format(str(lower[0])) + "\\end{tikzpicture}")
    else:
        upper = zip(*[iter(opt.upper_domain()[:-2])] * 2)

        txt = open(opt.textname, 'a')
        raw = opt.upper_domain(raw=True)[:-1] + [[(0, n)]]
        txt.write('Upper domain:\n' + str(raw) + '\n')
        txt.close()

        flist.insert(ind, "\n\n" + cleandoc("""
    We found the excluded regions shown above.

    \\newpage
    The complement of the lower
    excluded region (picture below) has the smallest eigenvalue above the
    threshold {}. Yet it must contain the nodal domain which has eigenvalue
    below that threshold. This gives a contradiction with the domain
    monotonicity.

    Vertices:

    {}

    Unprocessed eigenvalues:

    {}

    Residuals:

    {}

    Rescaled eigenvalues (but not postprocessed):

    {}

    Postprocessed eigenvalue (lower bound):

    {}

    """.format(threshold, raw, opt.testraweigs, opt.testresidual,
               np.array(opt.testraweigs)*left_side_denom**2,
               new_upper_eig)) + """
    \\begin{{tikzpicture}}[scale=10]
    \\begin{{scope}}
    \\draw[clip] (0,{0}) |- (1,0) --cycle;
    \\draw[help lines,yscale={0}] (0,0) grid[step={1}] (1,1);
    \\end{{scope}}
    \\draw[very thick, red]""".format(left_side_len, 1.0 / n) +
            "--".join([str(x) for x in upper]) +
            "coordinate (a); \draw[very thick,blue] (a) -- " +
            "(0,{}) -- (0,0);\n".format(left_side_len) + "\\end{tikzpicture}")
    f = open(opt.filename, "w")
    f.write("".join(flist))
    f.close()
    capture_command(opt.latex)
    exit(0)


# main loop
f = open(opt.genericname + '_lower_plots.tex', 'w')
f.close()
f = open(opt.genericname + '_upper_plots.tex', 'w')
f.close()
init_export()
old_upper_eig = new_upper_eig = 0
old_lower_eig = new_lower_eig = 0

counter = 1
while True:

    new_lower_eig, left_rect = opt.update_left_side(old_lower_eig)
    opt.save_matrices(counter)
    old_left_side = opt.left_side
    opt.left_side = augment(opt.left_side, n, *left_rect)
    test_side = augment(opt.bottom_side, n, left_rect[1], left_rect[0])
    export_lower(opt.lower_domain(), opt.upper_domain(test_side),
                 opt.upper_domain(test_side, raw=True), opt.upper_domain(),
                 opt.lower_domain(old_left_side), new_lower_eig, left_rect)
    print "Finished step: ", counter
    counter += 1
    if new_lower_eig > threshold:
        final("lower")

    new_upper_eig, bottom_rect = opt.update_bottom_side(old_upper_eig)
    opt.save_matrices(counter)
    old_bottom_side = opt.bottom_side
    opt.bottom_side = augment(opt.bottom_side, n, *bottom_rect)
    test_side = augment(opt.left_side, n, bottom_rect[1], bottom_rect[0])
    export_upper(opt.upper_domain(), opt.lower_domain(test_side),
                 opt.lower_domain(test_side, raw=True), opt.lower_domain(),
                 opt.upper_domain(old_bottom_side), new_upper_eig, bottom_rect)
    print "Finished step: ", counter
    counter += 1
    if new_upper_eig > threshold:
        final("upper")

    if abs(new_lower_eig - old_lower_eig) + \
            abs(new_upper_eig - old_upper_eig) < 1e-10:
        # algorithm failed
        exit(1)
    old_upper_eig = new_upper_eig
    old_lower_eig = new_lower_eig
