import math
from operator import mul


def get_substitution_matrix(t, alpha=(1/3)):
    '''
    :param t: time (length of edge)
    :param alpha: alpha for Jukes-Cantor matrix, defaults to 1/3
    :return: Jukes-Cantor matrix
    '''
    diagonal = (1 + 3*math.exp((-4.0/3.0)*alpha*t))/4
    other = (1 - math.exp((-4.0/3.0)*alpha*t))/4
    matrix = [[diagonal if i == j else other for i in range(4)] for j in range(4)]
    return matrix


def extract_alignment_data(alignment_file="cftr.txt"):
    '''
    :param alignment_file: path to alignment file
    :return:    names of animals contained in file,
                sequences translated from strings of "ACGTN-" to lists of numbers <0, 5>
    '''
    states = 'ACGTN-'
    names = []
    sequences = []
    with open(alignment_file) as seqs:
        for line in seqs:
            if line.startswith(">"):
                names.append(line[1:])
                sequences.append("")
            else:
                sequences[-1] += line.strip()

    num_sequences = []
    for seq in sequences:
        num_sequences.append([])  # 699013 bases each
        for base in seq:
            num_sequences[-1].append(states.index(base))
    return names, num_sequences


def get_columns(num_sequences):
    '''
    :param num_sequences: sequences translated from strings of "ACGTN-" to lists of numbers <0, 5>
    :return: columns of alignment data based on input sequences
    '''
    return list(zip(*num_sequences))  # cols to rows


def compute_args(alpha=(1/3), tree_file="strom.txt"):
    '''
    Computes arguments for Felsentein algorithm.
    :param alpha: alpha for Jukes-Cantor matrices
    :param tree_file: path to file containing tree
    :return:    list of ordered vertices with children before parents,
                map of each vertex to list of its children,
                map of edges to Jukes-Cantor matrices,
                (0.25, 0.25, 0.25, 0.25) equilibrium,
                list of leaves ordered from left to right
    '''
    ov_nums = {}
    vertex_to_children_nums = {}
    vertex_to_children = {}
    edge_to_matrix = {}
    q_a = [1/4, 1/4, 1/4, 1/4]  # values in root
    i = 0
    with open(tree_file) as tree:
        for line in tree:
            child, parent, _ = line.split()
            if child not in ov_nums:
                ov_nums[child] = i
                i += 1
            if parent not in ov_nums:
                ov_nums[parent] = i
                i += 1
        tree.seek(0)
        for line in tree:
            child, parent, t = line.split()
            t = float(t)
            if parent in vertex_to_children:
                vertex_to_children[parent].append(child)
                vertex_to_children_nums[ov_nums[parent]].append(ov_nums[child])
            else:
                vertex_to_children[parent] = [child]
                vertex_to_children_nums[ov_nums[parent]] = [ov_nums[child]]
            matrix = get_substitution_matrix(t, alpha)
            edge_to_matrix[(ov_nums[parent], ov_nums[child])] = matrix
    ov = [ov_nums["Root"]]
    i = 0
    while i < len(ov):
        values = vertex_to_children_nums.get(ov[i], None)
        ov = ov + values if values is not None else ov
        i += 1

    leaves = get_ordered_leaves(ov_nums["Root"], vertex_to_children_nums)

    return list(reversed(ov)), vertex_to_children_nums, edge_to_matrix, q_a, leaves


def get_ordered_leaves(root, tree):
    '''
    :param root: tree root number
    :param tree: map of nodes (parent) to list of nodes (choldren)
    :return: list of leaves ordered from left to right
    '''
    def subtree_preorder(node):
        values = tree.get(node, None)
        if values is None:
            leaves.append(node)
        else:
            for n in values:
                subtree_preorder(n)

    leaves = []
    subtree_preorder(root)
    return leaves


def felsenstein(ordered_vertices, vertex_to_children, edge_to_matrix, q_a, leaves, column):
    '''
    :param ordered_vertices: vertices ordered so that children come before parents
    :param vertex_to_children: map of each vertex to list of its children
    :param edge_to_matrix: map of edges to Jukes-Cantor matrices
    :param q_a: (0.25, 0.25, 0.25, 0.25) equilibrium
    :param leaves: leaf numbers ordered left to right
    :param column: one column of alignment (e.g. first, second, third... base from each animal)
    :return: likelihood
    '''
    column = list(column)  # one col of alignment ("ACGTN-" -> 012345)
    nvertices = len(ordered_vertices)
    nstates = len(q_a)  # number of bases: 4 -> ACGT
    root = ordered_vertices[-1]
    leaves_map = dict(zip(leaves, column))
    pattern = [leaves_map.get(num, -1) for num in range(nvertices)]
    A = [[1.0 for _ in range(nstates)] for _ in range(nvertices)]
    for v in ordered_vertices:
        state = pattern[v]
        if state != -1:  # leaf
            if state == 4 or state == 5:  # N or -
                for s in range(nstates):
                    A[v][s] = 1
            else:  # is leaf and has one of the bases - ACGT, one set to 1 other to 0
                for s in range(nstates):
                    if s != state:
                        A[v][s] = 0
        else:  # inner node
            for a in range(nstates):
                for y in vertex_to_children.get(v, []):
                    P = edge_to_matrix[v, y]
                    dot = sum(map(mul, P[a], A[y]))  # dot = \sum_{b} A[y, b] * P(b|a, t_y)
                    A[v][a] *= dot

    dot = sum(map(mul, q_a, A[root]))  # compute likelyhood in root
    return dot


def find_best_alpha(tree_file, columns):
    """
    :param tree_file: e.g. strom.txt
    :param columns: list of tuples, where each tuple is sequence of numbers in <0, 5>, where numbers correspond to
                    states ACTGN-, each tuple of lenght x is single column of alignment of x animals,
                    e.g. first column in cftr.txt: (0, 0, 0, 0, 0, 0, 5, 5, 0, 0, 0, 5, 5, 5)
    :return: best alpha
    """
    alpha_to_logll = {}
    for alpha in [x * 0.1 for x in range(1, 21)]:
        alpha = round(alpha, 1)
        args = compute_args(alpha, tree_file)
        alpha_to_logll[alpha] = 0
        for column in columns:
            logll = math.log(felsenstein(*args, column))
            alpha_to_logll[alpha] += logll
    return max(alpha_to_logll, key=lambda key: alpha_to_logll[key])


def find_best_alpha_for_window(columns, tree_file="strom.txt"):
    return find_best_alpha(tree_file, columns)


def divide_windows(sequences, exon_file):
    """
    Divides windows into "exon" windows and "non-exon" windows
    :param sequences: list sequences
    :param exon_file: str, path to file with exon positions
    :return:    list, where each element is tuple, where first element of tuple is starting position in sequence and
                second element is list of lists of sequences from each animal
    """
    sequence = sequences[0]
    exon_pos = []
    with open(exon_file) as exons:
        for line in exons:
            start, end = map(int, line.split())
            exon_pos.append((start, end))
    window_groups = {"exon": [], "non-exon": []}
    exon_indexes = []
    i = 0
    for base in sequence:
        if base != 5:
            i += 1
        exon_indexes.append(i)

    i = 0
    while i < len(sequence):
        window = [s[i:i + 100] for s in sequences]
        window_start = exon_indexes[i]
        window_end = exon_indexes[i+100] if len(exon_indexes) > i+100 else exon_indexes[len(exon_indexes)-1]
        exon = False
        for e_start, e_end in exon_pos:
            if e_start <= window_start <= e_end or e_start <= window_end <= e_end:
                window_groups["exon"].append(tuple((i, window)))
                exon = True
        if not exon:
            window_groups["non-exon"].append(tuple((i, window)))
        i += 100
    return window_groups


print("Column: {pat}, Likelihood: {f} for dashes".format(pat=tuple(5 for _ in range(14)), f=felsenstein(*compute_args(), tuple(5 for _ in range(14)))))
print("Column: {pat}, Likelihood: {f} for base A at Alpha 1.0".format(pat=tuple(0 for _ in range(14)), f=felsenstein(*compute_args(1.0), tuple(0 for _ in range(14)))))
print("Column: {pat}, Likelihood: {f} for base A at Alpha 0.2".format(pat=tuple(0 for _ in range(14)), f=felsenstein(*compute_args(0.2), tuple(0 for _ in range(14)))))


columns = get_columns(extract_alignment_data("cftr.txt")[1])
for i in range(10):
    print("Best alpha for window {w_num}: {alpha}".format(w_num=i+1, alpha=find_best_alpha_for_window(columns[i*100:i*100+100])))


groups = divide_windows(extract_alignment_data("cftr.txt")[1], "exony.txt")
alpha_freq_table = {(round(x*0.1, 2)): 0 for x in range(1, 21)}

for exon_window in groups["exon"]:
    alpha = find_best_alpha_for_window(get_columns(exon_window[1]))
    alpha_freq_table[alpha] += 1


# Diagrams
import matplotlib.pyplot as plt

plt.rcdefaults()
import numpy as np

objects = list(map(str, alpha_freq_table.keys()))
y_pos = np.arange(len(objects))
performance = list(alpha_freq_table.values())

fig, ax = plt.subplots()
plt.bar(y_pos, performance, align='center', alpha=0.5)
plt.xticks(y_pos, objects)
plt.ylabel('Count')
plt.xlabel('Alpha')
plt.title('Exon windows')
plt.show()













