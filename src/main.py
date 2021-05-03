import glob
import os
import subprocess
import sys

import click
import dendropy
from dendropy.utility import bitprocessing
import numpy as np

import CN3


def build_cn3_matrix(input_cnp, output_prefix):
    num_leaves = None
    cnp_dict = {}
    reading_profiles = False
    input_cnp_basename = os.path.basename(input_cnp)
    with open(input_cnp, "r") as f:
        for line in f:
            if(reading_profiles):
                current_line_arr = line.split()
                cnp_dict[current_line_arr[0]] = np.array(list(map(int, current_line_arr[2:])))
            else:
                # print(line)
                if("#number of leaves" in line):
                    num_leaves = int(line.split()[0])
                elif("#PROFILES" in line):
                    reading_profiles = True
    # print(num_leaves)
    # print(cnp_dict)
    with open(output_prefix + input_cnp_basename + "-cn3.mat", "w") as f:
        f.write(str(num_leaves) + "\n")
        for taxon_label,cnp in cnp_dict.items():
            f.write(taxon_label)
            padding_length = 10 - len(taxon_label)
            if(padding_length > 0):
                for padding_index in range(padding_length):
                    f.write(" ")
            for compare_to_taxon_label,compare_to_cnp in cnp_dict.items():
                try:
                    # valid_indices = np.where(np.array(cnp) + np.array(compare_to_cnp) != 0)[0]
                    valid_indices = np.where(cnp + compare_to_cnp != 0)[0]
                    # sanitized_cnp = np.array(cnp)[valid_indices]
                    sanitized_cnp = cnp[valid_indices]
                    # sanitized_compare_to_cnp = np.array(compare_to_cnp)[valid_indices]
                    sanitized_compare_to_cnp = compare_to_cnp[valid_indices]
                    current_distance = CN3.CN3(sanitized_cnp, sanitized_compare_to_cnp)[0]
                except:
                    return -1
                f.write(str(current_distance) + " ")
            f.write("\n")

def run_nj(input_cnp, output_prefix):
    input_cnp_basename = os.path.basename(input_cnp)
    distance_matrix_file = output_prefix + input_cnp_basename + "-cn3.mat"
    fastme_output_file = output_prefix + input_cnp_basename + "-fastme.tree"
    with open(output_prefix + input_cnp_basename + "-fastme.out", "w") as f_out:
        with open(output_prefix + input_cnp_basename + "-fastme.err", "w") as f_err:
            subprocess.call(["/usr/bin/time", "-v", "/home/min/git_repos/CS598MEBProject/fastme-2.1.6.2/binaries/fastme-2.1.6.2-linux64", "-i", distance_matrix_file, "-o", fastme_output_file], stdout=f_out, stderr=f_err)


def decompose_into_two(input_cnp, output_prefix):
    input_cnp_basename = os.path.basename(input_cnp)
    current_tree = output_prefix + input_cnp_basename + "-fastme.tree"
    full_tree = dendropy.Tree.get(path=current_tree, schema="newick")
    # full_tree.is_rooted = True
    full_tree.resolve_polytomies(limit=2)
    # full_tree.collapse_basal_bifurcation()
    full_tree.update_bipartitions(suppress_unifurcations=False)
    best_edge = get_centroid_edge(full_tree)
    left_tree,right_tree = bipartition_by_edge(full_tree, best_edge)
    left_tree.write(path=output_prefix + "left.tree", schema="newick")
    right_tree.write(path=output_prefix + "right.tree", schema="newick")
    left_leaves = bitprocessing.num_set_bits(left_tree.seed_node.tree_leafset_bitmask)
    right_leaves = bitprocessing.num_set_bits(right_tree.seed_node.tree_leafset_bitmask)

    # sys.stderr.write("left tree has " + str(left_leaves) + "\n")
    # sys.stderr.write("right tree has " + str(right_leaves) + "\n")
    if(left_leaves != right_leaves):
        return -1

    left_tree_set = set([n.taxon.label for n in left_tree.leaf_nodes()])
    right_tree_set = set([n.taxon.label for n in right_tree.leaf_nodes()])


    with open(input_cnp, "r") as f:
        with open(output_prefix + input_cnp_basename + "-left.cnp", "w") as f_left:
            with open(output_prefix + input_cnp_basename + "-right.cnp", "w") as f_right:
                is_profiles = False
                for line in f:
                    if("#PROFILES" in line):
                        is_profiles = True
                        f_left.write(line)
                        f_right.write(line)
                    elif(not is_profiles):
                        if("#number of leaves" in line):
                            num_leaves_suffix = " #number of leaves\n"
                            f_left.write(str(left_leaves) + num_leaves_suffix)
                            f_right.write(str(right_leaves) + num_leaves_suffix)
                        else:
                            f_left.write(line)
                            f_right.write(line)
                    else:
                        current_taxon = line.split()[0]
                        if(current_taxon in left_tree_set):
                            f_left.write(line)
                        elif(current_taxon in right_tree_set):
                            f_right.write(line)
                        else:
                            raise Exception("not in left or right")

def get_centroid_edge(tree):
    num_leaves = bitprocessing.num_set_bits(tree.seed_node.tree_leafset_bitmask)
    best_balance = float('inf')
    # sys.stderr.write("searching for best edge in num leaves:")
    # sys.stderr.write(str(num_leaves) + str("\n"))
    for edge in tree.postorder_edge_iter():
        if edge.tail_node is None:
            continue
        # sys.stderr.write(str(bitprocessing.num_set_bits(edge.bipartition.leafset_bitmask)) + "\n")
        balance = abs(num_leaves/2 - bitprocessing.num_set_bits(edge.bipartition.leafset_bitmask))
        # sys.stderr.write("current_balance:")
        # sys.stderr.write(str(balance) + "\n")
        if balance < best_balance:
            best_balance = balance
            best_edge = edge
    # sys.stderr.write(str(best_edge.head_node) + "\n")
    # sys.stderr.write(str(best_edge.length) + "\n")
    # sys.stderr.write(str(best_edge.head_node.label) + "\n")
    return best_edge


def bipartition_by_edge(tree, edge):
    new_root = edge.head_node
    edge.tail_node.remove_child(new_root)
    new_tree = dendropy.Tree(seed_node=new_root, taxon_namespace = tree.taxon_namespace)
    tree.update_bipartitions()
    new_tree.update_bipartitions()
    return tree,new_tree

def run_cnp_ilp(input_cnp, output_prefix):
    input_cnp_basename = os.path.basename(input_cnp)
    with open(output_prefix + input_cnp_basename + "-left.ilp.out", "w") as f_out:
        with open(output_prefix + input_cnp_basename + "-left.ilp.err", "w") as f_err:
            subprocess.call(["/usr/bin/time", "-v", "/home/min/git_repos/CS598MEBProject/CNT-ILP/build/cnt", "-s", str(2), output_prefix + input_cnp_basename + "-left.cnp"], stdout=f_out, stderr=f_err)
    with open(output_prefix + input_cnp_basename + "-right.ilp.out", "w") as f_out:
        with open(output_prefix + input_cnp_basename + "-right.ilp.err", "w") as f_err:
            subprocess.call(["/usr/bin/time", "-v", "/home/min/git_repos/CS598MEBProject/CNT-ILP/build/cnt", "-s", str(2), output_prefix + input_cnp_basename + "-right.cnp"], stdout=f_out, stderr=f_err)


def merge_files(input_cnp, output_prefix):
    input_cnp_basename = os.path.basename(input_cnp)
    root = None
    left_root = None
    right_root = None

    current_tree = output_prefix + input_cnp_basename + "-fastme.tree"
    full_tree = dendropy.Tree.get(path=current_tree, schema="newick")
    full_tree.is_rooted = False
    full_tree.resolve_polytomies(limit=2)
    full_tree.collapse_basal_bifurcation()
    full_tree.update_bipartitions()
    left_tree = dendropy.Tree.get(path=output_prefix + "left.tree", schema="newick")
    right_tree = dendropy.Tree.get(path=output_prefix + "right.tree", schema="newick")
    num_leaves = bitprocessing.num_set_bits(full_tree.seed_node.tree_leafset_bitmask)
    left_tree_set = sorted([n.taxon.label for n in left_tree.leaf_nodes()])
    right_tree_set = sorted([n.taxon.label for n in right_tree.leaf_nodes()])
    left_tree_num_leaves = len(left_tree_set)
    right_tree_num_leaves = len(right_tree_set)
    left_taxon_mapping = {}
    left_profiles = []
    right_profiles = []
    left_edges = []
    right_edges = []

    with open(output_prefix + input_cnp_basename + "-left.ilp.out", "r") as f_left:
        is_profiles = False
        is_edges = False
        is_events = False
        for line in f_left:
            if("#PROFILES" in line):
                is_profiles = True
                is_edges = False
                is_events = False
            elif("#EDGES" in line):
                is_profiles = False
                is_edges = True
                is_events = False
            elif("#EVENTS" in line):
                break
            elif(is_profiles):
                # print("adding profiles")
                if(left_root is None):
                    left_root = line.split()[2:]
                left_profiles.append(line.split()[2:])
            elif(is_edges):
                left_edges.append((int(line.split()[0]), line))

    with open(output_prefix + input_cnp_basename + "-right.ilp.out", "r") as f_right:
        is_profiles = False
        is_edges = False
        is_events = False
        for line in f_right:
            if("#PROFILES" in line):
                is_profiles = True
                is_edges = False
                is_events = False
            elif("#EDGES" in line):
                is_profiles = False
                is_edges = True
                is_events = False
            elif("#EVENTS" in line):
                break
            elif(is_profiles):
                if(right_root is None):
                    right_root = line.split()[2:]
                right_profiles.append(line.split()[2:])
            elif(is_edges):
                right_edges.append((int(line.split()[0]), line))

    left_edges = sorted(left_edges, key=lambda x: x[0])
    right_edges = sorted(right_edges, key=lambda x: x[0])
    left_offset = (2 * left_tree_num_leaves) - 1
    merged_edges = []
    for left_edge in left_edges:
        current_parent_node = left_edge[0]
        current_child_node = int(left_edge[1].split()[2])
        # fixed_parent = current_parent_node + 1
        # fixed_child = current_child_node + 1
        fixed_parent = (current_parent_node * 2) + 1
        fixed_child = (current_child_node * 2) + 1
        merged_edges.append((fixed_parent, fixed_child))
    for right_edge in right_edges:
        current_parent_node = right_edge[0]
        current_child_node = int(right_edge[1].split()[2])
        # fixed_parent = left_offset + current_parent_node + 1
        # fixed_child = left_offset + current_child_node + 1
        fixed_parent = (current_parent_node + 1) * 2
        fixed_child = (current_child_node + 1) * 2
        if((fixed_parent, fixed_child) not in merged_edges):
            merged_edges.append((fixed_parent, fixed_child))
    merged_edges = sorted(merged_edges, key=lambda x: x[0])

    # print("left profiles: " + str(left_profiles))
    # print("right profiles: " + str(right_profiles))
    # print("left edges: " + str(left_edges))
    # print("right edges: " + str(right_edges))
    # print("merged edges: " + str(merged_edges))

    # print(left_root)
    # print(right_root)
    left_cnp_input = [int(num) for num in left_root]
    right_cnp_input = [int(num) for num in right_root]
    root = CN3.CN3(left_cnp_input, right_cnp_input)[1]
    root_profile = [str(num) for num in root]
    # print(root_profile)

    left_counter = 0
    right_counter = 0
    with open(output_prefix + input_cnp_basename + "-merged.cnp", "w") as f_out:
        with open(input_cnp, "r") as f_in:
            for line in f_in:
                f_out.write(line)
                if("#PROFILES" in line):
                    break
        for i in range(2 * num_leaves - 1):
            f_out.write(str(i) + " : ")
            if(i == 0):
                f_out.write(" ".join(root_profile) + "\n")
            # elif(i < len(left_profiles) + 1):
            elif(i % 2 == 0):
                f_out.write(" ".join(left_profiles[left_counter]) + "\n")
                left_counter += 1
            else:
                f_out.write(" ".join(right_profiles[right_counter]) + "\n")
                right_counter += 1
        f_out.write("#EDGES\n")
        f_out.write("0 -> 1\n")
        f_out.write("0 -> 2\n")
        # f_out.write("0 -> " + str(left_offset + 1) + "\n")
        for merged_edge in merged_edges:
            f_out.write(str(merged_edge[0]) + " -> " + str(merged_edge[1]) + "\n")

def compare_output(input_cnp, output_prefix):
    input_cnp_basename = os.path.basename(input_cnp)
    true_cnp = "." + "".join(input_cnp.split(".")[:2]) + "." + input_cnp.split(".")[2].replace("input", "true")
    rf_rate = None
    with open(output_prefix + input_cnp_basename + "-merged.cnp.out", "w") as f_out:
        with open(output_prefix + input_cnp_basename + "-merged.cnp.err", "w") as f_err:
            subprocess.call(["/usr/bin/time", "-v", "/home/min/git_repos/CS598MEBProject/CNT-ILP/build/compare", true_cnp, output_prefix + input_cnp_basename + "-merged.cnp"], stdout=f_out, stderr=f_err)

    with open(output_prefix + input_cnp_basename + "-merged.cnp.out", "r") as f:
        for line in f:
            if ("RF =" in line):
                rf_rate = float(line.split()[2])
    return rf_rate




@click.command()
@click.option("--input-cnp", required=True, type=click.Path(exists=True), help="The input copy number profile file")
@click.option("--output-prefix", required=True, type=click.Path(), help="Output prefix")
def dac_cnt(input_cnp, output_prefix):
    build_cn3_matrix(input_cnp, output_prefix)
    # run_nj(input_cnp, output_prefix)
    # decompose_into_two(input_cnp, output_prefix)
    # run_cnp_ilp(input_cnp, output_prefix)
    # merge_files(input_cnp, output_prefix)
    # rf_rate = compare_output(input_cnp, output_prefix)
    # print("final rf rate is " + str(rf_rate))


@click.command()
@click.option("--input-cnp-dir", required=True, type=click.Path(exists=True), help="The input copy number profile directory")
@click.option("--output-prefix", required=True, type=click.Path(), help="Output prefix")
def dac_cnt_all(input_cnp_dir, output_prefix):
    cumulative_rf_rate = 0.0
    with open(output_prefix + "input_files.in", "w") as f:
        skip_count = 0
        for current_index,input_cnp in enumerate(glob.glob(input_cnp_dir + "*k4_n15*.input")):
            sys.stderr.write("currently on " + str(input_cnp) + ": ")
            if(build_cn3_matrix(input_cnp, output_prefix) == -1):
                sys.stderr.write("skipping due to cn3 \n")
                skip_count += 1
                continue
            run_nj(input_cnp, output_prefix)
            if (decompose_into_two(input_cnp, output_prefix) == -1):
                sys.stderr.write("skipping due to decompose \n")
                skip_count += 1
                continue
            run_cnp_ilp(input_cnp, output_prefix)
            merge_files(input_cnp, output_prefix)
            current_rf_rate = compare_output(input_cnp, output_prefix)
            cumulative_rf_rate += current_rf_rate
            sys.stderr.write("current rf rate is " + str(current_rf_rate) + " and ")
            sys.stderr.write("rf rate so far is " + str(cumulative_rf_rate / (current_index + 1 - skip_count)) + "\n")
            f.write(input_cnp + "\n")
    print("final mean rf rate is " + str(cumulative_rf_rate / (current_index + 1 - skip_count)))


if __name__ == "__main__":
    dac_cnt()
    # dac_cnt_all()

