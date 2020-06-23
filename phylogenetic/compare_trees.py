'''
Compare two trees by visualization and by 
Robinson-Foulds distance
'''

from Bio import Phylo
import matplotlib.pyplot as plt
import dendropy
from dendropy import Tree
from collections import defaultdict
import pandas as pd

num_seq = [10, 100, 200, 300, 400, 500, 600, 700, 900, 1000]
tree_dir = "../processed_data/phylogenetic/runtime_vs_num_seq/"
def draw_tree(tree_fn, out_fn):
	plt.close()
	tree = Phylo.read(tree_fn, "newick")
	Phylo.draw(tree)
	plt.tight_layout()
	plt.savefig(out_fn)
	plt.close()

def run_draw_tree(num_seq):
	for n in num_seq:
		print(f"Number of seq: {n}")
		iqtree_fn = f"{tree_dir}/{n}.np1.treefile"
		treebest_nj_fn = f"{tree_dir}/{n}.np1.nj.nhx"

		draw_tree(iqtree_fn, f"{tree_dir}/{n}.np1.png")
		draw_tree(treebest_nj_fn, f"{tree_dir}/{n}.np1.nj.png")


def get_robinson_foulds_distance(num_seq):
	container = defaultdict(list)
	for n in num_seq:
		iqtree_fn = f"{tree_dir}/{n}.np1.treefile"
		treebest_nj_fn = f"{tree_dir}/{n}.np1.nj.nhx"

		print(f"Number of seq: {n}")
		iqtree = Tree.get(path=iqtree_fn, schema="newick")
		treebest_nj = Tree.get(path=treebest_nj_fn, schema="newick", taxon_namespace=iqtree.taxon_namespace)
		dist = dendropy.calculate.treecompare.symmetric_difference(iqtree, treebest_nj)
		container["num"].append(n)
		container["robinson_foulds"].append(dist)

	robinson_foulds = pd.DataFrame(container)
	return robinson_foulds


def plot_robinson_foulds_distance(robinson_foulds):
	plt.close()
	fig, ax = plt.subplots(2,1, figsize=(4,6))
	ax[0].plot("num", "robinson_foulds", "go", data=robinson_foulds)
	ax[0].set_xlabel("number of sequence")
	ax[0].set_ylabel("Robinson-Foulds distance")

	ax[1].plot("num", "robinson_foulds", "go", data=robinson_foulds)
	ax[1].set_xscale("log")
	ax[1].set_yscale("log")
	ax[1].set_xlabel("log(number of sequence)")
	ax[1].set_ylabel("Robinson-Foulds distance")

	fig.tight_layout()
	fig.savefig(f"{tree_dir}/robinson_foulds_dist.png")


run_draw_tree(num_seq)
robinson_foulds = get_robinson_foulds_distance(num_seq)
plot_robinson_foulds_distance(robinson_foulds)
