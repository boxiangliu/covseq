import sys
sys.path.append(".")
from phylogenetic.construct_tree import msa, construct_tree
from Bio import SeqIO
import os
import time
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
import subprocess
import pickle

gisaid_fasta_fn = "../data/gisaid/fasta/sequences.fasta"
out_dir = "../processed_data/phylogenetic/runtime_vs_num_seq/"
beijing_fasta_fn = "../data/BJ_CDC/4seq.fasta"
os.makedirs(out_dir, exist_ok=True)
num_seq = [10, 100, 200, 300, 400, 500, 600, 700, 900, 1000]



def sample_fasta(gisaid_fasta_fn, out_dir, num_seq, lb=25000, ub=30000):
	'''
	lb: sequence length lower bound 
	ub: sequence length upper bound
	'''
	out_fn_list = []
	for n in num_seq:
		out_fn = f"{out_dir}/{n}.fasta"
		out_fn_list.append(out_fn)
		with open(out_fn, "w") as fout:
			records = SeqIO.parse(gisaid_fasta_fn, "fasta")
			i = 0 
			while i < n:
				record = next(records)
				if len(record.seq) > lb and len(record.seq) < ub:
					record.id = record.id.replace("/", "_")
					record.description = record.description.replace("/", "_")
					SeqIO.write(record, fout, "fasta")
					i += 1
				else:
					pass
	return out_fn_list


def get_msa_runtime(mode, num_seq, fasta_fn_list=None,\
	msa_fn_list=None, new_sequences=None):

	assert mode in ["denovo", "np1"], \
		"The argument mode must be denovo (from scratch) " \
		"or np1 (appending)."

	if mode == "denovo":
		print("mode = denovo.")
		assert fasta_fn_list is not None, \
			"The argument fasta_fn_list must be provided."
		fn_list = fasta_fn_list
	else:
		print("mode = add.")
		assert (msa_fn_list is not None) and \
			(new_sequences is not None), \
			"Arguments msa_fn_list and new_sequences must "\
			"be provided."
		fn_list = msa_fn_list

	container = defaultdict(list)
	out_fn_list = []

	for n, fn in zip(num_seq, fn_list):
		print(f"Number of sequence: {n}")
		container["num"].append(n)

		start = time.time()
		if mode == "denovo":
			msa_fn = fn.replace(".fasta", ".ali")
			msa(mode="denovo", fasta_fn=fn, align_fn=msa_fn)
		else:
			msa_fn = fn.replace(".ali", ".np1.ali")
			msa(mode="add", fasta_fn=new_sequences, align_fn=msa_fn,\
				existing_alignment=fn)

		out_fn_list.append(msa_fn)
		msa_time = time.time()
		container["msa_time"].append(msa_time - start)

	runtime = pd.DataFrame(container)
	return runtime, out_fn_list


def get_tree_runtime(msa_fn_list, num_seq, software):

	assert software in ["iqtree", "FastTree"], \
		"The software argument must be iqtree or FastTree."

	container = defaultdict(list)
	tree_fn_list = []
	for n, msa_fn in zip(num_seq, msa_fn_list):
		print(f"Number of sequence: {n}")
		container["num"].append(n)

		start = time.time()

		if software == "iqtree":
			tree_prefix = msa_fn.replace(".ali", "")
			construct_tree(software="iqtree", msa_fn=msa_fn, \
				out_prefix=tree_prefix)
			tree_fn_list.append(tree_prefix)
		else:
			tree_fn = msa_fn.replace(".ali", ".fastTree.nh")
			tree_log = msa_fn.replace(".ali", ".fastTree.log")
			construct_tree(software="FastTree", msa_fn=msa_fn, \
				out_fn=tree_fn, log_fn=tree_log)
			tree_fn_list.append(tree_fn)

		tree_time = time.time()
		container["tree_time"].append(tree_time - start)
	
	runtime = pd.DataFrame(container)
	return runtime, tree_fn_list


def plot_runtime_vs_num_seq(runtime, x, y, out_fn):
	'''
	x: a string scalar representing the column name for x axis.
	y: a string list representing column names for y axes.
	'''
	plt.close()
	width = 3 * len(y)
	height = 3 * 2
	fig, ax = plt.subplots(nrows=2, ncols=len(y), figsize=(width,height))
	
	for iy in y:
		ax[0,0].plot("num", "msa_time", "go", data=runtime)
		ax[0,0].set_title("MSA")
		ax[0,0].set_xlabel("Number of sequences")
		ax[0,0].set_ylabel("Runtime (in S)")

		ax[1,0].plot("num", "msa_time", "go", data=runtime)
		ax[1,0].set_yscale("log")
		ax[1,0].set_xscale("log")
		ax[1,0].set_title("MSA (log-log)")
		ax[1,0].set_xlabel("log(number of sequences) (in S)")
		ax[1,0].set_ylabel("log(runtime)")

	ax[0,1].plot("num", "tree_time", "go", data=runtime)
	ax[0,1].set_title("Phylogenetic Tree")
	ax[0,1].set_xlabel("Number of sequences")
	ax[0,1].set_ylabel("Runtime (in S)")

	ax[1,1].plot("num", "tree_time", "go", data=runtime)
	ax[1,1].set_yscale("log")
	ax[1,1].set_xscale("log")
	ax[1,1].set_title("Phylogenetic Tree (log-log)")
	ax[1,1].set_xlabel("log(number of sequences) (in S)")
	ax[1,1].set_ylabel("log(runtime)")

	plt.tight_layout()
	fig.savefig(out_fn)


def treebest(command, input_file, output_file, **kwargs):
	assert command in ["nj", "phyml"]
	cmd = ["treebest"]

	# Add command
	cmd.append(command)

	# Add options
	for k, v in kwargs.items():
		k = "-" + k
		cmd.append(k)
		cmd.append(v)

	# Add input_file
	cmd.append(input_file)
	print(cmd)

	with open(output_file, "w") as fout:
		subprocess.run(cmd, stdout=fout)


def get_treebest_nj_runtime(num_seq, msa_np1_fn_list):
	container = defaultdict(list)
	for n, msa_fn in zip(num_seq, msa_np1_fn_list):
		print(f"Number of sequence: {n}")
		container["num"].append(n)

		output_file = msa_fn.replace(".ali", ".nj.nhx")
		c = msa_fn.replace(".np1.ali", ".treefile")
		start = time.time()
		treebest(command="nj", input_file=msa_fn, output_file=output_file, c=c)
		tree_best = time.time()
		container["treebest_nj_time"].append(tree_best - start)
	runtime = pd.DataFrame(container)
	return runtime


def get_treebest_phyml_runtime(num_seq, msa_np1_fn_list):
	container = defaultdict(list)
	for n, msa_fn in zip(num_seq, msa_np1_fn_list):
		print(f"Number of sequence: {n}")
		container["num"].append(n)

		output_file = msa_fn.replace(".ali", ".phyml.nhx")
		C = msa_fn.replace(".np1.ali", ".treefile")
		start = time.time()
		treebest(command="phyml", input_file=msa_fn, output_file=output_file, C=C)
		tree_best = time.time()
		container["treebest_phyml_time"].append(tree_best - start)
	runtime = pd.DataFrame(container)
	return runtime


def plot_treebest_runtime_vs_num_seq(treebest_runtime, fig_fn):

	plt.close()
	fig, ax = plt.subplots(2,3, figsize=(9,6))
	ax[0,0].plot("num", "tree_time", "go", data=treebest_runtime)
	ax[0,0].set_title("iqtree")
	ax[0,0].set_xlabel("Num of seq")
	ax[0,0].set_ylabel("Runtime")

	ax[1,0].plot("num", "tree_time", "go", data=treebest_runtime)
	ax[1,0].set_title("iqtree (log-log)")
	ax[1,0].set_xscale("log")
	ax[1,0].set_yscale("log")
	ax[1,0].set_xlabel("log(num of seq)")
	ax[1,0].set_ylabel("log(untime)")

	ax[0,1].plot("num", "treebest_nj_time", "go", data=treebest_runtime)
	ax[0,1].set_title("Treebest NJ")
	ax[0,1].set_xlabel("Num of seq")
	ax[0,1].set_ylabel("Runtime")

	ax[1,1].plot("num", "treebest_nj_time", "go", data=treebest_runtime)
	ax[1,1].set_xscale("log")
	ax[1,1].set_yscale("log")
	ax[1,1].set_title("Treebest NJ (log-log)")
	ax[1,1].set_xlabel("log(num of seq)")
	ax[1,1].set_ylabel("log(runtime)")

	ax[0,2].plot("num", "treebest_phyml_time", "go", data=treebest_runtime)
	ax[0,2].set_title("Treebest PhyML")
	ax[0,2].set_xlabel("Num of seq")
	ax[0,2].set_ylabel("Runtime")

	ax[1,2].plot("num", "treebest_phyml_time", "go", data=treebest_runtime)
	ax[1,2].set_xscale("log")
	ax[1,2].set_yscale("log")
	ax[1,2].set_title("Treebest PhyML (log-log)")
	ax[1,2].set_xlabel("log(num of seq)")
	ax[1,2].set_ylabel("log(runtime)")

	fig.tight_layout()

	fig.savefig(fig_fn)


def main():
	print("Sample FASTA records.")
	fasta_fn_list = sample_fasta(gisaid_fasta_fn, out_dir, num_seq)

	print("Collect de novo runtime statistics.")
	print("Multiple sequence alignment.")
	denovo_msa_runtime, denovo_msa_fn_list = get_msa_runtime( \
		mode="denovo", num_seq=num_seq, fasta_fn_list=fasta_fn_list)
	with open(f"{out_dir}/denovo_msa_fn.pkl", "wb") as p: 
		pickle.dump(denovo_msa_fn_list, p)

	print("Phylogenetic tree with IQ-Tree.")
	denovo_iqtree_runtime, denovo_iqtree_fn_list = get_tree_runtime( \
		denovo_msa_fn_list, num_seq, software="iqtree")
	with open(f"{out_dir}/denovo_iqtree_fn.pkl", "wb") as p:
		pickle.dump(denovo_iqtree_fn_list, p)

	print("Phylogenetic tree with FastTree.")
	denovo_fastTree_runtime, denovo_fastTree_fn_list = get_tree_runtime( \
		denovo_msa_fn_list, num_seq, software="FastTree")
	with open(f"{out_dir}/denovo_fastTree_fn.pkl", "wb") as p:
		pickle.dump(denovo_fastTree_fn_list, p)

	print("Save runtime statistics.")
	denovo_runtime = pd.merge(denovo_msa_runtime, denovo_iqtree_runtime, on="num")
	denovo_runtime = pd.merge(denovo_runtime, denovo_fastTree_runtime, on="num")
	denovo_runtime_out_fn = f"{out_dir}/denovo_runtime.tsv"
	denovo_runtime.to_csv(denovo_runtime_out_fn, sep="\t", index=False)
	
	print("Plot runtime statistics.")
	denovo_fig_fn = f"{out_dir}/denovo_runtime_vs_num_seq.png"
	plot_runtime_vs_num_seq(denovo_runtime, denovo_fig_fn)

	np1_runtime, msa_np1_fn_list = get_np1_runtime(msa_fn_list, num_seq, beijing_fasta_fn)
	np1_runtime_out_fn = f"{out_dir}/np1_msa_runtime.tsv"
	np1_runtime.to_csv(np1_runtime_out_fn, sep="\t", index=False)
	msa_np1_fn_list_pkl = f"{out_dir}/np1_msa_fn.pkl"
	with open(msa_np1_fn_list_pkl, "wb") as p: 
		pickle.dump(msa_np1_fn_list, p)

	np1_fig_fn = f"{out_dir}/np1_runtime_vs_num_seq.png"
	plot_runtime_vs_num_seq(np1_runtime, np1_fig_fn)

	treebest_nj_runtime = get_treebest_nj_runtime(num_seq, msa_np1_fn_list)
	treebest_runtime = pd.merge(np1_runtime, treebest_nj_runtime)

	treebest_phyml_runtime = get_treebest_phyml_runtime(num_seq, msa_np1_fn_list)
	treebest_runtime = pd.merge(treebest_runtime, treebest_phyml_runtime)
	treebest_runtime_out_fn = f"{out_dir}/treebest_runtime.tsv"
	treebest_runtime.to_csv(treebest_runtime_out_fn, sep="\t", index=False)

	treebest_fig_fn = f"{out_dir}/treebest_runtime_vs_num_seq.png"
	plot_treebest_runtime_vs_num_seq(treebest_runtime, treebest_fig_fn)

if __name__ == "__main__":
	main()
