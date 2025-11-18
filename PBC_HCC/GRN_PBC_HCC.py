#region Stem-like Cells
# conda activate GRN
import os
import glob
import pickle
import pandas as pd
import numpy as np

# from dask.diagnostics import ProgressBar
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

os.chdir("/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/TFs/")
# expr_matrix: cells x genes (pandas DataFrame)
# tf_names: list of TF gene symbols
exp_matrix = pd.read_table("/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/Table/PBC_HCC_SM_normalize.txt", index_col=0)
tf_names = pd.read_table("/diazlab/data3/.abhinav/.immune/resources/GRN/TFs/hg38_TF_pyscenic_list.txt", header = None).iloc[:,0].to_list()
tf_names = ["BACH2", "BCL11B", "KLF12", "ZBTB20", "FOXN3","TCF7","LEF1"]
exp_matrix_trans = exp_matrix.T
tf_names_2 = list(set(tf_names) & set(exp_matrix_trans.columns))
network = grnboost2(expression_data=exp_matrix_trans, tf_names=tf_names_2)
savedir = "/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/TFs/unbiased/"
os.makedirs(savedir, exist_ok = True)
out_file = os.path.join(savedir,"networks.csv")
# exp_file = "/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/TFs/exp_file.csv"
network.to_csv(out_file, sep='\t', index=False, header=True)
# exp_matrix_trans.to_csv(exp_file, sep='\t', index=True, header=True)

# In addition, the transcription factor is added to the module and modules that have less than 20 genes are removed.
out_file = "/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/TFs/networks.csv"
exp_file = "/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/TFs/exp_file.csv" 
network = pd.read_table(out_file)
ex_matrix = pd.read_table(exp_file, index_col = 0)
# if you forgot to put the index_col =0 
# ex_matrix.index=ex_matrix.iloc[:,0] 
# ex_matrix=ex_matrix.iloc[:,1:]
# modules = list(modules_from_adjacencies(network, ex_matrix))   ### This make the regulons for different transcription factors
modules = list(modules_from_adjacencies(network, exp_matrix_trans))   ### This make the regulons for different transcription factors
# Calculate a list of enriched motifs and the corresponding target genes for all modules.
DATABASE_FOLDER="/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/TFs/"
DATABASES_GLOB=os.path.join(DATABASE_FOLDER, "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
MOTIF_ANNOTATIONS_FNAME=os.path.join(DATABASE_FOLDER,"motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl")
REGULONS_FNAME = os.path.join(DATABASE_FOLDER, "regulons.p")
MOTIFS_FNAME = os.path.join(DATABASE_FOLDER, "motifs.csv")

db_fnames = glob.glob(DATABASES_GLOB)
def name(fname):
    return os.path.splitext(os.path.basename(fname))[0]

dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]

with ProgressBar():
    df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)

# Create regulons from this table of enriched motifs.
regulons = df2regulons(df)

# Save the enriched motifs and the discovered regulons to disk.
df.to_csv(MOTIFS_FNAME)
with open(REGULONS_FNAME, "wb") as f:
    pickle.dump(regulons, f)

# We characterize the different cells in a single-cell transcriptomics experiment via the enrichment of the
# previously discovered regulons. Enrichment of a regulon is measured as the Area Under the recovery Curve (AUC) 
# of the genes that define this regulon.
auc_mtx = aucell(ex_matrix, regulons, num_workers=4)
import seaborn as sns

# Create the clustermap and assign it to a variable 'g'
g = sns.clustermap(auc_mtx, figsize=(8, 8))

# Save the plot using the savefig() method of the ClusterGrid object
os.chdir("/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/TFs/")
g.savefig("clustermap.pdf")

import networkx as nx
import matplotlib.pyplot as plt

# Example for one regulon
regulon = regulons[0]

G = nx.DiGraph()  # directed graph

# Add TF node
G.add_node(regulon.transcription_factor, type='TF')

# Add target genes with edge weights
for gene, weight in regulon.gene2weight.items():
    G.add_node(gene, type='target')
    G.add_edge(regulon.transcription_factor, gene, weight=weight)

# Draw graph
pos = nx.spring_layout(G, seed=42)
weights = [G[u][v]['weight'] for u,v in G.edges()]

# Node colors by type
node_color = ['skyblue' if G.nodes[n]['type']=='TF' else 'lightgreen' for n in G.nodes()]

plt.figure(figsize=(10,10))
nx.draw(G, pos, with_labels=True, node_color=node_color, width=[w/5 for w in weights],
        arrowsize=20, font_size=10)
plt.show()

plt.savefig("BACH2_network.pdf", bbox_inches='tight')

import networkx as nx
import matplotlib.pyplot as plt

# Create a directed graph
G = nx.DiGraph()

# Loop over all regulons
for regulon in regulons:
    tf = regulon.transcription_factor
    G.add_node(tf, type='TF')  # TF node
    for gene, weight in regulon.gene2weight.items():
        G.add_node(gene, type='target')
        G.add_edge(tf, gene, weight=weight)

# Layout
pos = nx.spring_layout(G, seed=42)  # deterministic layout
weights = [G[u][v]['weight'] for u,v in G.edges()]

# Node color: TF = skyblue, Target = lightgreen
node_color = ['skyblue' if G.nodes[n]['type']=='TF' else 'lightgreen' for n in G.nodes()]

# Node size by degree (optional)
node_size = [300 + 50*G.degree(n) for n in G.nodes()]

plt.figure(figsize=(20, 20))
nx.draw(
    G, pos, with_labels=True, node_color=node_color,
    width=[w/5 for w in weights], arrowsize=20,
    font_size=8, node_size=node_size
)

# Save figure
plt.savefig("all_regulons_network.png", dpi=300, bbox_inches='tight')
plt.savefig("all_regulons_network.pdf", bbox_inches='tight')
plt.show()

## making for each regulon
import networkx as nx
import matplotlib.pyplot as plt
import re

# Function to sanitize regulon names for filenames
def clean_name(name):
    return re.sub(r'[^A-Za-z0-9_-]', '_', name)

# Loop over all regulons
for regulon in regulons:
    G = nx.DiGraph()
    tf = regulon.transcription_factor
    G.add_node(tf, type='TF')  # TF node
    for gene, weight in regulon.gene2weight.items():
        G.add_node(gene, type='target')
        G.add_edge(tf, gene, weight=weight)
    # Layout
    pos = nx.spring_layout(G, seed=42)
    weights = [G[u][v]['weight'] for u,v in G.edges()]
    node_color = ['skyblue' if G.nodes[n]['type']=='TF' else 'lightgreen' for n in G.nodes()]
    node_size = [300 + 50*G.degree(n) for n in G.nodes()]
    plt.figure(figsize=(12, 12))
    nx.draw(
        G, pos, with_labels=True, node_color=node_color,
        width=[w/5 for w in weights], arrowsize=20,
        font_size=10, node_size=node_size
    )
    # Save figure
    filename = f"{clean_name(regulon.name)}_network.pdf"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()  # Close figure to free memory
    print(f"Saved figure for regulon: {regulon.name} -> {filename}")

#endregion

#region All data
# conda activate GRN
import os
import glob
import pickle
import pandas as pd
import numpy as np

# from dask.diagnostics import ProgressBar
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

savedir = "/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/TFs/all/"
os.makedirs(savedir, exist_ok = True)
# expr_matrix: cells x genes (pandas DataFrame)
# tf_names: list of TF gene symbols
exp_matrix = pd.read_table("/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/Table/PBC_HCC_normalize.txt", index_col=0)
tf_names = pd.read_table("/diazlab/data3/.abhinav/.immune/resources/GRN/TFs/hg38_TF_pyscenic_list.txt", header = None).iloc[:,0].to_list()
exp_matrix_trans = exp_matrix.T
tf_names_2 = list(set(tf_names) & set(exp_matrix_trans.columns))
network = grnboost2(expression_data=exp_matrix_trans, tf_names=tf_names_2)
savedir = "/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/TFs/unbiased/"
os.makedirs(savedir, exist_ok = True)
out_file = os.path.join(savedir,"networks.csv")
exp_file = os.path.join(savedir,"exp_file.csv")
network.to_csv(out_file, sep='\t', index=False, header=True)
exp_matrix_trans.to_csv(exp_file, sep='\t', index=True, header=True)

# In addition, the transcription factor is added to the module and modules that have less than 20 genes are removed.
out_file = "/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/TFs/unbiased/networks.csv"
exp_file = "/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/TFs/exp_file.csv" 
network = pd.read_table(out_file)
exp_matrix = pd.read_table(exp_file, index_col = 0)
# if you forgot to put the index_col =0 
# exp_matrix.index=exp_matrix.iloc[:,0] 
# exp_matrix=exp_matrix.iloc[:,1:]
modules = list(modules_from_adjacencies(network, ex_matrix))   ### This make the regulons for different transcription factors
# Calculate a list of enriched motifs and the corresponding target genes for all modules.
DATABASE_FOLDER="/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/TFs/"
DATABASES_GLOB=os.path.join(DATABASE_FOLDER, "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
MOTIF_ANNOTATIONS_FNAME=os.path.join(DATABASE_FOLDER,"motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl")
REGULONS_FNAME = os.path.join(DATABASE_FOLDER, "regulons.p")
MOTIFS_FNAME = os.path.join(DATABASE_FOLDER, "motifs.csv")

db_fnames = glob.glob(DATABASES_GLOB)
def name(fname):
    return os.path.splitext(os.path.basename(fname))[0]

dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]

with ProgressBar():
    df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)

# Create regulons from this table of enriched motifs.
regulons = df2regulons(df)

# Save the enriched motifs and the discovered regulons to disk.
df.to_csv(MOTIFS_FNAME)
with open(REGULONS_FNAME, "wb") as f:
    pickle.dump(regulons, f)

# We characterize the different cells in a single-cell transcriptomics experiment via the enrichment of the
# previously discovered regulons. Enrichment of a regulon is measured as the Area Under the recovery Curve (AUC) 
# of the genes that define this regulon.
auc_mtx = aucell(ex_matrix, regulons, num_workers=4)
import seaborn as sns

# Create the clustermap and assign it to a variable 'g'
g = sns.clustermap(auc_mtx, figsize=(8, 8))

# Save the plot using the savefig() method of the ClusterGrid object
os.chdir("/diazlab/data3/.abhinav/.immune/Kyoto/PBC_HCC/TFs/")
g.savefig("clustermap.pdf")

import networkx as nx
import matplotlib.pyplot as plt

# Example for one regulon
regulon = regulons[0]

G = nx.DiGraph()  # directed graph

# Add TF node
G.add_node(regulon.transcription_factor, type='TF')

# Add target genes with edge weights
for gene, weight in regulon.gene2weight.items():
    G.add_node(gene, type='target')
    G.add_edge(regulon.transcription_factor, gene, weight=weight)

# Draw graph
pos = nx.spring_layout(G, seed=42)
weights = [G[u][v]['weight'] for u,v in G.edges()]

# Node colors by type
node_color = ['skyblue' if G.nodes[n]['type']=='TF' else 'lightgreen' for n in G.nodes()]

plt.figure(figsize=(10,10))
nx.draw(G, pos, with_labels=True, node_color=node_color, width=[w/5 for w in weights],
        arrowsize=20, font_size=10)
plt.show()

plt.savefig("BACH2_network.pdf", bbox_inches='tight')

import networkx as nx
import matplotlib.pyplot as plt

# Create a directed graph
G = nx.DiGraph()

# Loop over all regulons
for regulon in regulons:
    tf = regulon.transcription_factor
    G.add_node(tf, type='TF')  # TF node
    for gene, weight in regulon.gene2weight.items():
        G.add_node(gene, type='target')
        G.add_edge(tf, gene, weight=weight)

# Layout
pos = nx.spring_layout(G, seed=42)  # deterministic layout
weights = [G[u][v]['weight'] for u,v in G.edges()]

# Node color: TF = skyblue, Target = lightgreen
node_color = ['skyblue' if G.nodes[n]['type']=='TF' else 'lightgreen' for n in G.nodes()]

# Node size by degree (optional)
node_size = [300 + 50*G.degree(n) for n in G.nodes()]

plt.figure(figsize=(20, 20))
nx.draw(
    G, pos, with_labels=True, node_color=node_color,
    width=[w/5 for w in weights], arrowsize=20,
    font_size=8, node_size=node_size
)

# Save figure
plt.savefig("all_regulons_network.png", dpi=300, bbox_inches='tight')
plt.savefig("all_regulons_network.pdf", bbox_inches='tight')
plt.show()

## making for each regulon
import networkx as nx
import matplotlib.pyplot as plt
import re

# Function to sanitize regulon names for filenames
def clean_name(name):
    return re.sub(r'[^A-Za-z0-9_-]', '_', name)

# Loop over all regulons
for regulon in regulons:
    G = nx.DiGraph()
    tf = regulon.transcription_factor
    G.add_node(tf, type='TF')  # TF node
    for gene, weight in regulon.gene2weight.items():
        G.add_node(gene, type='target')
        G.add_edge(tf, gene, weight=weight)
    # Layout
    pos = nx.spring_layout(G, seed=42)
    weights = [G[u][v]['weight'] for u,v in G.edges()]
    node_color = ['skyblue' if G.nodes[n]['type']=='TF' else 'lightgreen' for n in G.nodes()]
    node_size = [300 + 50*G.degree(n) for n in G.nodes()]
    plt.figure(figsize=(12, 12))
    nx.draw(
        G, pos, with_labels=True, node_color=node_color,
        width=[w/5 for w in weights], arrowsize=20,
        font_size=10, node_size=node_size
    )
    # Save figure
    filename = f"{clean_name(regulon.name)}_network.pdf"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()  # Close figure to free memory
    print(f"Saved figure for regulon: {regulon.name} -> {filename}")
