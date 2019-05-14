#! /usr/bin/env python3

import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import velocyto as vcy
import pandas as pd
import loompy
from sklearn.svm import SVR
from sklearn.linear_model import LinearRegression
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.interpolate import interp1d

sample_type = "Sorted"
sample_id = "Russ_GFP_Sort"
out_dir = "/Users/kriemo/Projects/sc_repos/sebeta/results/2018-12-19_figs/velocity/"

velo_data_dir = "/Users/kriemo/Projects/sc_repos/sebeta/data/velocyto"

loom_fn = os.path.join(velo_data_dir, sample_id, sample_id + ".loom")

metadata_fn = os.path.join("/Users/kriemo/Projects/sc_repos/sebeta/results/2018-12-19_figs", "GEO", sample_type + "_metadata.tsv.gz")
clustering_col_id = "res.0.5"
clustering_col_id = "cell_type"
hdf5_out_fn = sample_id + ".hdf5"
pdf_out_dir = sample_id

vlm = vcy.VelocytoLoom(loom_fn)
mdata = pd.read_csv(metadata_fn, sep="\t")

# get cell ids to match loom object
new_ids  = [sample_id + ":" + x + "x" for x in mdata["cell"] ]
mdata = mdata.assign(new_id = new_ids)

mdata = mdata[mdata.new_id.isin(list(vlm.ca["CellID"]))]

# reorder cell ids to match loom object
cids = pd.DataFrame({'new_id' : vlm.ca["CellID"]})
mdata = pd.merge(cids, mdata, how = 'left', on = 'new_id')
mdata = mdata.dropna()

#only keep cells found in seurat data
keep_idx = [x in list(mdata["new_id"]) for x in list(vlm.ca["CellID"])]
vlm.filter_cells(bool_array=keep_idx)

#add cluster annotations

col_df = pd.read_csv(os.path.join("/Users/kriemo/Projects/sc_repos/sebeta/results/2018-12-19_figs", sample_type + "_colorMap_with_commas.tsv"),
                     header = None, sep = "\t")

from matplotlib.colors import hex2color

cols = {}
for i in range(col_df.shape[0]):
  key =  col_df.iloc[i, 0]
  cols[key] = hex2color(col_df.iloc[i, 1])


vlm.set_clusters(np.array([str(x) for x in mdata[clustering_col_id]]), cluster_colors_dict = cols)

#add tSNE projections
vlm.ca["TSNE1"] = np.array(mdata["tSNE_1"])
vlm.ca["TSNE2"] = np.array(mdata["tSNE_2"])
vlm.ts = np.column_stack([vlm.ca["TSNE1"], vlm.ca["TSNE2"]])


vlm.normalize("S", size=True, log=True)


vlm.filter_cells(bool_array=vlm.initial_Ucell_size > np.percentile(vlm.initial_Ucell_size, 0.5))

vlm.score_detection_levels(min_expr_counts=30, min_cells_express=20)
vlm.filter_genes(by_detection_levels=True)
vlm.score_cv_vs_mean(3000, plot=True, max_expr_avg=35)

vlm.filter_genes(by_cv_vs_mean=True)

vlm._normalize_S(relative_size=vlm.S.sum(0),
             target_size=vlm.S.sum(0).mean())
vlm._normalize_U(relative_size=vlm.U.sum(0),
             target_size=vlm.U.sum(0).mean())
vlm.perform_PCA()
vlm.knn_imputation(n_pca_dims=20,
                   k = int(len(vlm.ca["CellID"]) * 0.075),
                   balanced=True,
                   b_sight=3000, b_maxl=1500, n_jobs=6)
vlm.fit_gammas()
vlm.predict_U()
vlm.calculate_velocity()
vlm.calculate_shift(assumption="constant_velocity")
vlm.extrapolate_cell_at_t(delta_t=1.)


vlm.estimate_transition_prob(hidim="Sx_sz", embed="ts", transform="sqrt",
                             psc=1, n_neighbors= int(vlm.S.shape[1] / 5),
                             knn_random=True, sampled_fraction=0.5)
vlm.calculate_embedding_shift(sigma_corr = 0.05, expression_scaling=True)
vlm.calculate_grid_arrows(smooth = 0.5,
                          steps=(30, 30),
                          n_neighbors=int(vlm.S.shape[1] * 0.10))


def despline():
    ax1 = plt.gca()
    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax1.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')

def minimal_xticks(start, end):
    end_ = np.around(end, -int(np.log10(end))+1)
    xlims = np.linspace(start, end_, 5)
    xlims_tx = [""]*len(xlims)
    xlims_tx[0], xlims_tx[-1] = f"{xlims[0]:.0f}", f"{xlims[-1]:.02f}"
    plt.xticks(xlims, xlims_tx)


def minimal_yticks(start, end):
    end_ = np.around(end, -int(np.log10(end))+1)
    ylims = np.linspace(start, end_, 5)
    ylims_tx = [""]*len(ylims)
    ylims_tx[0], ylims_tx[-1] = f"{ylims[0]:.0f}", f"{ylims[-1]:.02f}"
    plt.yticks(ylims, ylims_tx)



plt.figure(None,(14,14))
quiver_scale = 60
plt.scatter(vlm.embedding[:, 0], vlm.embedding[:, 1],
            c="0.8", alpha=0.2, s=10, edgecolor="")

ix_choice = np.random.choice(vlm.embedding.shape[0], size=int(vlm.embedding.shape[0]/1.), replace=False)
plt.scatter(vlm.embedding[ix_choice, 0], vlm.embedding[ix_choice, 1],
            c="0.8", alpha=0.4, s=10, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)

quiver_kwargs=dict(headaxislength=7, headlength=11, headwidth=8,linewidths=0.25, width=0.00045,edgecolors="k", color=vlm.colorandum[ix_choice], alpha=1)
plt.quiver(vlm.embedding[ix_choice, 0], vlm.embedding[ix_choice, 1],
           vlm.delta_embedding[ix_choice, 0], vlm.delta_embedding[ix_choice, 1],
           scale=quiver_scale, **quiver_kwargs)

plt.axis("off")

plt.savefig(os.path.join(out_dir, sample_id + "_full_arrows.pdf"))

# np.array(mdata[clustering_col_id])

#vlm.set_clusters(np.array([str(x) for x in mdata[clustering_col_id]]),
#                 colormap = plt.cm.get_cmap("Paired"))
#vlm

plt.figure(None,(20,10))
vlm.plot_grid_arrows(quiver_scale="auto",
                     headaxislength=2.75,
                     headlength=5,
                     headwidth=4.8,
                     minlength=1.5,
                     plot_random=True,
                     scale_type="relative")
              #       scatter_kwargs_dict = {'alpha' : 1})

plt.savefig(os.path.join(out_dir, sample_id + "_vectorfield.pdf"))

plt.figure(None,(7,7))
vlm.plot_grid_arrows(quiver_scale="auto",
                     headaxislength=2.75,
                     headlength=5,
                     headwidth=4.8,
                     minlength=1.5,
                     plot_random=False,
                     scale_type="relative",
                     scatter_kwargs_dict = {'alpha' : 1,
                                            's' : 2.5,
                                            'lw' : 1})
plt.savefig(os.path.join(out_dir, sample_id + "_vectorfield_pub.pdf"), transparent = True)


#genes = ["INS","IAPP", "CD9", "PCP4", "ERO1B",
#         "OVOS2", "TPH1", "FEV", "STC1",
#         "TTR", "PYY", "ENTPD3"]
#plt.figure(None, (17,24), dpi=80)
#gs = plt.GridSpec(10,6)
#for i, gn in enumerate(genes):
#    ax = plt.subplot(gs[i*3])
#    try:
#        ix=np.where(vlm.ra["Gene"] == gn)[0][0]
#    except:
#        continue
#    vcy.scatter_viz(vlm.Sx_sz[ix,:], vlm.Ux_sz[ix,:], c=vlm.colorandum, s=5, alpha=0.4, rasterized=True)
#    plt.title(gn)
#    xnew = np.linspace(0,vlm.Sx[ix,:].max())
#    plt.plot(xnew, vlm.gammas[ix] * xnew + vlm.q[ix], c="k")
#    plt.ylim(0, np.max(vlm.Ux_sz[ix,:])*1.02)
#    plt.xlim(0, np.max(vlm.Sx_sz[ix,:])*1.02)
#    minimal_yticks(0, np.max(vlm.Ux_sz[ix,:])*1.02)
#    minimal_xticks(0, np.max(vlm.Sx_sz[ix,:])*1.02)
#    despline()
#
#    vlm.plot_velocity_as_color(gene_name=gn, gs=gs[i*3+1], s=3, rasterized=True)
#
#    vlm.plot_expression_as_color(gene_name=gn, gs=gs[i*3+2], s=3, rasterized=True)

#plt.tight_layout()
#plt.savefig(sample_id + "_genes.pdf")


## Markov Chain analysis

# Sample uniformly the points to avoid density driven effects - Should reimplement as a method
steps = 100, 100
grs = []
for dim_i in range(vlm.embedding.shape[1]):
    m, M = np.min(vlm.embedding[:, dim_i]), np.max(vlm.embedding[:, dim_i])
    m = m - 0.025 * np.abs(M - m)
    M = M + 0.025 * np.abs(M - m)
    gr = np.linspace(m, M, steps[dim_i])
    grs.append(gr)

meshes_tuple = np.meshgrid(*grs)
gridpoints_coordinates = np.vstack([i.flat for i in meshes_tuple]).T

from sklearn.neighbors import NearestNeighbors
nn = NearestNeighbors()
nn.fit(vlm.embedding)
dist, ixs = nn.kneighbors(gridpoints_coordinates, 1)

diag_step_dist = np.sqrt((meshes_tuple[0][0,0] - meshes_tuple[0][0,1])**2 + (meshes_tuple[1][0,0] - meshes_tuple[1][1,0])**2)
min_dist = diag_step_dist / 2
ixs = ixs[dist < min_dist]
gridpoints_coordinates = gridpoints_coordinates[dist.flat[:]<min_dist,:]
dist = dist[dist < min_dist]

ixs = np.unique(ixs)

plt.figure(None,(8,8))
vcy.scatter_viz(vlm.embedding[ixs, 0], vlm.embedding[ixs, 1],
                c=vlm.colorandum[ixs], alpha=1, s=30, lw=0.4,
                edgecolor="0.4")


vlm.prepare_markov(sigma_D=diag_step_dist,
                   sigma_W=diag_step_dist/2.,
                   direction='forward', cells_ixs=ixs)

vlm.run_markov(starting_p=np.ones(len(ixs)), n_steps=2500)


diffused_n = vlm.diffused - np.percentile(vlm.diffused, 3)
diffused_n /= np.percentile(diffused_n, 97)
diffused_n = np.clip(diffused_n, 0, 1)

plt.figure(None,(7,7))
vcy.scatter_viz(vlm.embedding[ixs, 0], vlm.embedding[ixs, 1],
                c=diffused_n, alpha=0.5, s=50, lw=0.,
                edgecolor="", cmap="viridis_r", rasterized=False)
plt.axis("off")
plt.savefig(os.path.join(out_dir, sample_id + "_endpoint_distr.pdf"))


plt.figure(None,(7,7))
vcy.scatter_viz(vlm.embedding[ixs, 0], vlm.embedding[ixs, 1],
                c=diffused_n, alpha=0.5, s=25, lw=0.,
                edgecolor="", cmap="viridis_r", rasterized=False)
plt.axis("off")
plt.savefig(os.path.join(out_dir, sample_id + "_endpoint_distr_pub.pdf"), transparent = True)

out_df = pd.DataFrame({"CellID" : vlm.ca["CellID"][ixs],
                       "endpoints" : diffused_n})


vlm.prepare_markov(sigma_D=diag_step_dist,
                   sigma_W=diag_step_dist/2.,
                   direction='backwards', cells_ixs=ixs)
vlm.run_markov(starting_p=np.ones(len(ixs)), n_steps=2500)

diffused_n = vlm.diffused - np.percentile(vlm.diffused, 3)
diffused_n /= np.percentile(diffused_n, 97)
diffused_n = np.clip(diffused_n, 0, 1)

plt.figure(None,(7,7))
vcy.scatter_viz(vlm.embedding[ixs, 0], vlm.embedding[ixs, 1],
                c=diffused_n, alpha=0.5, s=50, lw=0.,
                edgecolor="", cmap="viridis_r", rasterized=False)
plt.axis("off")
plt.savefig(os.path.join(out_dir, sample_id + "_startpoint_distr.pdf"))

plt.figure(None,(7,7))
vcy.scatter_viz(vlm.embedding[ixs, 0], vlm.embedding[ixs, 1],
                c=diffused_n, alpha=0.5, s=25, lw=0.,
                edgecolor="", cmap="viridis_r", rasterized=False)
plt.axis("off")
plt.savefig(os.path.join(out_dir, sample_id + "_startpoint_distr_pub.pdf"), transparent = True)


out_df["startpoints"] = diffused_n

out_df.to_csv(os.path.join(out_dir, sample_id + "_markov_points.tsv"),  sep = "\t", index = False)

#vlm.to_hdf5(hdf5_out_fn)
