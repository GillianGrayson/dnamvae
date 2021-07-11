import pandas as pd
from src.dnam.routines.manifest import get_manifest
from src.dnam.routines.plot.save import save_figure
from src.dnam.routines.plot.scatter import add_scatter_trace
from src.dnam.routines.plot.layout import add_layout
from src.dnam.routines.datasets_features import *
import os
import plotly.graph_objects as go


dataset = "GSE168739"
platform = "GPL21145"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

is_recalc = True

age_pair = tuple([x.replace(' ','_') for x in get_age_pair(dataset)])
sex_pair = tuple([x.replace(' ','_') for x in get_sex_pair(dataset)])
sex_vals_pairs = get_sex_vals_pairs(dataset)

pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno_xtd.pkl")
pheno.columns = pheno.columns.str.replace(' ','_')
betas = pd.read_pickle(f"{path}/{platform}/{dataset}/betas.pkl")

targets = {
    "cg17178900": ["EPICOVID-098", "EPICOVID-236"],
    "cg13452062": ["EPICOVID-254", "EPICOVID-348"],
    "cg24795173": ["EPICOVID-198", "EPICOVID-406"],
    "cg02872426": ["EPICOVID-328", "EPICOVID-030"],
    "cg08309069": ["EPICOVID-309", "EPICOVID-377"],
    "cg04736673": ["EPICOVID-192", "EPICOVID-027"],
    "cg14859874": ["EPICOVID-168", "EPICOVID-125"],
    "cg07796016": ["EPICOVID-281", "EPICOVID-243"]
}
betas = betas[list(targets.keys())]

df = pd.merge(pheno, betas, left_index=True, right_index=True)
df_pass = df[df['Group'].notnull()]

manifest = get_manifest(platform)

save_path = f"{path}/{platform}/{dataset}/EWAS/status_extraction/figs"
if not os.path.exists(save_path):
    os.makedirs(save_path)

cmb = go.Figure()
for cpg_id, cpg in enumerate(targets.keys()):

    if cpg_id == 0:
        showlegend = True
    else:
        showlegend = False

    cmb.add_trace(go.Violin(x=[cpg] * len(df[cpg].values),
                            y=df[cpg].values,
                            points=False,
                            showlegend=False,
                            line=dict(color='grey'),
                            fillcolor='grey',
                            opacity=0.6
                            ))
    cmb.add_trace(go.Violin(x=[cpg] * len(df.loc[(df['Group'] == "G2/G1") & (~df['description'].isin(targets[cpg])), age_pair[0]].values),
                            y=df.loc[(df['Group'] == "G2/G1") & (~df['description'].isin(targets[cpg])), cpg].values,
                            points='all',
                            showlegend=False,
                            jitter=0.5,
                            pointpos=0.0,
                            marker=dict(
                                color='red',
                                size=4,
                                opacity=0.6,
                                line=dict(color='black', width=1)
                            ),
                            line=dict(color='rgba(0,0,0,0)'),
                            fillcolor='rgba(0,0,0,0)'
                            ))
    cmb.add_trace(go.Violin(x=[cpg] * len(df.loc[(df['Group'] == "G2/G1") & (df['description'].isin(targets[cpg])), age_pair[0]].values),
                            y=df.loc[(df['Group'] == "G2/G1") & (df['description'].isin(targets[cpg])), cpg].values,
                            name="G2/G1",
                            points='all',
                            showlegend=showlegend,
                            jitter=0.5,
                            pointpos=0.0,
                            marker=dict(
                                color='red',
                                size=8,
                                symbol='square',
                                opacity=0.6,
                                line=dict(color='black', width=1)
                            ),
                            line=dict(color='rgba(0,0,0,0)'),
                            fillcolor='rgba(0,0,0,0)'
                            ))
    cmb.add_trace(go.Violin(x=[cpg] * len(df.loc[(df['Group'] == "G3") & (~df['description'].isin(targets[cpg])), age_pair[0]].values),
                            y=df.loc[(df['Group'] == "G3") & (~df['description'].isin(targets[cpg])), cpg].values,
                            points='all',
                            showlegend=False,
                            jitter=0.5,
                            pointpos=0.0,
                            marker=dict(
                                color='blue',
                                size=4,
                                opacity=0.6,
                                line=dict(color='black', width=1)
                            ),
                            line=dict(color='rgba(0,0,0,0)'),
                            fillcolor='rgba(0,0,0,0)'
                            ))
    cmb.add_trace(go.Violin(x=[cpg] * len(df.loc[(df['Group'] == "G3") & (df['description'].isin(targets[cpg])), age_pair[0]].values),
                            y=df.loc[(df['Group'] == "G3") & (df['description'].isin(targets[cpg])), cpg].values,
                            name="G3",
                            points='all',
                            showlegend=showlegend,
                            jitter=0.5,
                            pointpos=0.0,
                            marker=dict(
                                color='blue',
                                size=8,
                                symbol='square',
                                opacity=0.6,
                                line=dict(color='black', width=1)
                            ),
                            line=dict(color='rgba(0,0,0,0)'),
                            fillcolor='rgba(0,0,0,0)'
                            ))

    fig = go.Figure()
    add_scatter_trace(fig, df.loc[df['Group'].isnull(), age_pair[0]].values, df.loc[df['Group'].isnull(), cpg].values, "No Info")
    add_scatter_trace(
        fig,
        df.loc[(df['Group'] == "G2/G1") & (~df['description'].isin(targets[cpg])), age_pair[0]].values,
        df.loc[(df['Group'] == "G2/G1") & (~df['description'].isin(targets[cpg])), cpg].values,
        "G2/G1"
    )
    add_scatter_trace(
        fig,
        df.loc[(df['Group'] == "G2/G1") & (df['description'].isin(targets[cpg])), age_pair[0]].values,
        df.loc[(df['Group'] == "G2/G1") & (df['description'].isin(targets[cpg])), cpg].values,
        "",
        size=15
    )
    add_scatter_trace(
        fig,
        df.loc[(df['Group'] == "G3") & (~df['description'].isin(targets[cpg])), age_pair[0]].values,
        df.loc[(df['Group'] == "G3") & (~df['description'].isin(targets[cpg])), cpg].values,
        "G3"
    )
    add_scatter_trace(
        fig,
        df.loc[(df['Group'] == "G3") & (df['description'].isin(targets[cpg])), age_pair[0]].values,
        df.loc[(df['Group'] == "G3") & (df['description'].isin(targets[cpg])), cpg].values,
        "",
        size=15
    )
    add_layout(fig, age_pair[1], 'Methylation Level', f"{cpg} ({manifest.loc[cpg, 'Gene']})")
    fig.update_layout({'colorway': ['grey', 'red', 'red', 'blue', 'blue']})
    save_figure(fig, f"{save_path}/{cpg_id}_{cpg}")

add_layout(cmb, "", 'Methylation Level', '')
cmb.update_yaxes(autorange=False, range=[-0.0, 1.0])
cmb.update_layout({'colorway': ['grey', 'red', 'red'] * len(targets)})
cmb.update_layout(margin=dict(l=80, r=20, t=60, b=160))
cmb.update_xaxes(tickangle=90)
save_figure(cmb, f"{save_path}/combo")
