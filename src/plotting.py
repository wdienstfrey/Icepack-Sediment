"""Figure 13 gate curves and common-domain maps; pending OOI-Image style copy."""
from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt
from .terminus import chronology, front, advance_rate


def style():
    # Provisional only. Replace with the exact OOI-Image style once located.
    plt.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False,
                         "savefig.dpi":300,"pdf.fonttype":42,"ps.fonttype":42})


def terminus_panel(output="results/figures",length=30000,width=3000):
    style()
    years,positions=chronology()
    fig,axes=plt.subplots(1,2,figsize=(7.2,2.6),layout="constrained")
    t=np.linspace(1953,2021,400)
    axes[0].plot(t,front(t),color="black",label="Smooth prescribed front")
    axes[0].plot(t,front(t,smooth=False),"--",color="0.5",label="Piecewise linear")
    axes[0].scatter(years,positions,color="black")
    axes[0].set(xlabel="Year",ylabel="Advance from 1953 (m)")
    axes[0].legend(fontsize=7)
    axes[1].fill_between([length/1000-4,(length+positions[-1])/1000],0,width/1000,color="0.95")
    for year,offset,color in zip(years,positions,plt.cm.viridis(np.linspace(0,1,4))):
        axes[1].plot([(length+offset)/1000]*2,[0,width/1000],color=color,label=str(int(year)))
    axes[1].axvline((length-1000)/1000,color="0.3",ls=":",label="Diagnostic gate")
    axes[1].set(xlabel="Downstream local x (km)",ylabel="Across-flow y (km)",
                title="Prescribed corridor fronts (schematic)")
    axes[1].legend(fontsize=7,loc="upper left",bbox_to_anchor=(1,1))
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    fig.savefig(out/"prescribed_terminus.pdf")
    fig.savefig(out/"prescribed_terminus.png")
    return fig


def comparison(fixed,advancing,output,smoke=False):
    style()
    paths=[Path(fixed),Path(advancing)]
    metadata=[json.loads((p/"metadata.json").read_text()) for p in paths]
    if metadata[0].get("input_kind")=="synthetic" and not smoke:
        raise ValueError("Synthetic inputs require --smoke labeling")
    if metadata[0]["input_sha256"]!=metadata[1]["input_sha256"]:
        raise ValueError("Experiments used different inputs")
    f,a=[np.genfromtxt(p/"history.csv",delimiter=",",names=True) for p in paths]
    spinup=metadata[0]["config"]["spinup_years"]
    mask=f["elapsed"]>=spinup-1e-8
    if not np.allclose(f["elapsed"],a["elapsed"]):
        raise ValueError("Output time grids differ")
    if not smoke and a["year"][-1]<2021-1e-6:
        raise ValueError("Historical run incomplete; use --smoke only for labeled synthetic outputs")
    f,a=f[mask],a[mask]
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    fig,axes=plt.subplots(3,1,sharex=True,figsize=(5.2,6),layout="constrained")
    t=a["year"]
    axes[0].plot(t,advance_rate(t,metadata[1]["config"]["smooth_front"]),color="black")
    axes[0].set_ylabel("Advance rate (m/year)")
    for ax,key,label in zip(axes[1:],["gate_flux_Gtyr","gate_h_m"],
                           ["Flux change (Gt/year)","Thickness change (m)"]):
        fixed_change=f[key]-f[key][0];advance_change=a[key]-a[key][0]
        ax.plot(t,fixed_change,color="0.5",label="Fixed-front change")
        ax.plot(t,advance_change,color="#0072B2",label="Advancing-front change")
        ax.plot(t,advance_change-fixed_change,color="black",label="Advance minus control")
        ax.set_ylabel(label)
    for ax,letter in zip(axes,"abc"):
        ax.text(0.02,0.92,letter,transform=ax.transAxes,fontweight="bold")
        for year in chronology()[0]: ax.axvline(year,color="0.7",ls="--",lw=.6)
    axes[-1].set_xlabel("Year")
    axes[1].legend(fontsize=7)
    fig.suptitle("Synthetic solver smoke test" if smoke else "Icepack prescribed-front experiment")
    fig.savefig(out/"figure13_comparison.pdf");fig.savefig(out/"figure13_comparison.png")
    snapshots=[[np.load(p/name) for name in ["initial.npz","final.npz"]] for p in paths]
    points=snapshots[0][0]["points"]/1000
    for initial,final in snapshots:
        np.testing.assert_allclose(initial["points"],final["points"])
    np.testing.assert_allclose(snapshots[0][0]["points"],snapshots[1][0]["points"])
    dh=[final["thickness"]-initial["thickness"] for initial,final in snapshots]
    du=[np.linalg.norm(final["velocity"],axis=1)-np.linalg.norm(initial["velocity"],axis=1)
        for initial,final in snapshots]
    fig,axes=plt.subplots(2,3,figsize=(8,4.8),layout="constrained")
    for row,(values,label) in enumerate([(dh,"Thickness change (m)"),(du,"Speed change (m/year)")]):
        arrays=[values[0],values[1],values[1]-values[0]]
        bound=max(1e-6,max(np.max(np.abs(v)) for v in arrays))
        for ax,values,title in zip(axes[row],arrays,["Fixed front","Advancing front","Advance minus control"]):
            im=ax.tricontourf(points[:,0],points[:,1],values,levels=np.linspace(-bound,bound,21),cmap="RdBu_r",extend="both")
            ax.set(title=title,xlabel="Downstream x (km)",ylabel="Across-flow y (km)")
        fig.colorbar(im,ax=list(axes[row]),label=label)
    fig.suptitle("Synthetic smoke: shared initial domain" if smoke else "Dynamic change on the shared 1953 domain")
    fig.savefig(out/"spatial_comparison.pdf");fig.savefig(out/"spatial_comparison.png")
    return fig
