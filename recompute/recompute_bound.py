"""Recompute the original Fisher/gamut illustration with explicit conventions.

The input arrays are unchanged Colour 0.4.7 Stockman--Sharpe 2-degree cone
fundamentals at 1-nm intervals, 390--830 nm. The primary boundary is the convex
hull of normalized responses from 390--700 nm; open-edge angular derivatives
are evaluated analytically. An additional control uses the entire spectral
locus, interpolated piecewise linearly in Cartesian chromaticity and closed
by its endpoint chord. Rays from the specified E white intersect each polygon.
The Fisher metric is evaluated at the intersection itself, not interpolated
from endpoint metric values. PCHIP controls resolve the original locus's
narrow violet feature. No empirical observations are synthesized or fitted.

The sampled inequality is a convention-dependent consistency calculation:
the simplex metric is assigned scale c=1, without psychophysical calibration.
"""
from pathlib import Path
import hashlib
import json
import platform
import numpy as np
import scipy
from scipy.interpolate import PchipInterpolator
from scipy.spatial import ConvexHull
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
WL = np.load(ROOT / 'inputs/wavelengths.npy')
LMS = np.load(ROOT / 'inputs/LMS.npy')
INTEGRALS = np.trapezoid(LMS, WL, axis=0)


def cross(a, b):
    return a[..., 0]*b[..., 1] - a[..., 1]*b[..., 0]


def boundary(theta, balanced, hull=False):
    sensitivity = LMS / INTEGRALS if balanced else LMS
    catches = np.trapezoid(sensitivity, WL, axis=0)
    white = catches / catches.sum()
    select = (WL >= 390) & (WL <= 700)
    p = sensitivity[select] / sensitivity[select].sum(axis=1)[:, None]
    vertices = ConvexHull(p[:,:2]).vertices if hull else np.arange(len(p))
    p = p[vertices]
    a = p[:, :2] - white[:2]
    delta = np.roll(a, -1, axis=0) - a
    u = np.column_stack([np.cos(theta), np.sin(theta)])
    den = cross(u[:, None, :], delta[None, :, :])
    with np.errstate(divide='ignore', invalid='ignore'):
        radius = cross(a, delta)[None, :] / den
        weight = cross(a[None, :, :], u[:, None, :]) / den
    valid = (radius > 0) & (weight >= -2e-12) & (weight <= 1+2e-12)
    assert np.all(valid.sum(axis=1) >= 1)
    radii = np.where(valid, radius, np.inf)
    segment = np.argmin(radii, axis=1)
    ii = np.arange(len(theta))
    radius = radii[ii, segment]
    weight = weight[ii, segment]
    chroma = p[segment]*(1-weight[:, None]) + p[(segment+1)%len(p)]*weight[:, None]
    assert np.all(chroma >= -1e-12)
    assert np.max(np.abs(chroma.sum(axis=1)-1)) < 1e-12
    numerator = np.column_stack([-np.sin(theta), np.cos(theta), np.sin(theta)-np.cos(theta)])**2
    with np.errstate(divide='ignore', invalid='ignore'):
        beta = np.sqrt((numerator / chroma).sum(axis=1))
    # Tabulated zero cone response is a singular Fisher-model boundary, not
    # an empirical inability to see the corresponding light.
    beta[np.any(chroma <= 0, axis=1)] = np.inf
    uprime = np.column_stack([-np.sin(theta), np.cos(theta)])
    exact_lhs = np.abs(cross(uprime, delta[segment]) / cross(u, delta[segment]))
    is_spectral = np.abs(np.roll(vertices,-1)-vertices)[segment] == 1
    return radius, beta, is_spectral, chroma, exact_lhs, white


def stats(theta, lhs, beta, mask):
    idx = np.flatnonzero(mask)
    if not len(idx):
        return {'count': 0}
    critical = idx[np.argmax((lhs[idx]/beta[idx])**2)]
    tight = idx[np.argmin(beta[idx]-lhs[idx])]
    peak = idx[np.argmax(lhs[idx])]
    return {
        'count': int(len(idx)),
        'satisfied_c1': int(np.sum(lhs[idx] < beta[idx])),
        'flat_beta1_violations': int(np.sum(lhs[idx] >= 1)),
        'beta_min': float(np.min(beta[idx])), 'beta_max': float(np.max(beta[idx])),
        'lhs_max': float(lhs[peak]), 'lhs_max_degrees': float(np.degrees(theta[peak])),
        'min_margin': float(beta[tight]-lhs[tight]),
        'min_margin_degrees': float(np.degrees(theta[tight])),
        'min_margin_lhs': float(lhs[tight]), 'min_margin_beta': float(beta[tight]),
        'c_star': float((lhs[critical]/beta[critical])**2),
        'c_star_degrees': float(np.degrees(theta[critical])),
    }


def run(balanced, n, hull=False):
    theta = np.linspace(-np.pi, np.pi, n, endpoint=False)
    dtheta = 2*np.pi/n
    radius, beta, spectral, chroma, exact, white = boundary(theta, balanced, hull=hull)
    lhs = np.abs(np.roll(radius, -1)-np.roll(radius, 1))/(2*dtheta*radius)
    retained = np.isfinite(beta) & (beta < 10)
    # Guard against derivatives crossing the two chosen closure junctions.
    spectral_stencil = spectral & np.roll(spectral, 1) & np.roll(spectral, -1)
    chord_stencil = ~spectral & ~np.roll(spectral, 1) & ~np.roll(spectral, -1)
    result = {
        'n_angles': n, 'degrees_step': 360/n, 'white_lms': white.tolist(),
        'alpha_min': float(radius.min()), 'alpha_max': float(radius.max()),
        'singular_samples': int(np.sum(~np.isfinite(beta))),
        'retained_all': stats(theta, lhs, beta, retained),
        'retained_spectral': stats(theta, lhs, beta, retained & spectral_stencil),
        'retained_chord': stats(theta, lhs, beta, retained & chord_stencil),
        'spectral_samples': int(spectral.sum()), 'chord_samples': int((~spectral).sum()),
        'excluded_junction_stencils': int((~(spectral_stencil|chord_stencil)).sum()),
        'analytic_segment_derivative_spectral': stats(theta, exact, beta, retained & spectral_stencil),
        'analytic_segment_derivative_all': stats(theta, exact, beta, retained),
        'analytic_segment_derivative_mixture_edges': stats(theta, exact, beta, retained & ~spectral),
        'metric_scale_sweep_fixed_retained_spectral': [
            {'c': c, 'pass': int(np.sum(lhs[retained & spectral_stencil] < np.sqrt(c)*beta[retained & spectral_stencil]))}
            for c in [.01, .025, .05, .1, .25, .5, 1, 2]
        ],
    }
    data = np.column_stack([np.degrees(theta), radius, lhs, beta, spectral, retained, exact])
    np.savetxt(ROOT / f"{'hull_' if hull else ''}{'E_balanced' if balanced else 'raw_E_centered'}_{n}.csv", data, delimiter=',',
               header='theta_degrees,alpha,central_lhs,beta_c1,is_spectral,beta_less_10,exact_polygon_lhs', comments='')
    return result, (theta, radius, lhs, beta, spectral, retained)


def wavelength_check(balanced, interpolation):
    sensitivity = LMS / INTEGRALS if balanced else LMS
    catches = np.trapezoid(sensitivity, WL, axis=0)
    white = catches / catches.sum()
    if interpolation == 'pchip_LMS':
        lam = np.linspace(390,700,31001)
        f = PchipInterpolator(WL,sensitivity)
        c = f(lam)
        v = f.derivative()(lam)
        total = c.sum(axis=1)
        p = c / total[:,None]
        dp = (v*total[:,None]-c*v.sum(axis=1)[:,None])/total[:,None]**2
    else:
        p0=sensitivity[:311]/sensitivity[:311].sum(axis=1)[:,None]
        t=np.tile(np.linspace(.00001,.99999,101),310)
        ids=np.repeat(np.arange(310),101)
        p=p0[ids]*(1-t[:,None])+p0[ids+1]*t[:,None]
        dp=p0[ids+1]-p0[ids]
        lam=390+ids+t
    xy=p[:,:2]-white[:2]
    theta=np.arctan2(xy[:,1],xy[:,0])
    with np.errstate(divide='ignore',invalid='ignore'):
        lhs=np.abs((xy*dp[:,:2]).sum(axis=1)/cross(xy,dp[:,:2]))
        beta=np.sqrt(np.sin(theta)**2/p[:,0]+np.cos(theta)**2/p[:,1]+(np.sin(theta)-np.cos(theta))**2/p[:,2])
    beta[np.any(p<=0,axis=1)]=np.inf
    keep=(beta<10)&np.isfinite(lhs)
    result=stats(theta,lhs,beta,keep)
    idx=np.flatnonzero(keep)
    critical=idx[np.argmax((lhs[idx]/beta[idx])**2)]
    result.update({'interpolation':interpolation,'total_wavelength_samples':len(lam),
                   'c_star_lambda_nm':float(lam[critical]),'c_star_lhs':float(lhs[critical]),
                   'c_star_beta':float(beta[critical])})
    failed=keep&(lhs>=beta)
    result['failed_wavelength_range_nm']=[float(lam[failed].min()),float(lam[failed].max())] if failed.any() else None
    name=f"{'E_balanced' if balanced else 'raw_E_centered'}_{interpolation}"
    np.savez_compressed(ROOT/f'{name}.npz',wavelength_nm=lam,theta=theta,lhs=lhs,beta=beta,retained=keep)
    return result,(lam,lhs,beta)


def hull_dense_check(balanced):
    sensitivity=LMS/INTEGRALS if balanced else LMS
    catches=np.trapezoid(sensitivity,WL,axis=0)
    white=catches/catches.sum()
    p=sensitivity[:311]/sensitivity[:311].sum(axis=1)[:,None]
    ids=ConvexHull(p[:,:2]).vertices
    p=p[ids]
    t=np.tile(np.linspace(.00001,.99999,101),len(p))
    seg=np.repeat(np.arange(len(p)),101)
    q=p[seg]*(1-t[:,None])+p[(seg+1)%len(p)]*t[:,None]
    d=p[(seg+1)%len(p),:2]-p[seg,:2]
    xy=q[:,:2]-white[:2]
    theta=np.arctan2(xy[:,1],xy[:,0])
    lhs=np.abs((xy*d).sum(axis=1)/cross(xy,d))
    with np.errstate(divide='ignore',invalid='ignore'):
        beta=np.sqrt(np.sin(theta)**2/q[:,0]+np.cos(theta)**2/q[:,1]+(np.sin(theta)-np.cos(theta))**2/q[:,2])
    beta[np.any(q<=0,axis=1)]=np.inf
    retained=(beta<10)&np.isfinite(lhs)
    result=stats(theta,lhs,beta,retained)
    result.update({'n_vertices':len(p),'samples_per_open_edge':101,
                   'total_samples':len(t),'mixture_edges_nm':[
                       [int(390+ids[i]),int(390+ids[(i+1)%len(p)])]
                       for i in range(len(p)) if abs(ids[i]-ids[(i+1)%len(p)])!=1]})
    return result


def main():
    report = {
        'input_description': 'Stockman & Sharpe 2 Degree Cone Fundamentals; saved Colour Science 0.4.7 arrays',
        'input_domain_nm': [390, 830], 'locus_domain_nm': [390, 700],
        'raw_E_integrals': INTEGRALS.tolist(),
        'metric_scale': 'c=1 simplex Fisher convention; not psychophysically calibrated',
        'boundary_interpolation': 'Cartesian straight segments between 1-nm spectral chromaticities; exact endpoint chord',
        'derivative': 'periodic central difference at grid spacing; exact within-segment derivative reported separately',
        'retention_rule': 'beta(c=1)<10; fixed model singular-tail cutoff; not empirical threshold',
        'software': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy':scipy.__version__, 'matplotlib': matplotlib.__version__},
        'hashes': {name: hashlib.sha256((ROOT/'inputs'/name).read_bytes()).hexdigest() for name in ['LMS.npy','wavelengths.npy']},
        'normalizations': {},
        'wavelength_derivative_checks': {},
        'convex_hull_checks': {},
    }
    main_data = []
    wavelength_data = []
    hull_data = []
    for balanced, key in [(True, 'E_balanced'), (False, 'raw_E_centered')]:
        rows=[]
        for n in [360,720,1440,3600]:
            result, data=run(balanced,n)
            rows.append(result)
            if n == 360:
                main_data.append((key,data))
        report['normalizations'][key]=rows
        report['wavelength_derivative_checks'][key]=[]
        for interpolation in ['polygon_chromaticity','pchip_LMS']:
            result,data=wavelength_check(balanced,interpolation)
            report['wavelength_derivative_checks'][key].append(result)
            if interpolation=='pchip_LMS':
                wavelength_data.append((key,data))
        result,data=run(balanced,360,hull=True)
        result['dense_edge_check']=hull_dense_check(balanced)
        report['convex_hull_checks'][key]=result
        hull_data.append((key,data))
    (ROOT/'bound_results.json').write_text(json.dumps(report, indent=2)+'\n')
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':8, 'pdf.fonttype':42, 'svg.fonttype':'none'})
    fig, axes=plt.subplots(2,1,figsize=(6.8,4.0),sharex=True,constrained_layout=True)
    for ax,(key,data) in zip(axes,main_data):
        theta,radius,lhs,beta,spectral,retained=data
        deg=np.degrees(theta)
        ax.fill_between(deg,0,10,where=~spectral,color='#bd3e93',alpha=.10,step='mid')
        ax.plot(deg,lhs,color='#263b57',lw=1.2,label=r"$|\alpha'/\alpha|$")
        ax.plot(deg,np.where(beta<=10,beta,np.nan),color='#b52480',lw=1.3,label=r'$\beta$, simplex scale $c=1$')
        ax.axhline(1,color='#707070',lw=.65,ls=':',label=r'Fixed Euclidean benchmark $\beta=1$')
        ax.set_ylim(0,10)
        ax.set_ylabel('Angular coefficients')
        ax.text(.015,.90,'E-balanced channels' if key=='E_balanced' else 'Raw channels, integrated E center',transform=ax.transAxes)
        ax.grid(alpha=.15)
    axes[0].legend(loc='upper right',frameon=False,fontsize=7)
    axes[-1].set_xlabel(r'Hue angle $\theta$ (degrees)')
    axes[-1].set_xlim(-180,180)
    axes[-1].set_xticks(np.arange(-180,181,60))
    for ext in ['pdf','svg','png']:
        fig.savefig(ROOT/f'fisher_bound_comparison.{ext}',dpi=180)
    plt.close(fig)
    # Publication figure: declared primary convention and the narrow feature
    # that coarse angular sampling misses in the raw-channel robustness test.
    fig,axes=plt.subplots(2,1,figsize=(6.8,4.1),constrained_layout=True)
    theta,radius,lhs,beta,spectral,retained=main_data[0][1]
    deg=np.degrees(theta)
    ax=axes[0]
    ax.fill_between(deg,0,10,where=~spectral,color='#bd3e93',alpha=.10,step='mid')
    ax.plot(deg,lhs,color='#263b57',lw=1.2,label=r"$|\alpha'/\alpha|$")
    ax.plot(deg,np.where(beta<=10,beta,np.nan),color='#b52480',lw=1.3,label=r'$\beta$, scale $c=1$')
    ax.axhline(1,color='#707070',lw=.65,ls=':')
    ax.set(xlim=(-180,180),ylim=(0,10),ylabel='Angular coefficients',xlabel=r'Hue angle $\theta$ (degrees)')
    ax.text(.015,.88,'(a) E-balanced channels, 1° sampling',transform=ax.transAxes)
    ax.legend(loc='upper right',frameon=False,fontsize=7)
    ax.grid(alpha=.15)
    ax=axes[1]
    for (key,(lam,lhs,beta)),color in zip(wavelength_data,['#b52480','#263b57']):
        ax.plot(lam,(lhs/beta)**2,color=color,lw=1.3,label='E-balanced' if key=='E_balanced' else 'Raw channels, E-centered')
    ax.axhline(1,color='#707070',lw=.8,ls='--',label='Compatibility threshold at c=1')
    ax.set(xlim=(395,415),ylim=(0,1.5),ylabel=r"$(|\alpha'/\alpha|/\beta)^2$",xlabel='Wavelength (nm)')
    ax.text(.015,.90,'(b) Resolved violet shoulder, PCHIP derivatives',transform=ax.transAxes)
    ax.legend(loc='upper right',frameon=False,fontsize=7)
    ax.grid(alpha=.15)
    for ext in ['pdf','svg','png']:
        fig.savefig(ROOT/f'fisher_bound_sensitivity.{ext}',dpi=180)
    plt.close(fig)
    fig,axes=plt.subplots(2,1,figsize=(6.8,4.1),sharex=True,constrained_layout=True)
    theta,radius,lhs,beta,spectral,retained=hull_data[0][1]
    radius,beta,spectral,chroma,exact,white=boundary(theta,True,hull=True)
    deg=np.degrees(theta)
    ax=axes[0]
    ax.fill_between(deg,0,10,where=~spectral,color='#bd3e93',alpha=.10,step='mid')
    ax.plot(deg,exact,color='#263b57',lw=1.2,label=r"$|\alpha'/\alpha|$")
    ax.plot(deg,np.where(beta<=10,beta,np.nan),color='#b52480',lw=1.3,label=r'$\beta$, scale $c=1$')
    ax.axhline(1,color='#707070',lw=.65,ls=':')
    ax.set(ylim=(0,10),ylabel='Angular coefficients')
    ax.text(.015,.89,'(a) Convex hull, E-balanced channels',transform=ax.transAxes)
    ax.legend(loc='upper right',frameon=False,fontsize=7)
    ax.grid(alpha=.15)
    ax=axes[1]
    for balanced,color,label in [(True,'#b52480','E-balanced'),(False,'#263b57','Raw channels, E-centered')]:
        radius,beta,spectral,chroma,exact,white=boundary(theta,balanced,hull=True)
        required=(exact/beta)**2
        ax.plot(deg,np.where(beta<10,required,np.nan),color=color,lw=1.25,label=label)
    ax.axhline(1,color='#707070',lw=.8,ls='--',label='Chosen metric scale c=1')
    ax.set(xlim=(-180,180),ylim=(0,1.12),ylabel=r"Required scale $(|\alpha'/\alpha|/\beta)^2$",xlabel=r'Hue angle $\theta$ (degrees)')
    ax.text(.015,.82,'(b) Metric-scale sensitivity',transform=ax.transAxes)
    ax.legend(loc='upper right',frameon=False,fontsize=7)
    ax.set_xticks(np.arange(-180,181,60))
    ax.grid(alpha=.15)
    for ext in ['pdf','svg','png']:
        fig.savefig(ROOT/f'fisher_bound_convex_hull.{ext}',dpi=180)
    plt.close(fig)
    print(json.dumps({k:v[0] for k,v in report['normalizations'].items()},indent=2))


if __name__=='__main__':
    main()
