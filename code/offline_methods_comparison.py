"""Offline comparison; all uncertainties are finite-shot only. No cloud access."""
import ast, csv, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize

ROOT=Path(__file__).resolve().parent.parent
P=ROOT/'data'
OUT=ROOT/'results'
state=json.loads((P/'experiment_state.json').read_text(encoding='utf-8'))
tree=ast.parse((P/'circuit_source.py').read_text(encoding='utf-8'))
C=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='COEFF' for t in n.targets))
bits=[f'{i:04b}' for i in range(16)]; bases=['ZZZZ','XXXX','XXYY','YYXX','YYYY']; B=200; rng=np.random.default_rng(20261031)
def vec(c):return np.array([c.get(b,0) for b in bits],dtype=float)
def parity(p,b):return (-1)**sum(int(b[i]) for i,x in enumerate(p) if x!='I')
projector=np.array([int(b.count('1')==2) for b in bits],float)
diag=np.array([sum(c*parity(p,b) for p,c in C.items() if set(p)<=set('IZ')) for b in bits])
off=np.array([[C[lab]*parity(lab,b) for b in bits] for lab in bases[1:]])
def matrix(c0,c1):
    out=np.ones((1,1))
    for i in range(4):
        z=sum(c0[j] for j,b in enumerate(bits) if b[i]=='0')/sum(c0)
        o=sum(c1[j] for j,b in enumerate(bits) if b[i]=='0')/sum(c1)
        out=np.kron(out,np.array([[z,o],[1-z,1-o]]))
    return out
tasks=state['tasks']; calibrations={}
for label,prefix in [('start','ControlStart'),('end','ControlEnd200_')]:
    c0=vec(tasks[prefix+'0000_ZZZZ']['counts']);c1=vec(tasks[prefix+'1111_ZZZZ']['counts'])
    matrices=[matrix(c0,c1)]
    for _ in range(B):matrices.append(matrix(rng.multinomial(int(c0.sum()),c0/c0.sum()),rng.multinomial(int(c1.sum()),c1/c1.sum())))
    mats=np.array(matrices);inv=np.linalg.inv(mats)
    # Fixed Lipschitz bound for projected least squares, including bootstrap samples.
    L=np.linalg.norm(mats,ord=2,axis=(1,2))**2
    calibrations[label]=(mats,inv,L)

def simplexes(x):
    z=np.sort(x,axis=-1)[...,::-1];cs=np.cumsum(z,axis=-1)-1
    ks=np.arange(1,17);rho=np.sum(z-cs/ks>0,axis=-1)-1
    tau=np.take_along_axis(cs,rho[...,None],axis=-1)[...,0]/(rho+1)
    return np.maximum(x-tau[...,None],0)
def forward(M,q):return np.einsum('nij,nbj->nbi',M,q,optimize=True)
def transpose(M,x):return np.einsum('nji,nbj->nbi',M,x,optimize=True)
def estimates(q):
    p=q[:,0,:]@projector
    a=q[:,0,:]@(projector*diag)+np.einsum('nbi,bi->n',q[:,1:,:],off)
    return a,p,np.divide(a,p,out=np.full_like(a,np.nan),where=p>1e-12)
def methods(M,inv,L,f):
    qinv=forward(inv,f)
    # Richardson-Lucy / iterative Bayesian unfolding, fixed uniform prior.
    q=np.full_like(f,1/16);ib={}
    for i in range(1,101):
        q*=transpose(M,f/np.maximum(forward(M,q),1e-15))
        q/=q.sum(axis=-1,keepdims=True)
        if i in (10,50,100):ib[i]=q.copy()
    # FISTA: constrained least squares on probability simplex for each basis.
    q=np.full_like(f,1/16); y=q.copy();t=1.
    for i in range(120):
        g=transpose(M,forward(M,y)-f)
        new=simplexes(y-g/L[:,None,None]);tn=(1+np.sqrt(1+4*t*t))/2
        y=new+(t-1)/tn*(new-q);q=new;t=tn
    return qinv,ib,q

rows=[]; max_constraint=0.; min_p=1.; invalid=0
for step in range(len(state['updates'])+1):
    e=state['evaluations'][str(step)];counts=np.array([vec(tasks[f'{step}_{b}']['counts']) for b in bases]); f0=counts/1024
    # Same resampled measurement data are paired for both calibration models.
    fs=np.concatenate([f0[None,:,:],np.stack([rng.multinomial(1024,c/1024,size=B) for c in counts],axis=1)/1024],axis=0)
    row={'step':step,'theta':e['theta'],'raw':e['raw_energy'],'raw_se':e['raw_se'],'n2':e['energy'],'n2_se':e['se'],'retention':e['p2'],'ideal':e['ideal_energy']}
    for label in ('start','end'):
        M,inv,L=calibrations[label];qinv,ib,qls=methods(M,inv,L,fs)
        max_constraint=max(max_constraint,float(np.max(np.abs(qls.sum(axis=-1)-1))),float(max(0,-qls.min())))
        for name,q in [('inverse',qinv),('ibu10',ib[10]),('ibu50',ib[50]),('ibu100',ib[100]),('constrained_ls',qls)]:
            A,p,E=estimates(q);good=np.isfinite(E[1:]);invalid+=int((~good).sum());min_p=min(min_p,float(p[0]))
            row[f'{label}_{name}']=float(E[0]);row[f'{label}_{name}_se']=float(np.std(E[1:][good],ddof=1)) if good.sum()>1 else None
            row[f'{label}_{name}_numerator']=float(A[0]);row[f'{label}_{name}_denominator']=float(p[0]);row[f'{label}_{name}_invalid_bootstrap']=int((~good).sum())
    rows.append(row)
with (OUT/'offline_methods_comparison.csv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)

# Compare FISTA against an independent SLSQP constrained solve for 10 representative cases.
checks=[]
for label in ('start','end'):
 M,inv,L=calibrations[label]
 for step in (0,50,100,150,200):
    f=np.array([vec(tasks[f'{step}_{b}']['counts'])/1024 for b in bases]);q=methods(M[:1],inv[:1],L[:1],f[None,:,:])[2][0]
    for j in (0,2):
      m=M[0];v=f[j];obj=lambda x:float(np.sum((m@x-v)**2)/2);jac=lambda x:m.T@(m@x-v)
      sol=minimize(obj,np.full(16,1/16),jac=jac,method='SLSQP',bounds=[(0,1)]*16,constraints={'type':'eq','fun':lambda x:x.sum()-1,'jac':lambda x:np.ones(16)},options={'ftol':1e-13,'maxiter':500})
      if not sol.success:raise RuntimeError('SLSQP validation failed: '+sol.message)
      checks.append(abs(obj(q[j])-sol.fun))
assert max(checks)<1e-9 and max_constraint<1e-10

x=np.array([r['step'] for r in rows]); fig,axes=plt.subplots(2,1,figsize=(13,10),sharex=True,layout='constrained')
for ax,label in zip(axes,['start','end']):
    ax.errorbar(x,[r['raw'] for r in rows],yerr=[r['raw_se'] for r in rows],fmt='-',color='gray',alpha=.5,lw=1,elinewidth=.6,errorevery=10,label='Raw (±1 shot SE)')
    for key,title,color in [('n2','Full N=2 (no REM)','C1'),(f'{label}_inverse','Inverse REM + N=2','C2'),(f'{label}_ibu50','IBU (50 iters) + N=2','C3'),(f'{label}_constrained_ls','Constrained LS + N=2','C4')]:
        ax.errorbar(x,[r[key] for r in rows],yerr=[r[key+'_se'] for r in rows],fmt='-',color=color,lw=1,alpha=.8,elinewidth=.6,capsize=1,errorevery=10,label=title+' (±1 SE)')
    ax.plot(x,[r['ideal'] for r in rows],'k--',lw=1.3,label='Ideal simulator');ax.axhline(-1.1361894542189095,color='k',ls=':',alpha=.55,label='Exact ground');ax.set(title=f'{label.capitalize()} readout calibration',ylabel='Energy (Ha)');ax.grid(alpha=.2);ax.legend(ncol=2,fontsize=8)
axes[-1].set_xlabel('Adam update step');fig.savefig(OUT/'offline_methods_comparison.png',dpi=180);fig.savefig(OUT/'offline_methods_comparison.pdf');plt.close(fig)
report={'bootstrap_replicates':B,'bootstrap_seed':20261031,'bootstrap_invalid':invalid,'min_projected_denominator':min_p,'max_probability_constraint_error':max_constraint,'max_ls_objective_difference_slsqp':max(checks),'slsqp_checks':len(checks),'details':'Endpoint calibrations are applied globally for sensitivity; bootstrap resamples joint calibration counts and each data basis. Nonparametric resampling cannot represent unobserved calibration errors. IBU uniform prior, 10/50/100 iterations. LS minimizes ||M p - f||^2/2 over the simplex; 120 FISTA iterations. Each basis is treated independently; no common density-matrix constraint. Bootstrap SE is not systematic uncertainty.','step0':rows[0],'step200':rows[-1]}
(OUT/'offline_methods_validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
