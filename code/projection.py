"""Full N=2 projected-energy estimator for the current grouped-spin H2 model.
Off-diagonal sum V must satisfy P2 V P2=V (verified for this Hamiltonian).
Uncertainty uses the delta method with shared Z numerator/denominator covariance.
"""
from math import sqrt

def parity(label,bits):
    return (-1)**sum(int(bits[i]) for i,p in enumerate(label) if p!='I')

def moment(values,counts):
    n=sum(counts.values())
    mean=sum(counts[b]*values[b] for b in counts)/n
    second=sum(counts[b]*values[b]**2 for b in counts)/n
    return mean,max(0.,second-mean*mean)/(n-1)

def projected_energy(coeff,groups):
    diag={p:c for p,c in coeff.items() if set(p)<=set('IZ')}
    off={p:c for p,c in coeff.items() if p not in diag}
    z=groups['ZZZZ'];total=sum(z.values());kept=sum(n for b,n in z.items() if b.count('1')==2)
    if kept<2:raise ValueError('Insufficient N=2 shots; cannot evaluate the projected objective')
    p2=kept/total
    hz={b:sum(c*parity(p,b) for p,c in diag.items()) for b in z}
    A,_=moment({b:hz[b] if b.count('1')==2 else 0. for b in z},z)
    raw_diag,raw_diag_var=moment(hz,z)
    off_mean=off_var=0.
    for label,c in off.items():
        counts=groups[label]
        mean,var=moment({b:c*parity(label,b) for b in counts},counts)
        off_mean+=mean;off_var+=var
    E=(A+off_mean)/p2
    influence={b:hz[b]-E if b.count('1')==2 else 0. for b in z}
    _,var=moment(influence,z)
    return {'energy':E,'se':sqrt((var+off_var)/p2**2),
            'raw_energy':raw_diag+off_mean,'raw_se':sqrt(raw_diag_var+off_var),
            'p2':p2,'z_shots_kept':kept,'z_shots_total':total,
            'diagonal_numerator':A,'offdiagonal_numerator':off_mean}
