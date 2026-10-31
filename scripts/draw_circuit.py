"""Draw the recorded build_circuit(theta, 'ZZZZ') at the logical gate level."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

root=Path(__file__).resolve().parent.parent
out=root/'results'
fig,ax=plt.subplots(figsize=(12.8,3.1),layout='constrained')
ax.set_xlim(-1.8,12.5);ax.set_ylim(-.65,3.7);ax.axis('off')
for q,physical in enumerate([6,7,12,17]):
    y=3-q
    ax.plot([0,11.8],[y,y],color='#334155',lw=1.4)
    ax.text(-.25,y,rf'$q_{q}$ (Q{physical})  $|0\rangle$',ha='right',va='center',fontsize=11)
def gate(x,q,label):
    y=3-q;ax.add_patch(Rectangle((x-.36,y-.26),.72,.52,facecolor='white',edgecolor='#1e40af',lw=1.3,zorder=3))
    ax.text(x,y,label,ha='center',va='center',fontsize=10,zorder=4)
def cz(x,a,b):
    ya,yb=3-a,3-b;ax.plot([x,x],[ya,yb],color='#111827',lw=1.5,zorder=3)
    ax.scatter([x,x],[ya,yb],s=35,color='#111827',zorder=4)
gate(1,0,r'$R_y(-4\theta)$');gate(1,1,'H')
cz(2.5,0,1);gate(4,1,'H');gate(4,2,'H')
cz(5.5,1,2);gate(7,2,'H');gate(7,3,'H')
cz(8.5,2,3);gate(10,0,'X');gate(10,2,'X');gate(10,3,'H')
for q in range(4):gate(11.4,q,'M')
ax.text(5.4,3.55,'Parameterized state preparation with ZZZZ measurement',ha='center',fontsize=12,fontweight='bold')
for ext in ('png','svg'):fig.savefig(out/f'parameterized_circuit_zzzz.{ext}',dpi=200)
plt.close(fig)
