"""Build the 201-step eight-method tables from the archived comparison data."""
from pathlib import Path
import csv, json, math
ROOT=Path(__file__).resolve().parent.parent
R=ROOT/'results'
with (R/'offline_methods_comparison.csv').open(encoding='utf-8-sig',newline='') as f: rows=list(csv.DictReader(f))
state=json.loads((ROOT/'data/experiment_state.json').read_text(encoding='utf8'))
methods=[('raw','Raw'),('n2','完整 N=2'),('start_inverse','直接求逆 REM + N=2，开始校准'),('start_ibu50','IBU + N=2，开始校准'),('start_constrained_ls','约束最小二乘 REM + N=2，开始校准'),('end_inverse','直接求逆 REM + N=2，结束校准'),('end_ibu50','IBU + N=2，结束校准'),('end_constrained_ls','约束最小二乘 REM + N=2，结束校准')]
assert len(rows)==201 and [int(x['step']) for x in rows]==list(range(201))
assert len(state['updates'])==200 and len(state['tasks'])==3017
errors=[]
for row in rows:
 step=int(row['step']);e=state['evaluations'][str(step)]
 for key,original in [('theta','theta'),('raw','raw_energy'),('raw_se','raw_se'),('n2','energy'),('n2_se','se'),('ideal','ideal_energy')]:
  if abs(float(row[key])-e[original])>1e-11:errors.append(f'step {step}: {key}')
 for basis in ['ZZZZ','XXXX','XXYY','YYXX','YYYY']:
  t=state['tasks'].get(f'{step}_{basis}')
  if not t or t['status']!='completed' or t['layout']!=[6,7,12,17] or t['shots']!=1024 or sum(t['counts'].values())!=1024 or abs(t['theta']-float(row['theta']))>1e-11:errors.append(f'step {step}: {basis}')
 for key,_ in methods:
  if not math.isfinite(float(row[key])) or not math.isfinite(float(row[key+'_se'])):errors.append(f'step {step}: {key} nonfinite')
if errors:raise ValueError('\n'.join(errors))
fields=['step','theta','ideal','retention']+[f for key,_ in methods for f in [key,key+'_se']]
with (R/'steps_000_200_eight_methods.csv').open('w',encoding='utf-8-sig',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows([{k:r[k] for k in fields} for r in rows])
with (R/'steps_000_200_long.csv').open('w',encoding='utf-8-sig',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['step','theta','method','energy_ha','sampling_se_ha','ideal_ha']);w.writeheader()
 for r in rows:
  for key,label in methods:w.writerow({'step':r['step'],'theta':r['theta'],'method':label,'energy_ha':r[key],'sampling_se_ha':r[key+'_se'],'ideal_ha':r['ideal']})
summary=[]
for key,label in methods:
 best=min(rows,key=lambda r:float(r[key]));summary.append({'method':label,'method_key':key,'initial':float(rows[0][key]),'initial_se':float(rows[0][key+'_se']),'final':float(rows[-1][key]),'final_se':float(rows[-1][key+'_se']),'minimum':float(best[key]),'minimum_se':float(best[key+'_se']),'minimum_step':int(best['step']),'theta_at_minimum':float(best['theta']),'ideal_at_minimum':float(best['ideal'])})
with (R/'method_summary.csv').open('w',encoding='utf-8-sig',newline='') as f:
 w=csv.DictWriter(f,fieldnames=summary[0]);w.writeheader();w.writerows(summary)
out=['# H₂ 0–200 步：八种方法结果','','单位 Ha；每项按“能量 ± 采样标准误”显示。理想值为当前 θ 的模拟结果。开始和结束校准各自作用于整条曲线，只用于敏感性分析。','', '| 步数 | θ (rad) | 理想模拟器 | '+' | '.join(label for _,label in methods)+' |', '|---:|---:|---:|'+ '|'.join(['---:']*8)+'|']
for r in rows:out.append('| '+r['step']+' | '+f"{float(r['theta']):.8f}"+' | '+f"{float(r['ideal']):.6f}"+' | '+' | '.join(f"{float(r[key]):.6f} ± {float(r[key+'_se']):.6f}" for key,_ in methods)+' |')
(R/'0-200逐步结果.md').write_text('\n'.join(out)+'\n',encoding='utf8')
(R/'integrity_check.json').write_text(json.dumps({'steps':len(rows),'methods':len(methods),'long_rows':len(rows)*len(methods),'issues':errors,'check':'theta, raw, N2 and ideal compared to saved state; five main-point counts, layout, status, shots validated; finite method results and SE','caveat':'does not establish statistical or physical accuracy of the mitigation'},ensure_ascii=False,indent=2),encoding='utf8')
print('Validated 201 steps, 8 methods and 1608 long-format rows; no inconsistencies.')
