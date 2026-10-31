"""H2 VQE, R=0.700 Angstrom, STO-3G, grouped-spin Jordan-Wigner.
Run this file to submit hardware circuits. No REM or post-selection.
"""
from math import pi
from pathlib import Path
from datetime import datetime
import csv
import json
import matplotlib.pyplot as plt
from spinqit import Circuit, X, H, Ry, Rx, Rz, CZ, get_compiler
from spinqit.backend import get_spinq_cloud, SpinQCloudConfig

THETA0 = 0.0
SHOTS = 1024
STEPS = 20
LEARNING_RATE = 0.05
PHYSICAL_LAYOUT = [19, 18, 13, 8]  # q0,q1,q2,q3，零起始编号
BASES = ['ZZZZ', 'XXXX', 'XXYY', 'YYXX', 'YYYY']

# 哈密顿量含核间排斥能，Pauli 字符串从左至右对应 q0,q1,q2,q3。
COEFF = {
    'IIII': -0.04207898609029992,
    'ZIII': 0.1777128749502261, 'IIZI': 0.17771287495022617,
    'IZII': -0.24274280022634193, 'IIIZ': -0.24274280022634187,
    'ZZII': 0.12293305042565142, 'IIZZ': 0.12293305042565142,
    'ZIZI': 0.1705973836281266, 'IZIZ': 0.17627640720169607,
    'ZIIZ': 0.1676831942978458, 'IZZI': 0.1676831942978458,
    'XXXX': 0.04475014387219442, 'XXYY': 0.04475014387219442,
    'YYXX': 0.04475014387219442, 'YYYY': 0.04475014387219442,
}


def build_circuit(theta, basis):
    circ = Circuit()
    q = circ.allocateQubits(4)
    # 固定门结构：theta=0 时也保留全部 3 个 CZ。
    circ << (Ry, q[0], -4 * theta)
    circ << (H, q[1])
    circ << (CZ, (q[0], q[1]))
    circ << (H, q[1])
    circ << (H, q[2])
    circ << (CZ, (q[1], q[2]))
    circ << (H, q[2])
    circ << (H, q[3])
    circ << (CZ, (q[2], q[3]))
    circ << (X, q[0])
    circ << (X, q[2])

    # q3 末尾 H 与测量变换合并：Z -> H，X -> 无门，Y -> Rx(-pi/2)。
    q3_gates = {'Z': [(H,)], 'X': [], 'Y': [(Rx, -pi / 2)]}
    for gate, *angle in q3_gates[basis[3]]:
        circ << (gate, q[3], *angle)

    # q1、q2 的 H 后仍有作用于它们的 CZ，不能直接与测量 H 抵消。
    rotations = {'Z': [], 'X': [(H,)], 'Y': [(Rz, -pi / 2), (H,)]}
    for wire in range(3):
        for gate, *angle in rotations[basis[wire]]:
            circ << (gate, q[wire], *angle)
    return circ


def expectation(pauli, counts):
    total = sum(counts.values())
    return sum(n * (-1) ** sum(int(bits[i]) for i, p in enumerate(pauli) if p != 'I')
               for bits, n in counts.items()) / total


def energy(theta, label):
    value = COEFF['IIII']
    for basis in BASES:
        ir = get_compiler('native').compile(build_circuit(theta, basis), 0)
        config = SpinQCloudConfig()
        config.configure_platform('sqc_25_vp')
        config.configure_shots(SHOTS)
        config.configure_keep_layout(PHYSICAL_LAYOUT)
        config.configure_task(f'h2_3cz_{label}_{basis}', f'theta={theta}')
        result = cloud.execute(ir, config)
        counts = {str(bits).replace(' ', ''): int(n) for bits, n in result.counts.items()}
        # 保存原始计数和任务编号，便于复查。
        with (output / 'counts.jsonl').open('a', encoding='utf-8') as f:
            f.write(json.dumps({'label': label, 'theta': theta, 'basis': basis,
                               'task_code': result.task_code, 'counts': counts}) + '\n')
        for pauli, coefficient in COEFF.items():
            if pauli != 'IIII' and all(p == 'I' or p == b for p, b in zip(pauli, basis)):
                value += coefficient * expectation(pauli, counts)
    return value


if __name__ == '__main__':
    output = Path(__file__).resolve().parent / f'H2_VQE_results_{datetime.now():%Y%m%d_%H%M%S}'
    output.mkdir()
    cloud = get_spinq_cloud('<CLOUD_ACCOUNT>', '<LOCAL_PATH_REDACTED>')
    theta = THETA0
    with (output / 'convergence.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['step', 'theta_radian', 'energy_hartree'])
        e = energy(theta, 'initial')
        rows = [(0, theta, e)]
        writer.writerow(rows[-1]); f.flush()
        print(f'step=0 theta={theta:.8f} energy={e:.8f} Ha', flush=True)
        for step in range(1, STEPS + 1):
            # 原生 Ry 角度 beta=-4theta：链式法则后的参数移位梯度。
            e_plus = energy(theta + pi / 8, f'{step}_plus')
            e_minus = energy(theta - pi / 8, f'{step}_minus')
            gradient = 2 * (e_plus - e_minus)
            theta -= LEARNING_RATE * gradient
            e = energy(theta, f'{step}_updated')
            rows.append((step, theta, e))
            writer.writerow(rows[-1]); f.flush()
            print(f'step={step} theta={theta:.8f} energy={e:.8f} Ha', flush=True)
    plt.plot([r[0] for r in rows], [r[2] for r in rows], 'o-', label='Hardware VQE')
    plt.axhline(-1.1361894542189095, linestyle='--', color='red', label='Exact STO-3G')
    plt.xlabel('VQE step'); plt.ylabel('Energy (Ha)')
    plt.legend(); plt.grid(alpha=0.3); plt.tight_layout()
    plt.savefig(output / 'convergence.png', dpi=200)
    print('Results:', output)
