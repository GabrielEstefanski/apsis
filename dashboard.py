import sys, json
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from collections import defaultdict

plt.ion()
fig = plt.figure(figsize=(12, 7), tight_layout=True)
gs = gridspec.GridSpec(2, 3, figure=fig)

ax_trail  = fig.add_subplot(gs[:, 0])   # trajetórias XY (coluna inteira)
ax_energy = fig.add_subplot(gs[0, 1])   # E, K, U
ax_err    = fig.add_subplot(gs[1, 1])   # dE/E0
ax_vel    = fig.add_subplot(gs[0, 2])   # velocidades body 0
ax_info   = fig.add_subplot(gs[1, 2])   # texto: theta, dt, Lz

COLORS = ['#4C9BE8', '#2ECC71', '#E74C3C', '#F39C12', '#9B59B6']

steps = []
trails = defaultdict(lambda: {'x': [], 'y': []})

for line in sys.stdin:
    line = line.strip()
    if not line:
        continue
    try:
        d = json.loads(line)
    except json.JSONDecodeError:
        continue

    steps.append(d)
    xs = [s['step'] for s in steps]

    for ax in [ax_trail, ax_energy, ax_err, ax_vel, ax_info]:
        ax.cla()

    # trajetórias
    ax_trail.set_title('Trajetórias XY')
    ax_trail.set_xlabel('x'); ax_trail.set_ylabel('y')
    for bi, b in enumerate(d['bodies']):
        trails[bi]['x'].append(b['x'])
        trails[bi]['y'].append(b['y'])
        ax_trail.plot(trails[bi]['x'], trails[bi]['y'],
                      color=COLORS[bi % len(COLORS)], lw=1.2,
                      label=f"m={b['mass']}")
        ax_trail.scatter([b['x']], [b['y']], color=COLORS[bi % len(COLORS)], s=20 + b['mass'] * 8, zorder=5)
    ax_trail.legend(fontsize=8)

    # energia
    ax_energy.set_title('Energia')
    ax_energy.plot(xs, [s['E'] for s in steps], label='E', lw=1.5)
    ax_energy.plot(xs, [s['K'] for s in steps], label='K', lw=1.2, linestyle='--')
    ax_energy.plot(xs, [s['U'] for s in steps], label='U', lw=1.2, linestyle=':')
    ax_energy.legend(fontsize=8); ax_energy.set_xlabel('step')

    # erro relativo
    ax_err.set_title('|dE/E₀|')
    ax_err.semilogy(xs, [abs(s['dE']) if s['dE'] != 0 else 1e-20 for s in steps],
                    label='dE', lw=1.5)
    ax_err.semilogy(xs, [s['max_dE'] if s['max_dE'] != 0 else 1e-20 for s in steps],
                    label='max', lw=1.2, linestyle='--')
    ax_err.legend(fontsize=8); ax_err.set_xlabel('step')

    # velocidade body 0
    b0 = [s['bodies'][0] for s in steps]
    ax_vel.set_title('Velocidade — body 0')
    ax_vel.plot(xs, [b['vx'] for b in b0], label='vx', lw=1.5)
    ax_vel.plot(xs, [b['vy'] for b in b0], label='vy', lw=1.5, linestyle='--')
    ax_vel.legend(fontsize=8); ax_vel.set_xlabel('step')

    # painel de texto
    ax_info.axis('off')
    info = (
        f"step:    {d['step']}\n"
        f"theta:   {d['theta']:.4f}\n"
        f"dt:      {d['dt']:.4e}\n"
        f"Lz:      {d['Lz']:.6f}\n"
        f"max|dE|: {d['max_dE']:.3e}\n"
        f"E:       {d['E']:.6f}"
    )
    ax_info.text(0.1, 0.5, info, transform=ax_info.transAxes,
                 fontsize=11, verticalalignment='center', fontfamily='monospace')

    plt.pause(0.05)

plt.ioff()
plt.show()