import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import AutoMinorLocator
from matplotlib.patches import Patch
from matplotlib.legend_handler import HandlerPatch
from matplotlib.transforms import Bbox

models = ['ebv_model2_dm15', 'ebv_model2_st', 'max_model_st']

model_styles = {
    'ebv_model2_dm15': {'color': 'royalblue',   'hatch': '///////////', 'label': r'$EBV\_model2$, $dm15$'},
    'ebv_model2_st':   {'color': 'darkorange',  'hatch': '///////////', 'label': r'$EBV\_model2$, $s_{BV}$'},
    'max_model_st':    {'color': 'forestgreen', 'hatch': '///////////', 'label': r'$max\_model$, $s_{BV}$'},
}

Tmax_2019np = 0
plt.rcParams['hatch.linewidth'] = 0.5

# ┌─────────────────────────────────────────────────────────────────┐
# │  INSET SIZE CONTROLS — adjust these two values only            │
# │  f_lc  : fraction of the combined inset height for light curve │
# │  f_res : fraction of the combined inset height for residuals   │
# │  They must sum to 1.                                           │
f_lc  = 0.65   # ← change this to redistribute height
f_res = 0.35  # ← change this (= 1 - f_lc)
# │                                                                 │
# │  Total inset footprint (in axes-fraction units):               │
inset_width  = '35%'   # ← width of both panels
inset_height = '50%'   # ← total combined height
inset_x      = 0.15    # ← left anchor in axes fraction
inset_y      = 0.02    # ← bottom anchor in axes fraction
gap_fig      = 0.004   # ← gap between the two panels (figure coords)
# └─────────────────────────────────────────────────────────────────┘


def parse_snoopy_file(filepath, target_filter):
    data = []
    in_target = False
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('filter'):
                in_target = (line.split()[1] == target_filter)
                continue
            if in_target:
                parts = line.split()
                data.append([float(x) for x in parts[:3]])
    return np.array(data)


fig, ax = plt.subplots(figsize=(8, 6))
fig.subplots_adjust(left=0.10, right=0.97, top=0.96, bottom=0.11)

# ── Observed data ─────────────────────────────────────────────────────────────
obs     = parse_snoopy_file('2019np.txt', 'Is')
t_obs   = obs[:, 0]
mag_obs = obs[:, 1]
err_obs = obs[:, 2]

ax.errorbar(
    t_obs - Tmax_2019np, mag_obs, yerr=err_obs,
    fmt='o', color='k', ecolor='k',
    elinewidth=1.0, capsize=2.5, markersize=4,
    label='SN2019np (observed)', zorder=10,
)

# ── Create a single dummy inset just to get the total Bbox ────────────────────
# We immediately compute its figure-coordinate Bbox, then remove it and
# manually place two stacked axes covering the same footprint.
_dummy = inset_axes(ax, width=inset_width, height=inset_height,
                    loc='lower left',
                    bbox_to_anchor=(inset_x, inset_y, 1, 1),
                    bbox_transform=ax.transAxes)
fig.canvas.draw()
bbox_total = _dummy.get_position()   # Bbox in figure coords
_dummy.remove()

# Split the total Bbox into light-curve (top) and residuals (bottom)
y0 = bbox_total.y0
y1 = bbox_total.y1
h  = y1 - y0

y_split = y0 + f_res * h   # boundary between the two panels

bbox_res = Bbox([[bbox_total.x0, y0],
                 [bbox_total.x1, y_split - gap_fig / 2]])

bbox_lc  = Bbox([[bbox_total.x0, y_split + gap_fig / 2],
                 [bbox_total.x1, y1]])

ax_inset = fig.add_axes(bbox_lc)
ax_res   = fig.add_axes(bbox_res)

t_zoom_lo,   t_zoom_hi   = 58500 - Tmax_2019np, 58517 - Tmax_2019np
mag_zoom_lo, mag_zoom_hi = 12.94, 13.8

ax_inset.errorbar(
    t_obs - Tmax_2019np, mag_obs, yerr=err_obs,
    fmt='o', color='k', ecolor='k',
    elinewidth=1.0, capsize=2.5, markersize=3,
    zorder=10,
)

# ── Model fits + collect residuals ────────────────────────────────────────────
splines_mid   = {}
spline_ranges = {}

for model in models:
    style  = model_styles[model]
    colour = style['color']
    hatch  = style['hatch']
    label  = style['label']

    filename = f'{model}/SN2019np_lc_Is_model.dat'
    data = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            parts = line.replace(',', '').split()
            data.append([float(x) for x in parts])

    data    = np.array(data)
    t       = data[:, 0]
    mag     = data[:, 1]
    mag_err = data[:, 2]

    cs_mid = CubicSpline(t - Tmax_2019np, mag)
    cs_lo  = CubicSpline(t - Tmax_2019np, mag - mag_err)
    cs_hi  = CubicSpline(t - Tmax_2019np, mag + mag_err)

    splines_mid[model]   = cs_mid
    spline_ranges[model] = (t.min() - Tmax_2019np, t.max() - Tmax_2019np)

    t_dense   = np.linspace(t.min() - Tmax_2019np, t.max() - Tmax_2019np, 500)
    mag_dense = cs_mid(t_dense)
    lo_dense  = cs_lo(t_dense)
    hi_dense  = cs_hi(t_dense)

    for a in [ax, ax_inset]:
        a.fill_between(
            t_dense, lo_dense, hi_dense,
            facecolor='none', edgecolor=colour,
            hatch=hatch, linewidth=0.0, alpha=0.65,
            label=f'{label} ($\\pm1\\sigma$)' if a is ax else None,
        )
        a.plot(t_dense, lo_dense, color=colour, lw=0.8, ls='--')
        a.plot(t_dense, hi_dense, color=colour, lw=0.8, ls='--')
        a.plot(t_dense, mag_dense, color=colour, lw=1.2)

    # Residuals within zoom window only
    t_shifted = t_obs - Tmax_2019np
    in_zoom   = (t_shifted >= t_zoom_lo) & (t_shifted <= t_zoom_hi) & \
                (t_shifted >= spline_ranges[model][0]) & \
                (t_shifted <= spline_ranges[model][1])

    t_res   = t_shifted[in_zoom]
    res     = mag_obs[in_zoom] - cs_mid(t_res)
    err_res = err_obs[in_zoom]

    ax_res.errorbar(
        t_res, res/err_res, yerr=err_res,
        fmt='o', color=colour, ecolor=colour,
        elinewidth=1.0, capsize=2.5, markersize=4,
        zorder=10,
    )

ax_res.axhline(0, color='k', lw=0.8, ls='--', zorder=1)

# ── Zoom inset formatting ─────────────────────────────────────────────────────
ax_inset.set_xlim(t_zoom_lo, t_zoom_hi)
ax_inset.set_ylim(mag_zoom_hi, mag_zoom_lo)   # inverted
ax_inset.set_facecolor('white')
ax_inset.minorticks_on()
ax_inset.tick_params(which='both', direction='in', labelsize=12)
ax_inset.tick_params(which='major', length=5, width=1.5)
ax_inset.tick_params(which='minor', length=2, width=1.0)
ax_inset.xaxis.set_minor_locator(AutoMinorLocator())
ax_inset.yaxis.set_minor_locator(AutoMinorLocator())
ax_inset.xaxis.set_major_locator(plt.MaxNLocator(3))
ax_inset.yaxis.set_major_locator(plt.MaxNLocator(4))
plt.setp(ax_inset.get_xticklabels(), visible=False)
#plt.setp(ax_inset.get_yticklabels(), visible=False)
for spine in ax_inset.spines.values():
    spine.set_linewidth(0.8)
    spine.set_edgecolor('0.3')

# ── Residuals inset formatting ────────────────────────────────────────────────
ax_res.set_xlim(t_zoom_lo, t_zoom_hi)
ax_res.set_facecolor('white')
ax_res.minorticks_on()
ax_res.tick_params(which='both', direction='in', labelsize=12)
ax_res.tick_params(which='major', length=5, width=1.5)
ax_res.tick_params(which='minor', length=2, width=1.0)
ax_res.xaxis.set_minor_locator(AutoMinorLocator())
ax_res.yaxis.set_minor_locator(AutoMinorLocator())
ax_res.xaxis.set_major_locator(plt.MaxNLocator(3))
ax_res.yaxis.set_major_locator(plt.MaxNLocator(3))
plt.setp(ax_res.get_xticklabels(), visible=False)
#ax_res.set_xlabel('MJD', fontsize=10)
ax_res.set_ylabel('Norm.Res.', fontsize=10)
for spine in ax_res.spines.values():
    spine.set_linewidth(0.8)
    spine.set_edgecolor('0.3')

# ── Main axes formatting ──────────────────────────────────────────────────────
ax.invert_yaxis()
ax.set_xlabel('Time (Modified Julian Date)', fontsize=18)
ax.set_ylabel('apparent magnitude', fontsize=18)
ax.minorticks_on()
ax.tick_params(which='both', direction='in', labelsize=17)
ax.tick_params(which='major', length=5, width=1.5)
ax.tick_params(which='minor', length=2, width=1.0)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.set_major_locator(plt.MaxNLocator(6, integer=True))
ax.axvspan(t_zoom_lo, t_zoom_hi, color='silver', alpha=0.5, lw=0)


ax.text(0.015, 0.93, 'Cousins $I$', transform=ax.transAxes,
        fontsize=16, va='bottom', ha='left', color='black', fontweight='bold')

# ── Legend ────────────────────────────────────────────────────────────────────
class ThickHatchHandler(HandlerPatch):
    def create_artists(self, legend, orig_handle, xdescent, ydescent,
                       width, height, fontsize, trans):
        p = super().create_artists(legend, orig_handle, xdescent, ydescent,
                                   width, height, fontsize, trans)[0]
        p.set_linewidth(2.5)
        return [p]

legend_handles = [
    Patch(facecolor='none', edgecolor=style['color'],
          hatch=style['hatch'], label=f"{style['label']} ($\\pm1\\sigma$)")
    for style in model_styles.values()
]
legend_handles.append(
    plt.Line2D([0], [0], marker='o', color='k', linestyle='None',
               markersize=4, label='2019np (observed)')
)

ax.legend(
    handles=legend_handles,
    handler_map={Patch: ThickHatchHandler()},
    fontsize=11,
    handlelength=3.4, handleheight=2,
    labelspacing=0.6, borderpad=0.4, frameon=False,
    loc='upper right',
    bbox_to_anchor=(0.995, 0.995),
    bbox_transform=ax.transAxes,
)

plt.savefig('../Figures/model_comparator_eek.png', dpi=150)
plt.show()