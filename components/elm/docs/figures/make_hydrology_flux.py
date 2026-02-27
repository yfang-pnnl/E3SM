# -*- coding: utf-8 -*-
"""
Draw the ELM hydrology flux diagram (precipitation -> soil column recharge)
and save as hydrology_flux.png in the figures directory.
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

fig, ax = plt.subplots(figsize=(10, 16))
ax.set_xlim(0, 10)
ax.set_ylim(0, 16)
ax.axis('off')

C = {
    'precip': ('#cce5ff', '#3399ff'),
    'canopy': ('#d4edda', '#28a745'),
    'snow':   ('#e8f4fd', '#5bc0de'),
    'surf':   ('#fff3cd', '#e6ac00'),
    'soil':   ('#fdebd0', '#e67e22'),
    'charge': ('#e2d9f3', '#6f42c1'),
    'loss':   ('#f8d7da', '#dc3545'),
}

def box(ax, x, y, w, h, label, ck, fontsize=8.5):
    fc, ec = C[ck]
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.1",
                                fc=fc, ec=ec, lw=1.5, zorder=3))
    ax.text(x + w/2.0, y + h/2.0, label, ha='center', va='center',
            fontsize=fontsize, zorder=4, multialignment='center')

def arrow(ax, x1, y1, x2, y2, color='#555555'):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=1.4), zorder=2)

# boxes
box(ax, 3.5, 14.5, 3.0, 1.2,
    'Precipitation\nqflx_prec_grnd\n(rain + snow)', 'precip', 9)

box(ax, 0.3, 12.0, 4.8, 2.1,
    'Canopy  (H2OCAN)\n'
    'throughfall: qflx_through_rain\n'
    '   qflx_through_snow\n'
    'drip: qflx_drain_cansnow', 'canopy', 8.5)

box(ax, 5.7, 12.5, 3.8, 1.0,
    'Canopy evaporation\nqflx_evap_can', 'loss', 8.5)

box(ax, 0.3, 9.5, 4.8, 2.1,
    'Snow Pack  (H2OSNO)\n'
    'snowmelt: qflx_snow_melt\n'
    'sublimation: qflx_sub_snow', 'snow', 8.5)

box(ax, 0.3, 7.2, 4.8, 1.9,
    'Soil Surface  (fh2osfc)\n'
    'surface runoff: qflx_surf\n'
    'infiltration: qflx_infl', 'surf', 8.5)

box(ax, 0.3, 1.8, 4.8, 5.1,
    'Soil Column  (j = 1 ... nlevsoi)\n\n'
    '  Layer 1   H2OSOI,1\n'
    '  Layer 2   H2OSOI,2\n'
    '       ...\n'
    '  Layer n   H2OSOI,n\n', 'soil', 8.5)

box(ax, 5.7, 5.7, 3.8, 1.0,
    'Soil evaporation  qflx_evap_soi\n(top layer only)', 'loss', 8.5)

box(ax, 5.7, 3.6, 3.8, 1.7,
    'Transpiration sink\nqflx_rootsoi,j\n(per root-containing layer)', 'loss', 8.5)

box(ax, 3.5, 0.2, 3.0, 1.2,
    'Aquifer Recharge\nqflx_recharge', 'charge', 9)

# arrows
arrow(ax, 5.0, 14.5, 3.3, 14.1)
arrow(ax, 2.7, 12.0, 2.7, 11.6, '#28a745')
ax.annotate('', xy=(2.7, 9.1), xytext=(2.7, 12.0),
            arrowprops=dict(arrowstyle='->', color='#28a745', lw=1.2,
                            connectionstyle='arc3,rad=0.38'), zorder=2)
arrow(ax, 5.1, 13.0, 5.7, 13.0, '#dc3545')
arrow(ax, 2.7, 9.5, 2.7, 9.1, '#5bc0de')
arrow(ax, 2.7, 7.2, 2.7, 6.9, '#e6ac00')
arrow(ax, 5.1, 6.2, 5.7, 6.2, '#dc3545')
arrow(ax, 5.1, 4.5, 5.7, 4.5, '#dc3545')
arrow(ax, 2.7, 1.8, 5.0, 1.4, '#e67e22')

# section labels
for txt, yy in [('CANOPY', 13.05), ('SNOW', 10.55),
                ('SURFACE', 8.15), ('SOIL COLUMN', 4.35)]:
    ax.text(9.75, yy, txt, fontsize=8, color='#444444', rotation=90,
            va='center', ha='center', style='italic')

ax.set_title('ELM Hydrology: Precipitation to Soil Column Recharge',
             fontsize=12, fontweight='bold', pad=10)

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hydrology_flux.png')
plt.savefig(out, dpi=150, bbox_inches='tight')
print('Saved: ' + out)
