
import json
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.colors import LogNorm
from numpy import linspace, asarray

rc('text', usetex=True)

def read_data(data_path):
	with open(data_path, 'r') as f:
		data = f.read()
	data = json.loads(data)
	return asarray(data['d_bp']).T, asarray(data['d_su']).T


D2 = 2
D1 = 2*D2*D2 + 10

data_path = 'result/periodic_xxx_rdm2_m_10_D2%s_D1%s.json'%(D2, D1)


f_bp, f_su = read_data(data_path)

print(f_bp.shape)

fontsize = 20
labelsize = 16
linewidth = 3.
markersize = 10

colors = ['c', 'g', 'b', 'y', 'r', 'k']
markers = ['o', 's', '*', '+', '^', '*']


cdict = {
  'red'  :  ( (0.0, 0.25, .25), (0.02, .59, .59), (1., 1., 1.)),
  'green':  ( (0.0, 0.0, 0.0), (0.02, .45, .45), (1., .97, .97)),
  'blue' :  ( (0.0, 1.0, 1.0), (0.02, .75, .75), (1., 0.45, 0.45))
}

cm = matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

fig, ax = plt.subplots(1, 2, figsize=(12, 6))

D = 2

alphas = [1., 0.6, 0.2]

color_a = 'g'
color_b = 'c'

y = [i for i in range(1,11)]
x = [i for i in range(1,11)]

# print(f_bp)
print(f_bp.min(), f_bp.max(), f_bp.mean())
pcm = ax[0].pcolormesh(x, y, f_bp, cmap='PuBu_r', norm=LogNorm(vmin=0.0001, vmax=0.25))
ax[0].set_xlabel(r'Horizontal bonds', fontsize=fontsize)

ax[0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0].locator_params(nbins=6, integer=True)
ax[0].set_title(r'$\mathcal{D}(\rho_{bp}, \rho_{bmps})$ (D=%s)'%(D), fontsize=fontsize)

cbar = fig.colorbar(pcm, ax=ax[0])
cbar.ax.tick_params(labelsize=labelsize) 


print(f_su.min(), f_su.max(), f_su.mean())
pcm = ax[1].pcolormesh(x, y, f_su, cmap='PuBu_r', norm=LogNorm(vmin=0.0001, vmax=0.25))
ax[1].set_xlabel(r'Horizontal bonds', fontsize=fontsize)

ax[1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1].locator_params(nbins=6, integer=True)
ax[1].set_title(r'$\mathcal{D}(\rho_{su}, \rho_{bmps})$ (D=%s)'%(D), fontsize=fontsize)

cbar = fig.colorbar(pcm, ax=ax[1])
cbar.ax.tick_params(labelsize=labelsize) 



# im2 = ax[1].imshow(f_su, extent=(1, 21, 1, 20), aspect='auto')
# divider = make_axes_locatable(ax[1])
# cax = divider.append_axes("right", size="5%", pad=0.05)
# cb = plt.colorbar(im2, cax=cax)
# cb.ax.tick_params(labelsize=fontsize)
# cb.ax.locator_params(nbins=6)
# ax[1].set_xlabel(r'Horizontal bonds', fontsize=fontsize)

# ax[1].tick_params(axis='both', which='major', labelsize=labelsize)
# ax[1].locator_params(nbins=6, integer=True)
# ax[1].set_title(r'$\mathcal{D}(\rho_{su}, \rho_{bmps})$ (D=%s)'%(D), fontsize=fontsize)


# ax[0].tick_params(axis='both', which='major', labelsize=labelsize)
# ax[0].tick_params(axis='y', which='both')
# ax[0].locator_params(nbins=6)
# ax[0].set_ylabel(r'$E$', fontsize=fontsize, rotation=0)
# ax[0].set_xlabel(r'$iterations$', fontsize=fontsize)
# ax[0].set_xlim(0,1)
# # ax[0].set_ylim(0,1.02)
# # ax[0].yaxis.set_label_coords(-0.2, 0.55)
# # ax[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
# # ax[1, 1].set_ylim(0)
# # ax1.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
# # ax[1,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
# ax[0].annotate(r'(a)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
# ax[0].legend(loc = 'upper right', fontsize=12)
# ax[0].set_title("Heisenberg model (10*10)", fontsize=12)



plt.tight_layout(pad=0.5)

# plt.savefig('compare_rdm2_D%s.pdf'%(D), dpi=150)

plt.show()
