
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

fig, ax = plt.subplots(2, 3, figsize=(12, 6), sharex=True, sharey=True)

D2 = 2
D1 = 2*D2*D2 + 10
B = 3.5

data_path = 'result/ising_rdm2_D2%s_D1%s_B%s.json'%(D2, D1, B)

f_bp, f_su = read_data(data_path)


alphas = [1., 0.6, 0.2]

color_a = 'g'
color_b = 'c'

y = [i for i in range(1,22)]
x = [i for i in range(1,21)]

# print(f_bp)
print(f_bp.min(), f_bp.max(), f_bp.mean())

cmap1 = 'PuBu_r'
cmap2 = cm

pcm = ax[0,0].pcolormesh(x, y, f_bp, cmap=cmap1, norm=LogNorm(vmin=1.0e-5, vmax=5e-3))
# ax[0,0].set_xlabel(r'Horizontal bonds', fontsize=fontsize)
ax[0,0].set_ylabel(r'$\mathcal{D}(\rho_{bp}, \rho_{bmps})$', fontsize=fontsize)


ax[0,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,0].locator_params(nbins=6, integer=True)
ax[0,0].set_title(r'D=%s'%(D2), fontsize=fontsize)

cbar = fig.colorbar(pcm, ax=ax[0,0])
cbar.ax.tick_params(labelsize=labelsize) 


print(f_su.min(), f_su.max(), f_su.mean())
pcm = ax[1,0].pcolormesh(x, y, f_su, cmap=cmap2, norm=LogNorm(vmin=5e-3, vmax=5e-1))
ax[1,0].set_xlabel(r'Horizontal bonds', fontsize=fontsize)
ax[1,0].set_ylabel(r'$\mathcal{D}(\rho_{su}, \rho_{bmps})$', fontsize=fontsize)

ax[1,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,0].locator_params(nbins=6, integer=True)
# ax[1,0].set_title(r'(D=%s)'%(D2), fontsize=fontsize)

cbar = fig.colorbar(pcm, ax=ax[1,0])
cbar.ax.tick_params(labelsize=labelsize) 


D2 = 3
D1 = 2*D2*D2 + 10
data_path = 'result/ising_rdm2_D2%s_D1%s_B%s.json'%(D2, D1, B)
f_bp, f_su = read_data(data_path)

print(f_bp.min(), f_bp.max(), f_bp.mean())
pcm = ax[0,1].pcolormesh(x, y, f_bp, cmap=cmap1, norm=LogNorm(vmin=1.0e-5, vmax=5e-3))
# ax[0,0].set_xlabel(r'Horizontal bonds', fontsize=fontsize)
# ax[0,0].set_ylabel(r'$\mathcal{D}(\rho_{bp}, \rho_{bmps})$', fontsize=fontsize)


ax[0,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,1].locator_params(nbins=6, integer=True)
ax[0,1].set_title(r'D=%s'%(D2), fontsize=fontsize)

cbar = fig.colorbar(pcm, ax=ax[0,1])
cbar.ax.tick_params(labelsize=labelsize) 


print(f_su.min(), f_su.max(), f_su.mean())
pcm = ax[1,1].pcolormesh(x, y, f_su, cmap=cmap2, norm=LogNorm(vmin=5e-3, vmax=5e-1))
ax[1,1].set_xlabel(r'Horizontal bonds', fontsize=fontsize)
# ax[1,1].set_ylabel(r'$\mathcal{D}(\rho_{su}, \rho_{bmps})$', fontsize=fontsize)

ax[1,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,1].locator_params(nbins=6, integer=True)
# ax[1,0].set_title(r'(D=%s)'%(D2), fontsize=fontsize)

cbar = fig.colorbar(pcm, ax=ax[1,1])
cbar.ax.tick_params(labelsize=labelsize) 


B = 3.

D2 = 4
D1 = 2*D2*D2 + 10
data_path = 'result/ising_rdm2_D2%s_D1%s_B%s.json'%(D2, D1, B)
f_bp, f_su = read_data(data_path)

print(f_bp.min(), f_bp.max(), f_bp.mean())
pcm = ax[0,2].pcolormesh(x, y, f_bp, cmap=cmap1, norm=LogNorm(vmin=1.0e-5, vmax=5e-3))
# ax[0,0].set_xlabel(r'Horizontal bonds', fontsize=fontsize)
# ax[0,0].set_ylabel(r'$\mathcal{D}(\rho_{bp}, \rho_{bmps})$', fontsize=fontsize)


ax[0,2].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,2].locator_params(nbins=6, integer=True)
ax[0,2].set_title(r'D=%s'%(D2), fontsize=fontsize)

cbar = fig.colorbar(pcm, ax=ax[0,2])
cbar.ax.tick_params(labelsize=labelsize) 


print(f_su.min(), f_su.max(), f_su.mean())
pcm = ax[1,2].pcolormesh(x, y, f_su, cmap=cmap2, norm=LogNorm(vmin=5e-3, vmax=5e-1))
ax[1,2].set_xlabel(r'Horizontal bonds', fontsize=fontsize)
# ax[1,1].set_ylabel(r'$\mathcal{D}(\rho_{su}, \rho_{bmps})$', fontsize=fontsize)

ax[1,2].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,2].locator_params(nbins=6, integer=True)
# ax[1,0].set_title(r'(D=%s)'%(D2), fontsize=fontsize)

cbar = fig.colorbar(pcm, ax=ax[1,2])
cbar.ax.tick_params(labelsize=labelsize) 


plt.tight_layout(pad=0.5)

plt.savefig('compare_rdm_ising.pdf', dpi=150)

plt.show()
