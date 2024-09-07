import json
from matplotlib import pyplot as plt
from matplotlib import rc
from numpy import linspace

rc('text', usetex=True)

def read_data(data_path):
	with open(data_path, 'r') as f:
		data = f.read()
		data = json.loads(data)
	return data['energies'], data['final_energy']

def gen_bp_result_path(m, n, D1, D2, block_size):
	return "result/Table2_bp_m%s_n%s_blocksize%s_%s_D1%s_D2%s.json"%(m, n, block_size[0], block_size[1], D1, D2)



m = 8
n = 8
block_size = (4, 4)

fig_8_paths = [gen_bp_result_path(m, n, 2*D2*D2+10, D2, block_size) for D2 in range(2,5)]
fig_8_results = [read_data(path) for path in fig_8_paths]

fig_10_paths = [gen_bp_result_path(10, 10, 2*D2*D2+10, D2, (5,5)) for D2 in range(2,5)]
fig_10_results = [read_data(path) for path in fig_10_paths]


mc_value_8 = -0.6724
mc_value_10 = -0.6699


fontsize = 20
labelsize = 16
linewidth = 3.
markersize = 10

colors = ['c', 'g', 'b', 'y', 'r', 'k']
markers = ['o', 's', '*', '+', '^', '*']



fig, ax = plt.subplots(1, 2, figsize=(10, 4.5))


alphas = [1., 0.6, 0.2]

color_a = 'g'
color_b = 'c'


Ds = [2,3,4]


for (i, item) in enumerate(fig_8_results):
	energies, final_energy = item
	energies = energies[10:]
	xs = linspace(0, 1, len(energies))
	ax[0].plot(xs, energies, color = colors[i],  ls='-', linewidth=linewidth, label='$BP, D=%s$'%(Ds[i]))
	ax[0].plot(xs, [final_energy for x in xs], color = colors[i],  ls='--', linewidth=linewidth, label='$Final, D=%s$'%(Ds[i]))

ax[0].plot(xs, [mc_value_8 for x in xs], color = 'r',  ls='--', linewidth=linewidth, label='$MC$')

ax[0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0].tick_params(axis='y', which='both')
ax[0].locator_params(nbins=6)
ax[0].set_ylabel(r'$E$', fontsize=fontsize, rotation=0)
ax[0].set_xlabel(r'$iterations$', fontsize=fontsize)
ax[0].set_xlim(0,1)
# ax[0].set_ylim(0,1.02)
# ax[0].yaxis.set_label_coords(-0.2, 0.55)
# ax[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
# ax[1, 1].set_ylim(0)
# ax1.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
# ax[1,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0].annotate(r'(a)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[0].legend(loc = 'upper right', fontsize=12)
ax[0].set_title("Heisenberg model (8*8)", fontsize=12)


for (i, item) in enumerate(fig_10_results):
	energies, final_energy = item
	energies = energies[10:]
	xs = linspace(0, 1, len(energies))
	ax[1].plot(xs, energies, color = colors[i],  ls='-', linewidth=linewidth, label='$BP, D=%s$'%(Ds[i]))
	ax[1].plot(xs, [final_energy for x in xs], color = colors[i],  ls='--', linewidth=linewidth, label='$Final, D=%s$'%(Ds[i]))

ax[1].plot(xs, [mc_value_10 for x in xs], color = 'r',  ls='--', linewidth=linewidth, label='$MC$')


ax[1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1].tick_params(axis='y', which='both')
ax[1].locator_params(nbins=6)
ax[1].set_ylabel(r'$E$', fontsize=fontsize, rotation=0)
ax[1].set_xlabel(r'$iterations$', fontsize=fontsize)
ax[1].set_xlim(0,1)
# ax[1].set_ylim(0,0.035)
# ax[1].yaxis.set_label_coords(-0.2, 0.55)
ax[1].annotate(r'(b)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[1].legend(loc = 'upper right', fontsize=12)
ax[1].set_title("Heisenberg model (10*10)", fontsize=12)

plt.tight_layout(pad=0.2)

plt.savefig('reproduce_table2.pdf', dpi=150)

plt.show()


