import json
from matplotlib import pyplot as plt
from matplotlib import rc
from numpy import linspace

rc('text', usetex=True)

def read_data(data_path):
	with open(data_path, 'r') as f:
		data = f.read()
		data = json.loads(data)
	return data['energies']

def gen_result_path(m, n, D1, D2):
	return "result/Table_iv_m%s_n%s_D1%s_D2%s.json"%(m, n, D1, D2)

def gen_bp_result_path(m, n, D1, D2, block_size):
	return "result/Table_iv_bp_m%s_n%s_blocksize%s_%s_D1%s_D2%s.json"%(m, n, block_size[0], block_size[1], D1, D2)



m = 10
n = 10

figa_fu_paths = [gen_result_path(m, n, 20, 2), gen_result_path(m, n, 50, 3), gen_result_path(m, n, 50, 4)]
figa_fu_results = [read_data(path) for path in figa_fu_paths]

block_size = (5, 5)
figa_bp_paths = [gen_bp_result_path(m, n, 20, 2, block_size), gen_bp_result_path(m, n, 50, 3, block_size)]
figa_bp_results = [read_data(path) for path in figa_bp_paths]

m = 14
n = 14

figb_fu_paths = [gen_result_path(m, n, 20, 2), gen_result_path(m, n, 50, 3), gen_result_path(m, n, 50, 4)]
figb_fu_results = [read_data(path) for path in figb_fu_paths]

block_size = (7, 7)
figb_bp_paths = [gen_bp_result_path(m, n, 20, 2, block_size), gen_bp_result_path(m, n, 50, 3, block_size)]
figb_bp_results = [read_data(path) for path in figb_bp_paths]


fontsize = 20
labelsize = 16
linewidth = 3.
markersize = 10

colors = ['c', 'g', 'b', 'y', 'r', 'k']
markers = ['o', 's', '*', '+', '^', '*']



fig, ax = plt.subplots(1, 2, figsize=(10, 5))


alphas = [1., 0.6, 0.2]

color_a = 'g'
color_b = 'c'


Ds = [2,3,4,5]


for (i, energies) in enumerate(figa_fu_results):
	energies = energies[10:]
	xs = linspace(0, 1, len(energies))
	ax[0].plot(xs, energies, color = colors[i],  ls='--', linewidth=3.5, label='$FU, D=%s$'%(Ds[i]))


for (i, energies) in enumerate(figa_bp_results):
	energies = energies[10:]
	xs = linspace(0, 1, len(energies))
	ax[0].plot(xs, energies, color = colors[i], ls='-', linewidth=2.5, label='$BP(5 * 5), D=%s$'%(Ds[i]))


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
ax[0].set_title("Heisenberg model (10*10)", fontsize=12)


for (i, energies) in enumerate(figb_fu_results):
	energies = energies[10:]
	xs = linspace(0, 1, len(energies))
	ax[1].plot(xs, energies, color = colors[i], ls='--', linewidth=3.5, label='$FU, D=%s$'%(Ds[i]))
	# ax[0].plot(a, color = colors[i], ls='-', linewidth=linewidth, label='$\Gamma=%s$'%(Gamma))

for (i, energies) in enumerate(figb_bp_results):
	energies = energies[10:]
	xs = linspace(0, 1, len(energies))
	ax[1].plot(xs, energies, color = colors[i], ls='-', linewidth=2.5, label='$BP(7 * 7), D=%s$'%(Ds[i]))


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
ax[1].set_title("Heisenberg model (14*14)", fontsize=12)

plt.tight_layout(pad=0.2)

plt.savefig('reproduce_table_iv.pdf', dpi=150)

plt.show()


