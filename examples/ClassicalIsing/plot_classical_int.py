import json
from matplotlib import pyplot as plt
from matplotlib import rc
from numpy import linspace, asarray

rc('text', usetex=True)

def read_data(data_path, isperiodic):
	with open(data_path, 'r') as f:
		data = f.read()
		data = json.loads(data)
	if isperiodic:
		return asarray(data['boundarympo']), abs(asarray(data['bp'])), asarray(data["Onsager"]), asarray(data["betas"])
	else:
		return asarray(data['boundarymps']), abs(asarray(data['bp'])), asarray(data["Onsager"]), asarray(data["betas"])

def gen_classical_tn_path(N, h, block_size, D, isperiodic):
	if isperiodic:
		return "result/classical_periodic_tn_int_N%s_h%s_blocksize%s_D%s.json"%(N, h, block_size, D)
	else:
		return "result/classical_non_periodic_tn_int_N%s_h%s_blocksize%s_D%s.json"%(N, h, block_size, D)




block_size = 10

fontsize = 20
labelsize = 16
linewidth = 3.
markersize = 10

colors = ['c', 'g', 'b', 'y', 'r', 'k']
markers = ['o', 's', '*', '+', '^', '*']



fig, ax = plt.subplots(1, 1, figsize=(6, 5))


N = 20
D = 13

h = 0.001
data_path = gen_classical_tn_path(N, h, block_size, D, isperiodic=False)
boundarymps, bp, Onsager, betas = read_data(data_path, isperiodic=False)
diff = abs((bp - boundarymps) / boundarymps)

i = 0
ax.semilogy(betas, diff, color = colors[i], marker=markers[i], markersize=markersize,  ls='--',  linewidth=linewidth, markerfacecolor='none', label='h=%s'%(h))

h = 0.01
data_path = gen_classical_tn_path(N, h, block_size, D, isperiodic=False)
boundarymps, bp, Onsager, betas = read_data(data_path, isperiodic=False)
diff = abs((bp - boundarymps) / boundarymps)

i += 1
ax.semilogy(betas, diff, color = colors[i], marker=markers[i], markersize=markersize,  ls='--',  linewidth=linewidth, markerfacecolor='none', label='h=%s'%(h))


ax.tick_params(axis='both', which='major', labelsize=labelsize)
ax.tick_params(axis='y', which='both')
# ax.locator_params(nbins=6)
ax.set_ylabel(r'$|m_{BP}/m_{BMPS} - 1|$', fontsize=fontsize, rotation=90)
ax.set_xlabel(r'$\beta$', fontsize=fontsize)
ax.set_ylim(1.0e-10, 10)

# ax[0].set_ylim(0,1.02)
# ax[0].yaxis.set_label_coords(-0.2, 0.55)
# ax[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
# ax[1, 1].set_ylim(0)
# ax1.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
# ax[1,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
# ax.annotate(r'(a)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax.legend(fontsize=fontsize)
# ax[0,0].set_title("ising model (OBC-20)", fontsize=fontsize)





plt.tight_layout(pad=0.5)

plt.savefig('classical_tn_interaction.pdf', dpi=150)

plt.show()
