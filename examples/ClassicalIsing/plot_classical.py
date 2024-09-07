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
		return "result/classical_periodic_tn_N%s_h%s_blocksize%s_D%s.json"%(N, h, block_size, D)
	else:
		return "result/classical_non_periodic_tn_N%s_h%s_blocksize%s_D%s.json"%(N, h, block_size, D)




block_size = 10
h = 0.001

fontsize = 16
labelsize = 14
linewidth = 3.
markersize = 8

colors = ['c', 'g', 'b', 'y', 'r', 'k']
markers = ['o', 's', '*', '+', '^', '*']



fig, ax = plt.subplots(2, 2, figsize=(8, 7))


N = 20
D = 13

data_path = gen_classical_tn_path(N, h, block_size, D, isperiodic=False)
boundarymps, bp, Onsager, betas = read_data(data_path, isperiodic=False)


i = 0
ax[0,0].plot(betas, boundarymps, color = colors[i], marker=markers[i], markersize=markersize,  ls='--',  linewidth=linewidth, markerfacecolor='none', label='BMPS-20')

i = 1
ax[0,0].plot(betas, bp, color = colors[i], marker=markers[i], markersize=markersize,  ls='--',  linewidth=linewidth, markerfacecolor='none', label='BP-20-10')

i = 2
ax[0,0].plot(betas, Onsager, color = colors[i], ls='-',  linewidth=linewidth, markerfacecolor='none', label='Onsager')

ax[0,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,0].tick_params(axis='y', which='both')
ax[0,0].locator_params(nbins=6)
ax[0,0].set_ylabel(r'$\langle m\rangle$', fontsize=fontsize, rotation=90)
ax[0,0].set_xlabel(r'$\beta$', fontsize=fontsize)

# ax[0].set_ylim(0,1.02)
# ax[0].yaxis.set_label_coords(-0.2, 0.55)
# ax[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
# ax[1, 1].set_ylim(0)
# ax1.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
# ax[1,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,0].annotate(r'(a)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[0,0].legend(loc = 'lower right', fontsize=fontsize)
# ax[0,0].set_title("ising model (OBC-20)", fontsize=fontsize)


data_path = gen_classical_tn_path(N, h, block_size, D, isperiodic=True)
boundarymps, bp, Onsager, betas = read_data(data_path, isperiodic=True)

# print(Onsager)

i = 0
ax[0,1].plot(betas, boundarymps, color = colors[i], marker=markers[i], markersize=markersize,  ls='--',  linewidth=linewidth, markerfacecolor='none', label='BMPO-20')

i = 1
ax[0,1].plot(betas, bp, color = colors[i], marker=markers[i], markersize=markersize,  ls='--',  linewidth=linewidth, markerfacecolor='none', label='BP-20-10')

i = 2
ax[0,1].plot(betas, Onsager, color = colors[i], ls='-',  linewidth=linewidth, markerfacecolor='none', label='Onsager')

ax[0,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,1].tick_params(axis='y', which='both')
ax[0,1].locator_params(nbins=6)
ax[0,1].set_ylabel(r'$\langle m\rangle$', fontsize=fontsize, rotation=90)
ax[0,1].set_xlabel(r'$\beta$', fontsize=fontsize)

# ax[0].set_ylim(0,1.02)
# ax[0].yaxis.set_label_coords(-0.2, 0.55)
# ax[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
# ax[1, 1].set_ylim(0)
# ax1.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
# ax[1,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,1].annotate(r'(b)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[0,1].legend(loc = 'lower right', fontsize=fontsize)
# ax[0,1].set_title("ising model (PBC-20)", fontsize=fontsize)


# 

N = 40

data_path = gen_classical_tn_path(N, h, block_size, D, isperiodic=False)
boundarymps, bp, Onsager, betas = read_data(data_path, isperiodic=False)


i = 0
ax[1,0].plot(betas, boundarymps, color = colors[i], marker=markers[i], markersize=markersize,  ls='--',  linewidth=linewidth, markerfacecolor='none', label='BMPS-40')

i = 1
ax[1,0].plot(betas, bp, color = colors[i], marker=markers[i], markersize=markersize,  ls='--',  linewidth=linewidth, markerfacecolor='none', label='BP-40-10')

i = 2
ax[1,0].plot(betas, Onsager, color = colors[i], ls='-',  linewidth=linewidth, markerfacecolor='none', label='Onsager')

ax[1,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,0].tick_params(axis='y', which='both')
ax[1,0].locator_params(nbins=6)
ax[1,0].set_ylabel(r'$\langle m\rangle$', fontsize=fontsize, rotation=90)
ax[1,0].set_xlabel(r'$\beta$', fontsize=fontsize)

# ax[0].set_ylim(0,1.02)
# ax[0].yaxis.set_label_coords(-0.2, 0.55)
# ax[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
# ax[1, 1].set_ylim(0)
# ax1.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
# ax[1,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,0].annotate(r'(a)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[1,0].legend(loc = 'lower right', fontsize=fontsize)
# ax[1,0].set_title("ising model (OBC-40)", fontsize=fontsize)


data_path = gen_classical_tn_path(N, h, block_size, D, isperiodic=True)
boundarymps, bp, Onsager, betas = read_data(data_path, isperiodic=True)

# print(Onsager)

i = 0
ax[1,1].plot(betas, boundarymps, color = colors[i], marker=markers[i], markersize=markersize,  ls='--',  linewidth=linewidth, markerfacecolor='none', label='BMPO-40')

i = 1
ax[1,1].plot(betas, bp, color = colors[i], marker=markers[i], markersize=markersize,  ls='--',  linewidth=linewidth, markerfacecolor='none', label='BP-40-10')

i = 2
ax[1,1].plot(betas, Onsager, color = colors[i], ls='-',  linewidth=linewidth, markerfacecolor='none', label='Onsager')

ax[1,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,1].tick_params(axis='y', which='both')
ax[1,1].locator_params(nbins=6)
ax[1,1].set_ylabel(r'$\langle m\rangle$', fontsize=fontsize, rotation=90)
ax[1,1].set_xlabel(r'$\beta$', fontsize=fontsize)

# ax[0].set_ylim(0,1.02)
# ax[0].yaxis.set_label_coords(-0.2, 0.55)
# ax[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
# ax[1, 1].set_ylim(0)
# ax1.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
# ax[1,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,1].annotate(r'(b)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[1,1].legend(loc = 'lower right', fontsize=fontsize)
# ax[0,1].set_title("ising model (PBC-20)", fontsize=fontsize)


plt.tight_layout(pad=0.5)

plt.savefig('classical_tn.pdf', dpi=150)

plt.show()
