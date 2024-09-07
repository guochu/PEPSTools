import json
from matplotlib import pyplot as plt
from matplotlib import rc
from numpy import linspace, asarray

rc('text', usetex=True)

def read_data(data_path):
	with open(data_path, 'r') as f:
		data = f.read()
		data = json.loads(data)
	return asarray(data['bp_energies']), asarray(data['bp_energies_central']), asarray(data["bo_energies"])

def gen_bp_result_path(m, n, D2):
	return "result/periodic_expectations_comparison_m%s_n%s_D%s.json"%(m, n, D2)


m = 10
n = 10
D2 = 2

fontsize = 16
labelsize = 14
linewidth = 3.
markersize = 10

colors = ['c', 'g', 'b', 'y', 'r', 'k']
markers = ['o', 's', '*', '+', '^', '*']



fig, ax = plt.subplots(2, 2, figsize=(8, 7))


data_path = gen_bp_result_path(m, n, D2)

bp_energies, bp_energies_central, bo_energies = read_data(data_path)

xs = range(1, len(bp_energies)+1)

i = 0
ax[0,0].plot(xs, bp_energies, color = colors[i], marker=markers[i], markersize=markersize,  ls='--',  linewidth=linewidth, markerfacecolor='none', label='BP-Naive')

i = 1
ax[0,0].plot(xs, bp_energies_central, color = colors[i], marker=markers[i], markersize=markersize,  ls='--',  linewidth=linewidth, markerfacecolor='none', label='BP-Central')

i = 2
ax[0,0].plot(xs, bo_energies, color = colors[i], marker=markers[i], markersize=markersize,  ls='--',  linewidth=linewidth, markerfacecolor='none', label='BMPO')

ax[0,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,0].tick_params(axis='y', which='both')
ax[0,0].locator_params(nbins=6)
ax[0,0].set_ylabel(r'$E$', fontsize=fontsize, rotation=90)
ax[0,0].set_xlabel(r'$sweeps$', fontsize=fontsize)

# ax[0].set_ylim(0,1.02)
# ax[0].yaxis.set_label_coords(-0.2, 0.55)
# ax[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
# ax[1, 1].set_ylim(0)
# ax1.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
# ax[1,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,0].annotate(r'(a)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[0,0].legend(loc = 'upper right', fontsize=fontsize)
ax[0,0].set_title("Heisenberg model (D=2)", fontsize=fontsize)

diff_bp = abs((bp_energies - bo_energies) /  bo_energies)
diff_bp_central = abs((bp_energies_central - bo_energies) /  bo_energies)


i = 0
ax[0,1].plot(xs, diff_bp, color = colors[i], marker=markers[i], markersize=markersize,  ls='--',  linewidth=linewidth, markerfacecolor='none', label='BP-Naive')

i = 1
ax[0,1].plot(xs, diff_bp_central, color = colors[i], marker=markers[i], markersize=markersize,  ls='--',  linewidth=linewidth, markerfacecolor='none', label='BP-Central')


ax[0,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,1].tick_params(axis='y', which='both')
ax[0,1].locator_params(nbins=6)
ax[0,1].set_ylabel(r'$|E/E_{BM}-1|$', fontsize=fontsize, rotation=90)
ax[0,1].set_xlabel(r'$sweeps$', fontsize=fontsize)
# ax[1].set_xlim(0,1)
# ax[1].set_ylim(0,0.035)
# ax[1].yaxis.set_label_coords(-0.2, 0.55)
ax[0,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,1].annotate(r'(b)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[0,1].legend(loc = 'upper right', fontsize=fontsize)
ax[0,1].set_title("Heisenberg model (D=2)", fontsize=fontsize)



# D = 3
bp_energies = asarray([-0.661668,-0.6630726, -0.663929, -0.664472, -0.6648176, -0.665035, -0.6651669])
bp_energies_central = asarray([-0.661296679, -0.6625885, -0.66334454, -0.663799, -0.6640684, -0.6642193, -0.6642933])
bo_energies = asarray([-0.66127246, -0.66254316, -0.6632779886, -0.6637164969, -0.6639687, -0.6640997, -0.66416177])

xs = range(1, len(bp_energies)+1)

i = 0
ax[1,0].plot(xs, bp_energies, color = colors[i], marker=markers[i], markersize=markersize,  ls='--',  linewidth=linewidth, markerfacecolor='none', label='BP-Naive')

i = 1
ax[1,0].plot(xs, bp_energies_central, color = colors[i], marker=markers[i], markersize=markersize,  ls='--',  linewidth=linewidth, markerfacecolor='none', label='BP-Central')

i = 2
ax[1,0].plot(xs, bo_energies, color = colors[i], marker=markers[i], markersize=markersize,  ls='--',  linewidth=linewidth, markerfacecolor='none', label='BMPO')

ax[1,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,0].tick_params(axis='y', which='both')
ax[1,0].locator_params(nbins=6)
ax[1,0].set_ylabel(r'$E$', fontsize=fontsize, rotation=90)
ax[1,0].set_xlabel(r'$sweeps$', fontsize=fontsize)

ax[1,0].annotate(r'(c)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[1,0].legend(loc = 'upper right', fontsize=fontsize)
ax[1,0].set_title("Heisenberg model (D=3)", fontsize=fontsize)


diff_bp = abs((bp_energies - bo_energies) /  bo_energies)
diff_bp_central = abs((bp_energies_central - bo_energies) /  bo_energies)


i = 0
ax[1,1].plot(xs, diff_bp, color = colors[i], marker=markers[i], markersize=markersize,  ls='--',  linewidth=linewidth, markerfacecolor='none', label='BP-Naive')

i = 1
ax[1,1].plot(xs, diff_bp_central, color = colors[i], marker=markers[i], markersize=markersize,  ls='--',  linewidth=linewidth, markerfacecolor='none', label='BP-Central')


ax[1,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,1].tick_params(axis='y', which='both')
ax[1,1].locator_params(nbins=6)
ax[1,1].set_ylabel(r'$|E/E_{BM}-1|$', fontsize=fontsize, rotation=90)
ax[1,1].set_xlabel(r'$sweeps$', fontsize=fontsize)
# ax[1].set_xlim(0,1)
# ax[1].set_ylim(0,0.035)
# ax[1].yaxis.set_label_coords(-0.2, 0.55)
ax[1,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,1].annotate(r'(d)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[1,1].legend(loc = 'upper right', fontsize=fontsize)
ax[1,1].set_title("Heisenberg model (D=3)", fontsize=fontsize)


plt.tight_layout(pad=0.5)

plt.savefig('expec_comparison.pdf', dpi=150)

plt.show()
