import json
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
from numpy import linspace, asarray
from math import log, sqrt, sinh

rc('text', usetex=True)


def read_critical_data(data_path):
	mpath = '/Users/guochu/Documents/QuantumSimulator/RectPEPS/examples/MezzacapoCirac2009/'
	with open(mpath+data_path, 'r') as f:
		data = f.read()
	data = json.loads(data)
	return data['Bs'], data['blocksizes'], data['ms']



linewidth = 3

fontsize = 16
labelsize = 12
markersize = 8

colors = ['lightgray', 'darkgray', 'dimgray', 'black', 'r', 'k']
markers = ['x', 's', '>', 'o', '^', '*']


fig, ax = plt.subplots(1, 1, figsize=(7, 7))


alphas = [0.2, 0.4, 0.7, 1.]


D2 = 3
D1 = 2*D2*D2 + 10
critical_data_path = 'result/ising_pbc_critical_mag_obc_m2_D1_%s_D2_%s.json'%(D1, D2)
Bs3_critical, blocksizes, ms3_critical_full = read_critical_data(critical_data_path)
ms3_critical_full = asarray(ms3_critical_full)

print(Bs3_critical)

Bs3_critical = asarray(Bs3_critical)


for i in range(len(blocksizes)):
	block_size = blocksizes[i]
	ax.plot(Bs3_critical, -ms3_critical_full[:,i], color='g', alpha=alphas[i], ls='--', linewidth=linewidth, markersize=markersize, marker='o', markerfacecolor='none', label=r'BMPS $%s\times %s$'%(block_size,block_size))


critical_data_path = 'result/ising_pbc_critical_mag_pbc_m2_D1_%s_D2_%s.json'%(D1, D2)
Bs3_critical, blocksizes, ms3_critical_full = read_critical_data(critical_data_path)
ms3_critical_full = asarray(ms3_critical_full)

for i in range(len(blocksizes)):
	block_size = blocksizes[i]
	ax.plot(Bs3_critical, -ms3_critical_full[:,i], color='b', alpha=alphas[i], ls=':', linewidth=linewidth, markersize=markersize, marker='x', markerfacecolor='none', label=r'BlockBP $%s\times %s$'%(block_size,block_size))


# ax.set_ylim(top=0.6)

ax.tick_params(axis='both', which='major', labelsize=labelsize)
ax.tick_params(axis='y', which='both')
ax.locator_params(axis='x', nbins=8)

ax.legend(fontsize=labelsize)


plt.tight_layout(pad=0.5)

plt.savefig('check_critical_magnetization.pdf', dpi=150)



plt.show()


