import json
from matplotlib import pyplot as plt
from matplotlib import rc
from numpy import linspace, asarray, log

rc('text', usetex=True)


def read_data(data_path):
	with open(data_path, 'r') as f:
		data = f.read()
	data = json.loads(data)
	return data['Bs'][:-4], data['ms'][:-4]

def read_critical_data(data_path):
	with open(data_path, 'r') as f:
		data = f.read()
	data = json.loads(data)
	return data['Bs'], data['blocksizes'], data['ms']

D2 = 2
D1 = 2*D2*D2 + 10
data_path = 'result/ising_pbc_mag_obc_m2_D1_%s_D2_%s.json'%(D1, D2)


Bs2, ms2_full = read_data(data_path)
ms2 = [-item[-1] for item in ms2_full]

critical_data_path = 'result/ising_pbc_critical_mag_obc_m2_D1_%s_D2_%s.json'%(D1, D2)


Bs2_critical, blocksizes, ms2_critical_full = read_critical_data(critical_data_path)
ms2_critical_full = asarray(ms2_critical_full)

# print(Bs2_critical)

# data_path = 'result/ising_pbc_mag_obc_m4_D1_%s_D2_%s.json'%(D1, D2)


# Bs3, ms3 = read_data(data_path)
# ms3 = [-item[-1] for item in ms3]

D2 = 3
D1 = 2*D2*D2 + 10
data_path = 'result/ising_pbc_mag_obc_m2_D1_%s_D2_%s.json'%(D1, D2)


Bs3, ms3_full = read_data(data_path)
ms3 = [-item[-1] for item in ms3_full]

# print(ms)

linewidth = 3
fontsize = 18
labelsize = 14
markersize = 8

colors = ['c', 'g', 'b', 'y', 'r', 'k']
markers = ['o', 's', '>', '+', '^', '*']



fig, ax = plt.subplots(1, 1, figsize=(7, 3))


ax.plot(Bs2, ms2, color='c', ls='--', linewidth=linewidth, markersize=markersize, marker='o', markerfacecolor='none')

ax.plot(Bs3, ms3, color='b', ls='--', linewidth=linewidth, markersize=markersize, marker='s', markerfacecolor='none')

ax.tick_params(axis='both', which='major', labelsize=labelsize)
ax.tick_params(axis='y', which='both')
ax.locator_params(axis='x', nbins=10)


ax1 = ax.inset_axes([0.1, 0.2, 0.4, 0.4])

marker_size_s = 5
fontsize_s = 14
labelsize_s = 12
linewidth_s = 2.

print(Bs2_critical)

Bs2_critical_pos = [7]

print([Bs2_critical[item] for item in Bs2_critical_pos])

for i in range(len(Bs2_critical_pos)):
	ax1.plot(blocksizes, -ms2_critical_full[Bs2_critical_pos[i], :], color='c', ls='--', linewidth=linewidth_s, markersize=marker_size_s, marker=markers[i], markerfacecolor='none', label=r'B=%s'%(Bs2_critical[Bs2_critical_pos[i]]))

# for i in range(len(blocksizes)):
# 	ax1.plot(Bs2_critical, -ms2_critical_full[:,i], ls='--', linewidth=linewidth_s, markersize=marker_size_s, marker=markers[i], markerfacecolor='none')


# ax.set_ylim(top=0.6)

ax1.tick_params(axis='both', which='major', labelsize=labelsize_s)
ax1.tick_params(axis='y', which='both')
ax1.locator_params(axis='x', nbins=10)
# ax1.set_title(r'$B=%s$'%(Bs2_critical_pos[0]), fontsize=fontsize_s)

ax1.legend(fontsize=fontsize_s)

plt.tight_layout(pad=0.5)

plt.show()


