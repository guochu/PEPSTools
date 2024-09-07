import json
from matplotlib import pyplot as plt
from matplotlib import rc
from numpy import linspace, asarray
from math import log, sqrt, sinh

rc('text', usetex=True)

def Onsager_m(beta):
	beta_C =  log(1+sqrt(2.0))/2
	if beta <= beta_C:
		exact=0.0
	else:
		exact = ( 1-1/(sinh(2*beta))**4)**(1.0/8)
	return exact


def read_data(data_path):
	with open(data_path, 'r') as f:
		data = f.read()
		data = json.loads(data)

	final_point = 31
	return abs(asarray(data['bp']))[:final_point], asarray(data["Onsager"])[:final_point], asarray(data["betas"])[:final_point]
	

def gen_classical_tn_path(N, h, block_size, D):
	return "result/infinite_ising_tn_N%s_h%s_blocksize%s_D%s.json"%(N, h, block_size, D)

def gen_classical_critical_tn_path(N, h, block_size, D):
	return "result/infinite_ising_critical_tn_N%s_h%s_blocksize%s_D%s.json"%(N, h, block_size, D)

fontsize = 16
labelsize = 14
linewidth = 2.5
markersize = 8

colors = ['c', 'g', 'b', 'y', 'r', 'k']
markers = ['o', 's', '>', '+', '^', '*']


h = 0.0

D = 18


fig, ax = plt.subplots(1, 1, figsize=(7, 3))

Ns = [5,15,101]
Ds = [18,18,18]

i = 0
data_path = gen_classical_tn_path(Ns[i], h, Ns[i], Ds[i])
bp1, Onsager, betas = read_data(data_path)


print(betas)


# critical_start = 6
# critical_end = 17

# critical_betas = betas[critical_start:critical_end]

# print(critical_betas)

ax.plot(betas, bp1, color = colors[i], marker=markers[i], ls='--', markersize=markersize, linewidth=linewidth, markerfacecolor='none', label=r'$%s\times %s$'%(Ns[i], Ns[i]))

i = 1
data_path = gen_classical_tn_path(Ns[i], h, Ns[i], Ds[i])
bp2, Onsager, betas = read_data(data_path)


ax.plot(betas, bp2, color = colors[i], marker=markers[i], ls='--', markersize=markersize, linewidth=linewidth, markerfacecolor='none', label=r'$%s\times %s$'%(Ns[i], Ns[i]))


i = 2
data_path = gen_classical_tn_path(Ns[i], h, Ns[i], 8)
bp3, Onsager, betas = read_data(data_path)


ax.plot(betas, bp3, color = colors[i], marker=markers[i], ls='--', markersize=markersize, linewidth=linewidth, markerfacecolor='none', label=r'$%s\times %s$'%(Ns[i], Ns[i]))

Onsager_betas = linspace(betas[0], betas[-1], 100)
Onsager = [Onsager_m(beta) for beta in Onsager_betas]

ax.plot(Onsager_betas, Onsager, color = 'r', ls='-', markersize=markersize, linewidth=1.5, markerfacecolor='none', label=r'Onsager')


ax.tick_params(axis='both', which='major', labelsize=labelsize)
ax.tick_params(axis='y', which='both')
ax.locator_params(axis='x', nbins=10)

ax.set_ylabel(r'$\bar{m}$', fontsize=fontsize)
ax.set_xlabel(r'$\beta$', fontsize=fontsize)



ax1 = ax.inset_axes([0.5, 0.3, 0.45, 0.45])

marker_size_s = 5
fontsize_s = 14
labelsize_s = 12
linewidth_s = 2.


i = 0
data_path = gen_classical_critical_tn_path(Ns[i], h, Ns[i], Ds[i])
bp1, Onsager, betas = read_data(data_path)

print(betas)

ax1.plot(betas, bp1, color = colors[i], marker=markers[i], ls='--', markersize=marker_size_s, linewidth=linewidth_s, markerfacecolor='none')

i = 1
data_path = gen_classical_critical_tn_path(Ns[i], h, Ns[i], Ds[i])
bp2, Onsager, betas = read_data(data_path)

ax1.plot(betas, bp2, color = colors[i], marker=markers[i], ls='--', markersize=marker_size_s, linewidth=linewidth_s, markerfacecolor='none')

i = 2
data_path = gen_classical_critical_tn_path(Ns[i], h, Ns[i], Ds[i])
bp3, Onsager, betas = read_data(data_path)

ax1.plot(betas, bp3, color = colors[i], marker=markers[i], ls='--', markersize=marker_size_s, linewidth=linewidth_s, markerfacecolor='none')

Onsager_betas = linspace(betas[0], betas[-1], 100)
Onsager = [Onsager_m(beta) for beta in Onsager_betas]

ax1.plot(Onsager_betas, Onsager, color = 'r', ls='-', markersize=marker_size_s, linewidth=1.5, markerfacecolor='none')

ax1.tick_params(axis='both', which='major', labelsize=labelsize_s)
ax1.tick_params(axis='y', which='both')
ax1.locator_params(axis='x', nbins=6)



ax.legend(loc='upper left', fontsize=labelsize)

plt.tight_layout(pad=0.5)

# plt.savefig('ising_magnetization.pdf', dpi=150)

plt.show()

