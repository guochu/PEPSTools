import json
from matplotlib import pyplot as plt
from matplotlib import rc
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy import linspace, asarray

rc('text', usetex=True)

def read_data(data_path):
	with open(data_path, "r") as f:
		data = f.read()
	data = json.loads(data)
	return data["final_energies"]

def get_HV(data):
	H = asarray(data[0]['H'])
	V = asarray(data[0]['V'])
	for i in range(1, len(data)):
		H += asarray(data[i]['H'])
		V += asarray(data[i]['V'])

	return H, V


fontsize = 14
labelsize = 12


D2 = 4
D1 = 2*D2*D2 + 10
data_path = "result/FinalEnergies_m16_n16_blocksize8_8_D1%s_D2%s.json"%(D1, D2)

data = read_data(data_path)

H, V = get_HV(data)

fig, ax = plt.subplots(1, 2, figsize=(14, 6))



im1 = ax[0].imshow(H, extent=(1, 16, 1, 16), aspect='auto')
divider = make_axes_locatable(ax[0])
cax = divider.append_axes("right", size="10%", pad=0.05)
cb = plt.colorbar(im1, cax=cax)
cb.ax.tick_params(labelsize=fontsize)
cb.ax.locator_params(nbins=6)
ax[0].set_title(r'Horizontal bonds', fontsize=fontsize)

ax[0].tick_params(axis='both', which='major', labelsize=labelsize)
# ax[0].locator_params(nbins=8)
# ax[0].set_ylabel(r'$\epsilon$', fontsize=fontsize)
# ax[0].set_xlabel(r'$W$', fontsize=fontsize)



im2 = ax[1].imshow(V, extent=(1, 16, 1, 16), aspect='auto')
divider = make_axes_locatable(ax[1])
cax = divider.append_axes("right", size="10%", pad=0.05)
cb = plt.colorbar(im2, cax=cax)
cb.ax.tick_params(labelsize=fontsize)
cb.ax.locator_params(nbins=6)
ax[1].set_title(r'Vertical bonds', fontsize=fontsize)


ax[1].tick_params(axis='both', which='major', labelsize=labelsize)

plt.tight_layout(pad=0.5)

# plt.savefig('classical_tn.pdf', dpi=150)

plt.show()
