import json
from matplotlib import pyplot as plt
from matplotlib import rc
from numpy import linspace, asarray

rc('text', usetex=True)

def read_data(data_path):
	with open(data_path, "r") as f:
		data = f.read()
	data = json.loads(data)
	return data["final_energies"]

def get_nonzeros_v(v):
	return [item for item in v if abs(item) > 1.0e-12]

def get_nonzeros_vv(vv):
	r = []
	for item in vv:
		rj = get_nonzeros_v(item)
		r.extend(rj)
	return r

def get_nonzeros_item(item):
	return get_nonzeros_vv(item['H']) + get_nonzeros_vv(item['V'])

def get_nonzeros(d):
	r = []
	for item in d:
		r.extend(get_nonzeros_item(item))
	return r


D2 = 4
D1 = 2*D2*D2 + 10
data_path = "result/FinalEnergies_m16_n16_blocksize8_8_D1%s_D2%s.json"%(D1, D2)


data = read_data(data_path)

nonzeros = get_nonzeros(data)

fig, ax = plt.subplots(1, 1, figsize=(8, 7))


xs = [i for i in range(1, len(nonzeros)+1)]

ax.scatter(xs, nonzeros, label=r'bond energy')
mean_value = asarray(nonzeros).mean()
std = asarray(nonzeros).std()
# ax.plot(xs, [mean_value for i in range(len(xs))], ls = '--', linewidth=3, label='')

textstr = '\n'.join((r'$\bar{E}=%s$'%(mean_value), r'$\sigma =%s$'%(std)))

ax.text(250, -0.33, textstr, fontsize=14, ma='left',  \
va='center', ha='center', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.2))

ax.legend(fontsize=14)

plt.tight_layout(pad=0.5)

# plt.savefig('classical_tn.pdf', dpi=150)

plt.show()
