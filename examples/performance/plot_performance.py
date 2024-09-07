import json
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.colors import LogNorm
from numpy import linspace, asarray

rc('text', usetex=True)

mpath = 'result/'

def read_data(data_path):
	with open(data_path, 'r') as f:
		data = f.read()
	data = json.loads(data)
	return asarray(data['times']).T, asarray(data['energies'])


def gen_bmps_path(D2):
	D1 = 2*D2*D2 + 10
	return mpath + 'bmps_m40_D2_%s_D1_%s.json'%(D2, D1)

def gen_serial_bp_path(D2):
	D1 = 2*D2*D2 + 10
	return mpath + 'serial_bp_m40_blocksize_5_D2_%s_D1_%s.json'%(D2, D1)

def gen_parallel_bp_path(D2, np):
	D1 = 2*D2*D2 + 10
	return mpath + 'parallel_bp_np_%s_m40_blocksize_5_D2_%s_D1_%s.json'%(np, D2, D1)


def average_time(ts):
	return asarray(ts[1:]).mean()

D2 = 2

data_path = gen_bmps_path(D2)

ts, energies = read_data(data_path)

bmps_t = average_time(ts)
# print(energies)


data_path = gen_serial_bp_path(D2)

ts, energies = read_data(data_path)

bp_t = average_time(ts)

# print(energies)


np = [2,4,8,16,32]

parallal_data = [read_data(gen_parallel_bp_path(D2, p)) for p in np]
parallal_ts = [average_time(ts) for (ts, energies) in parallal_data]

np = asarray([1] + np)
parallal_ts = [bp_t] + parallal_ts

print(np)
print(parallal_ts)

linewidth = 3

fontsize = 18
labelsize = 14
markersize = 8

fig, ax = plt.subplots(1, 1, figsize=(6, 5))

ax.plot(1 / np, parallal_ts, color = 'g', ls='--', marker='s', markersize=markersize, linewidth=linewidth)

# xs = [i for i in range(1, 33)]

xs = linspace(0, 1, 100)

ax.plot(xs, [bmps_t for x in xs], color='r', ls='-', linewidth=linewidth)

ax.set_xlim(left=0)

ax.tick_params(axis='both', which='major', labelsize=labelsize)
ax.tick_params(axis='y', which='both')
# ax.locator_params(nbins=8)
ax.set_ylabel(r'$t (s)$', fontsize=fontsize, rotation=90)
ax.set_xlabel(r'$1 /n$', fontsize=fontsize)

plt.tight_layout(pad=0.5)


plt.show()



