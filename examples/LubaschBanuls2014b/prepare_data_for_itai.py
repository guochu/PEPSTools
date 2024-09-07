from julia import Serialization
from julia import Main

from numpy import savez, empty, ndarray, load

Main.eval('push!(LOAD_PATH, "../../src")')
Main.eval("using PEPSTools")

def gen_peps_path(D2, B):
	return "data/Table_vii_D2%s_B%s.peps"%(D2, B)

def read_data(peps_path):
	# PEPSTools_path = "/Users/guochu/Documents/QuantumSimulator/PEPSTools/src"
	# Main.eval("push!(LOAD_PATH, %s)"%(PEPSTools_path))
	peps = Serialization.deserialize(peps_path)
	peps = peps.data
	r = empty((21,21), dtype = ndarray)
	for i in range(21):
		for j in range(21):
			r[i,j] = peps[i][j]
	return r



peps_path = gen_peps_path(2, 2.5)

peps = read_data(peps_path)

savez("pepe_D_2_B_25", peps=peps)



# data = load('pepe_D_2_B_25.npz', allow_pickle=True)
# peps = data['peps']

# print(peps.shape)

# print(peps[0,0])