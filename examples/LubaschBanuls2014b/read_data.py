from numpy import load

data = load('pepe_D_2_B_25.npz', allow_pickle=True)
peps = data['peps']


print(peps.shape)
print(peps[0,0].shape)

print(peps[0,20].shape)