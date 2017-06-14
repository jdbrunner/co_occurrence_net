from pylab import *

aa = randint(10,size = (4,4))
print(aa)
print(sum(aa))

while aa.shape[0] >1:
	mask = array(aa)
	fill_diagonal(mask, 0)
	maxA = unravel_index(mask.argmax(),aa.shape)
	print(maxA)
	aa[maxA[0]] = aa[maxA[0]]+aa[maxA[1]]
	aa = delete(aa,maxA[1],0)
	print(aa)
	aa[:,maxA[0]] = aa[:,maxA[0]] + aa[:,maxA[1]]
	aa = delete(aa,maxA[1],1)
	print(aa)