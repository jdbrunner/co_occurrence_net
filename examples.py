############ Small examples

## Jim Brunner
## modules needed
from pylab import *
import pandas as pd
###pandas DataFrame usage: df['column_name'] or df.loc['row_name' (index),'column'] or df.iloc[row_number, 'column']
import sys
import matplotlib.pyplot as plt
from co_occ_funs import *

if __name__ == '__main__':


	tiny = False
	if tiny:
		A1 = array([[0,1,1,1,1,0],
					[1,0,1,0,0,0],
					[1,1,0,0,1,1],
					[1,0,0,0,0,0],
					[1,0,1,0,0,0],
					[0,0,1,0,0,0]])
			
		A2 = array([[0,1,1,1,0,0],
					[1,0,1,0,0,0],
					[1,1,0,0,0,0],
					[1,0,0,0,1,1],
					[0,0,0,1,0,1],
					[0,0,0,1,1,0]])
			
		A3 = array([[0,1,1,1,0,0,0],
					[1,0,1,1,1,0,0],
					[1,1,0,1,1,0,0],
					[1,1,1,0,0,1,0],
					[0,1,1,0,0,0,1],
					[0,0,0,1,0,0,1],
					[0,0,0,0,1,1,0]])
			
		x1 = 1/6
		x2 = 1/2
		x3 = 1/3
		u1 = array([x1,x2,x3,0,0,0])
		tot_ab1 = sum(u1)
		u2 = array([x1,x2,0,0,x3,0])
		tot_ab2 = sum(u2)
		u3 = array([0,0,0,x1,x2,x3])
		tot_ab3 = sum(u3)

		con = 0.8
		geoms = con**array(range(len(u1)))



		u1_rkd_A1 = diffusion_ivp([],[], A1, sample = u1, all = True)
		u1_rkd_A2 = diffusion_ivp([],[], A2, sample = u1, all = True)
	

	
		u1A1 = u1[flat_two_deep(u1_rkd_A1[0])]
		#print(u1_rkd_A1[0])
		val11 = dot(u1A1/tot_ab1,geoms)
		u1A2 = u1[flat_two_deep(u1_rkd_A2[0])]
		#print(u1_rkd_A2[0])
		val12 = dot(u1A2/tot_ab1,geoms)


		u2_rkd_A1 = diffusion_ivp([],[], A1, sample = u2, all = True)
		u2_rkd_A2 = diffusion_ivp([],[], A2, sample = u2, all = True)

		u2A1 = u2[flat_two_deep(u2_rkd_A1[0])]
		#print(u2_rkd_A1[0])
		val21 = dot(u2A1/tot_ab2,geoms)
		u2A2 = u2[flat_two_deep(u2_rkd_A2[0])]
		#print(u2_rkd_A2[0])
		val22 = dot(u2A2/tot_ab2,geoms)

		u3_rkd_A1 = diffusion_ivp([],[], A1, sample = u3, all = True)
		u3_rkd_A2 = diffusion_ivp([],[], A2, sample = u3, all = True)

		u3A1 = u3[flat_two_deep(u3_rkd_A1[0])]
		#print(u3_rkd_A1[0])
		val31 = dot(u3A1/tot_ab3,geoms)
		u3A2 = u3[flat_two_deep(u3_rkd_A2[0])]
		#print(u3_rkd_A2[0])
		#print(u3A2)
		val32 = dot(u3A2/tot_ab3,geoms)


		f, axarr = plt.subplots(3, 2)
	
		axarr[0, 0].plot(range(len(u1)),u1A1)
		axarr[0, 0].set_title('u1,A1', fontsize=10)
		axarr[0, 0].text(0.6, 0.6, val11, fontsize = 10, transform = axarr[0,0].transAxes)
	
		axarr[0, 1].plot(range(len(u1)),u1A2)
		axarr[0, 1].set_title('u1,A2', fontsize=10)
		axarr[0, 1].text(0.6, 0.6, val12, fontsize = 10, transform = axarr[0,1].transAxes)
	
		axarr[1, 0].plot(range(len(u2)),u2A1)
		axarr[1, 0].set_title('u2,A1', fontsize=10)
		axarr[1, 0].text(0.6, 0.6, val21, fontsize = 10, transform = axarr[1,0].transAxes)
	
		axarr[1, 1].plot(range(len(u2)),u2A2)
		axarr[1, 1].set_title('u2,A2', fontsize=10)
		axarr[1, 1].text(0.6, 0.6, val22, fontsize = 10, transform = axarr[1,1].transAxes)
	
		axarr[2, 0].plot(range(len(u3)),u3A1)
		axarr[2, 0].set_title('u3,A1', fontsize=10)
		axarr[2, 0].text(0.6, 0.6, val31, fontsize = 10, transform = axarr[2,0].transAxes)
	
		axarr[2, 1].plot(range(len(u3)),u3A2)
		axarr[2, 1].set_title('u3,A2', fontsize=10)
		axarr[2, 1].text(0.6, 0.6, val32, fontsize = 10, transform = axarr[2,1].transAxes)
	
		f.suptitle('Abundance v Rank', fontsize=20)

		show()
	
	
	
	else:

		sample_name = sys.argv[1]
		rkd_sample = pd.read_csv(sample_name,sep = '\t')
		cols = rkd_sample.columns.values
		ivp = list(where(['ivp' in col for col in cols])[0])
		eq = list(where(['eq' in col for col in cols])[0])
		colo = list(where(['ivp_col' in col for col in cols])[0])
		trans = list(where(['transient' in col for col in cols])[0])
		for eqq in eq:
			ivp.remove(eqq)
		for tr in trans:
			ivp.remove(tr)
		for c in colo:
			ivp.remove(c)
		#print(mean(rkd_sample['abundance'].values[nonzero(rkd_sample['abundance'].values)]))
		#print(len(nonzero(rkd_sample['abundance'].values)[0]))
		fit = True
		if fit:
			tot_ab = sum(rkd_sample['abundance'].values)
			num_netws = len(ivp)
			num_netws2 = num_netws//2
			f2, axarr2 = plt.subplots(num_netws2,2)
			if num_netws2 > 1:
				i = 0
				for x in range(num_netws2):
					for j in range(2):
						ord = rkd_sample[cols[ivp[i]]].values.argsort()
						axarr2[x,j].plot(rkd_sample[cols[ivp[i]]].values[ord], rkd_sample['abundance'].values[ord])
						axarr2[x,j].text(0, 0.7, cols[ivp[i]], fontsize = 10, transform = axarr2[x,j].transAxes)
						geo = 0.5**rkd_sample[cols[ivp[i]]].values[ord]
						valu = dot(geo,rkd_sample['abundance'].values[ord]/tot_ab)
						axarr2[x,j].text(.8, 0.7, valu, fontsize = 10, transform = axarr2[x,j].transAxes)
						print(cols[ivp[i]])
						print(valu)
						i = i+1
			else:
				for x in range(num_netws):
					ord = rkd_sample[cols[ivp[x]]].values.argsort()
					axarr2[x].plot(rkd_sample[cols[ivp[x]]].values[ord], rkd_sample['abundance'].values[ord])
					axarr2[x].text(0.6, 0.7, cols[ivp[x]], fontsize = 10, transform = axarr2[x].transAxes)
					geo = 0.5**rkd_sample[cols[ivp[x]]].values[ord]
					valu = dot(geo,rkd_sample['abundance'].values[ord]/tot_ab)
					axarr2[x].text(0.6, 0.6, valu, fontsize = 10, transform = axarr2[x].transAxes)
					print(cols[ivp[x]])
					print(valu)
			slsh = sample_name.rfind('/')
			f2.suptitle(sample_name[slsh+1:-4])
			show()
	
		choose_sub = True
		if choose_sub:
	# 		c_clusts = unique(rkd_sample['community_cluster'].values)
	# 		s_clusts = unique(rkd_sample['spectral_cluster'].values)
	# 		sample_types = unique(rkd_sample['common_location'].values)	
			ccls = list(where(['community_cluster' in col for col in cols])[0])
			scls = list(where(['spectral_cluster' in col for col in cols])[0])
			types = list(where(['common_location' in col for col in cols])[0])
			rkd_sample[cols[ccls]] = rkd_sample[cols[ccls]].astype(float)
			most_c_clusts = max([len(unique(rkd_sample[net].values[~isnan(rkd_sample[net].values)])) for net in cols[ccls]])
			rkd_sample[cols[scls]] = rkd_sample[cols[scls]].astype(float)
			most_s_clusts = max([len(unique(rkd_sample[net].values[~isnan(rkd_sample[net].values)])) for net in cols[scls]])
			rkd_sample[cols[types]] = rkd_sample[cols[types]].astype(str)
			most_s_types = max([len(unique(rkd_sample[net].values[~(rkd_sample[net].values == 'nan')])) for net in cols[types]])

			community_ranks = pd.DataFrame(index = range(most_c_clusts), columns = cols[ccls+ivp])
			spectral_ranks = pd.DataFrame(index = range(most_s_clusts), columns = cols[scls+ivp])
			sample_type_ranks = pd.DataFrame(index = range(most_s_types), columns = cols[types+ivp])
		
			nets = [word[4:] for word in cols[ivp]]
		
			for net in nets:
				unccls = unique(rkd_sample['community_cluster_'+net].values[~isnan(rkd_sample['community_cluster_'+net].values)])
				community_ranks.loc[:len(unccls)-1,'community_cluster_'+net] = unccls
				for cl in community_ranks['community_cluster_'+net]:
					in_clus = where([cls == cl for cls in rkd_sample['community_cluster_'+net]])
					clus_loc = where(community_ranks['community_cluster_'+net] == cl)
					community_ranks.loc[community_ranks.index[clus_loc],'ivp_'+net] = mean(rkd_sample.loc[rkd_sample.index[in_clus],'ivp_'+net])
			
				unscls = unique(rkd_sample['spectral_cluster_'+net].values[~isnan(rkd_sample['spectral_cluster_'+net].values)])
				spectral_ranks.loc[:len(unscls)-1,'spectral_cluster_'+net] = unscls
				for cl in spectral_ranks['spectral_cluster_'+net]:
					in_clus = where([cls == cl for cls in rkd_sample['spectral_cluster_'+net]])
					clus_loc = where(spectral_ranks['spectral_cluster_'+net] == cl)
					spectral_ranks.loc[spectral_ranks.index[clus_loc],'ivp_'+net] = mean(rkd_sample.loc[rkd_sample.index[in_clus],'ivp_'+net])
			
			
				untypes = unique(rkd_sample['common_location_'+net].values[~(rkd_sample['common_location_'+net].values == 'nan')])
				sample_type_ranks.loc[:len(untypes)-1,'common_location_'+net]= untypes
				for ty in sample_type_ranks['common_location_'+net]:
					this_type = where([typ == ty for typ in rkd_sample['common_location_'+net]])
					type_loc = where(sample_type_ranks['common_location_'+net] == ty)
					rks_here = rkd_sample.loc[rkd_sample.index[this_type],'ivp_'+net]
					sample_type_ranks.loc[sample_type_ranks.index[type_loc],'ivp_'+net] = mean(sort(rks_here)[:int(0.75*len(rks_here))])
				
			for net in nets:
	# 			print(community_ranks.sort_values('ivp_'+net))
	# 			print(spectral_ranks.sort_values('ivp_'+net))
				print(sample_type_ranks[['common_location_'+net,'ivp_'+net]].sort_values('ivp_'+net))
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	