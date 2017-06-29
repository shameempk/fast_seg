import numpy as np



def cie_de_2000(lab1,lab2,k_L,k_C,k_H):
	L_1_star,a_1_star,b_1_star=lab1
	L_2_star,a_2_star,b_2_star=lab2
	C_1_star=np.sqrt(np.power(a_1_star,2)+np.power(b_1_star,2))
	C_2_star=np.sqrt(np.power(a_2_star,2)+np.power(b_2_star,2))
	C_bar_star=np.average([C_1_star,C_2_star])
	
	G=0.5*(1-np.sqrt(np.power(C_bar_star,7)/(np.power(C_bar_star,7)+np.power(25,7))))
	
	a_1_dash=(1+G)*a_1_star
	a_2_dash=(1+G)*a_2_star
	C_1_dash=np.sqrt(np.power(a_1_dash,2)+np.power(b_1_star,2))
	C_2_dash=np.sqrt(np.power(a_2_dash,2)+np.power(b_2_star,2))
	h_1_dash=np.degrees(np.arctan2(b_1_star,a_1_dash))
	h_1_dash += (h_1_dash < 0) * 360
	h_2_dash=np.degrees(np.arctan2(b_2_star,a_2_dash))
	h_2_dash += (h_2_dash < 0) * 360
	
	delta_L_dash=L_2_star-L_1_star
	delta_C_dash=C_2_dash-C_1_dash
	delta_h_dash=0.0
	
	if(C_1_dash*C_2_dash):
		if(np.fabs(h_2_dash-h_1_dash)<=180):
			delta_h_dash=h_2_dash-h_1_dash
		elif(h_2_dash-h_1_dash>180):
			delta_h_dash=(h_2_dash-h_1_dash)-360
		elif(h_2_dash-h_1_dash)<-180:
			delta_h_dash=(h_2_dash-h_1_dash)+360
	
	delta_H_dash=2*np.sqrt(C_1_dash*C_2_dash)*np.sin(np.radians(delta_h_dash)/2.0)
	
	L_bar_dash=np.average([L_1_star,L_2_star])
	C_bar_dash=np.average([C_1_dash,C_2_dash])
	h_bar_dash=h_1_dash+h_2_dash
	
	if(C_1_dash*C_2_dash):
		if(np.fabs(h_1_dash-h_2_dash)<=180):
			h_bar_dash=np.average([h_1_dash,h_2_dash])
		else:
			if(h_1_dash+h_2_dash)<360:
				h_bar_dash=(h_1_dash+h_2_dash+360)/2
			else:
				h_bar_dash=(h_1_dash+h_2_dash-360)/2

	T=1-0.17*np.cos(np.radians(h_bar_dash-30))+0.24*np.cos(np.radians(2*h_bar_dash))\
	+0.32*np.cos(np.radians(3*h_bar_dash+6))-0.20*np.cos(np.radians(4*h_bar_dash-63))
	
	delta_theta=30 * np.exp(- np.power( (h_bar_dash-275) / 25, 2))
	
	R_c=2*np.sqrt( np.power(C_bar_dash,7) / (np.power(C_bar_dash,7)+np.power(25,7)) )
	
	S_L=1+((0.015*np.power(L_bar_dash-50,2))/np.sqrt(20+np.power(L_bar_dash-50,2)))
	S_C=1+0.045*C_bar_dash
	S_H=1+0.015*C_bar_dash*T
	R_T=-R_c * np.sin(2*np.radians(delta_theta))
	
	return delta_e=np.sqrt(np.power(delta_L_dash/(k_L*S_L),2)+\
		np.power(delta_C_dash/(k_C*S_C),2)+\
		np.power(delta_H_dash/(k_H*S_H),2)+\
		R_T*(delta_C_dash/(k_C*S_C))*(delta_H_dash/(k_H*S_H))\
		)
def main():
	# cie_de_2000([50,2.6772,-79.7751],[50,0,-82.7485],1,1,1)
	# cie_de_2000([50,3.1571,-77.2803],[50,0,-82.7485],1,1,1)
	cie_de_2000([50.0000,2.8361,-74.0200],[50.0000,0.0000,-82.7485],1,1,1)
if __name__ == '__main__':
	main()
