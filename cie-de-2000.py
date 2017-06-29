'''
This cie_de_2000 implementation is based on the papers: 
1. 
@article{Sharma2005TheObservations,
    title = {{The CIEDE2000 color-difference formula: Implementation notes, supplementary test data, and mathematical observations}},
    year = {2005},
    journal = {Color Research {\&} Application},
    author = {Sharma, Gaurav and Wu, Wencheng and Dalal, Edul N},
    number = {1},
    month = {2},
    pages = {21--30},
    volume = {30},
    publisher = {Wiley Subscription Services, Inc., A Wiley Company},
    url = {http://dx.doi.org/10.1002/col.20070},
    doi = {10.1002/col.20070},
    issn = {1520-6378},
    keywords = {CIE, CIE94, CIEDE2000, CIELAB, CMC, color-difference metrics}
}

2. 
@article{Luo2001TheCIEDE2000,
    title = {{The development of the CIE 2000 colour-difference formula: CIEDE2000}},
    year = {2001},
    journal = {Color Research {\&} Application},
    author = {Luo, M R and Cui, G and Rigg, B},
    number = {5},
    month = {10},
    pages = {340--350},
    volume = {26},
    publisher = {John Wiley {\&} Sons, Inc.},
    url = {http://dx.doi.org/10.1002/col.1049},
    doi = {10.1002/col.1049},
    issn = {1520-6378},
    keywords = {BFP, CIE, CIE94, CIEDE2000, CIELAB, CMC, LCD, color difference metrics}
}

	Copyright 2017 Muhammed Shameem mailtoshameempk@gmail.com

    This file is part of Fast_seg.

    Fast_seg is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Fast_seg is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Fast_seg.  If not, see <http://www.gnu.org/licenses/>.
    
'''


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
