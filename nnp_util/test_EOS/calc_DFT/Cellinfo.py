# Poscar.py
# class for use poscar

from numpy import *
import CVT, math, mPlane

class Cellinfo:
	def __init__(self, lparam, atoms, cotype, atposi):
		self.lparam = array(lparam)
		self.atoms = atoms
		self.cotype = cotype
		self.atposi = []
		for i in range(len(atposi)):
			type = self.idxToAtom(i+1)
			self.atposi.append([type, array(atposi[i])])

	def printforTest(self):
		print self.lparam
		print self.atoms
		print self.cotype
		print self.atposi

	def idxToAtom(self, idx):
		tot = 0; at_idx = -1;
		for i in range(len(self.atoms)):
			tot = tot + self.atoms[i][1]
			if idx <= tot:
				at_idx = i; break;

		return self.atoms[at_idx][0]

	def getTotnum(self):
		tot = 0
		for item in self.atoms:
			tot = tot + item[1]

		return tot

	def changeCandD(self, type):
		if self.cotype[0] != type[0]:
			if self.cotype[0] in ['c','C','k','K']:
				inv_lparam = linalg.inv(self.lparam)
				for item in self.atposi:
					item[1] = dot(item[1],inv_lparam)
				self.cotype = 'Direct'
			else:
				for item in self.atposi:
       	                        	item[1] = dot(item[1],self.lparam)
				self.cotype = 'Cartesian'

	def latMove(self, lvec):
		return lvec[0]*self.lparam[0] + lvec[1]*self.lparam[1] + lvec[2]*self.lparam[2]

	def checkNeighbor(self, core_posi, cutoff):
		#if self.cotype == 'Direct':
		self.changeCandD('Cartesian')
		p_end = [[0,1],[0,1],[0,1]]

		for a in range(6):
		#	print 'plane search: ' + str(a)
			small_tag = True
			if a%2 == 0:
				extrac = a + 1
			else:
				extrac = a - 1

			while small_tag:
				c_p_end = []
	                        for item in p_end:
	                                c_p_end.append(item[:])
	                        del c_p_end[extrac/2][extrac%2]

				pforp_list = []
				for item in c_p_end[0]:
					for jtem in c_p_end[1]:
						for ktem in c_p_end[2]:
							pforp_list.append(self.latMove([item,jtem,ktem]))

				lplane = mPlane.setPlane(pforp_list)
				#print lplane
				pp_dist = mPlane.distToPoint(core_posi, lplane)
				if pp_dist < cutoff:
					if a%2 == 0:
						p_end[a/2][a%2] = p_end[a/2][a%2] - 1
					else:
						p_end[a/2][a%2] = p_end[a/2][a%2] + 1
				else:
					small_tag = False
			#	print 'plane status: ' + str(p_end)

		nei_list = []		
		for i in range(p_end[0][0], p_end[0][1]):
			for j in range(p_end[1][0], p_end[1][1]):
				for k in range(p_end[2][0], p_end[2][1]):
					for a in range(len(self.atposi)):
						pbc_posi = self.atposi[a][1] + self.latMove([i,j,k])
						tmp_dist = linalg.norm(core_posi - pbc_posi)
						if tmp_dist < cutoff and tmp_dist != 0:
							nei_list.append([a+1, self.atposi[a][0], pbc_posi, tmp_dist])

		for i in range(len(nei_list)-1):
			for j in range(len(nei_list)-1-i):
				if nei_list[j][1] > nei_list[j+1][1]:
					tmp = nei_list[j]
					nei_list[j] = nei_list[j+1]
					nei_list[j+1] = tmp
				elif nei_list[j][1] == nei_list[j+1][1]:
					if nei_list[j][3] > nei_list[j+1][3]:
						tmp = nei_list[j]
                                        	nei_list[j] = nei_list[j+1]
                                        	nei_list[j+1] = tmp
		return nei_list

	# Function for writing other files	
	def writePOS(self, name):
		POS = open(name, 'w')
		POS.write('generated_poscar\n  1.0\n')
		
		for row in self.lparam:
			POS.write('     '+str(row[0])+'   '+str(row[1])+'   '+str(row[2])+'\n')

		str_tmp1 = '  '; str_tmp2 = '  '
		for item in self.atoms:
			str_tmp1 = str_tmp1 + item[0] + '   '
			str_tmp2 = str_tmp2 + str(item[1]) + '   '

		POS.write(str_tmp1 + '\n' + str_tmp2 + '\n')

		POS.write('Selective dynamics\n')
		POS.write(self.cotype + '\n')

		for item in self.atposi:
			POS.write('    '+str(item[1][0])+'   '+str(item[1][1])+'   '+str(item[1][2])+'   T T T\n')

		POS.close()

	def writemsi(self, name):
		print 'not support'

class Poscar(Cellinfo):
	def __init__(self, pos_name):
		POS = open(pos_name, 'r')
		POS.readline()

		lp_coeff = float(POS.readline().replace('\n','').strip())
		self.lparam = []
		for i in [0,1,2]:
			temp = POS.readline().replace('\n','').strip().split()
			self.lparam.append([float(temp[0]), float(temp[1]), float(temp[2])])
		self.lparam = lp_coeff * array(self.lparam)

		atlist = POS.readline().replace('\n','').strip().split()
		atnum = POS.readline().replace('\n','').strip().split()

		self.atoms = []
		for i in range(len(atlist)):
			self.atoms.append([atlist[i], int(atnum[i])])

		temp = POS.readline().replace('\n','').strip()
		if 's' ==  temp[0] or 'S' == temp[0]:
			self.cotype = POS.readline().replace('\n','').strip()
		else:
			self.cotype = temp

		self.atposi = []; totnum = self.getTotnum()
		for i in range(totnum):
			temp = POS.readline().replace('\n','').strip().split()
			self.atposi.append([self.idxToAtom(i+1), array([float(temp[0]), float(temp[1]), float(temp[2])])])

		POS.close()

class Cif(Cellinfo):
	def __init__(self, cif_name):
		k_table = {'1' : '19', '2' : '19', '3' : '15', '4' : '15', '5' : '16', '6' : '15', '7' : '15', '8' : '16', '9' : '16', '10' : '15', '11' : '15', '12' : '16', '13' : '15', '14' : '15', '15' : '16', '16' : '7', '17' : '7', '18' : '7', '19' : '7', '19' : '11', '21' : '11', '22' : '8', '23' : '10', '24' : '10', '25' : '7', '26' : '7', '27' : '7', '28' : '7', '29' : '7', '30' : '7', '31' : '7', '32' : '7', '33' : '7', '34' : '7', '35' : '11', '36' : '11', '37' : '11', '38' : '11', '39' : '11', '40' : '11', '41' : '11', '42' : '9', '43' : '8', '44' : '10', '45' : '10', '46' : '10', '47' : '7', '48' : '7', '49' : '7', '50' : '7', '51' : '7', '52' : '7', '53' : '7', '54' : '7', '55' : '7', '56' : '7', '57' : '7', '58' : '7', '59' : '7', '60' : '7', '61' : '7', '62' : '7', '63' : '11', '64' : '11', '65' : '11', '66' : '11', '67' : '11', '68' : '11', '69' : '9', '70' : '8', '71' : '10', '72' : '10', '73' : '10', '74' : '10', '75' : '4', '76' : '4', '77' : '4', '78' : '4', '79' : '6', '80' : '6', '81' : '4', '82' : '5', '83' : '4', '84' : '4', '85' : '4', '86' : '4', '87' : '6', '88' : '6', '89' : '4', '90' : '4', '91' : '4', '92' : '4', '93' : '4', '94' : '4', '95' : '4', '96' : '4', '97' : '6', '98' : '6', '99' : '4', '100' : '4', '101' : '4', '102' : '4', '103' : '4', '104' : '4', '105' : '4', '106' : '4', '107' : '6', '108' : '6', '109' : '6', '110' : '6', '111' : '4', '112' : '4', '113' : '4', '114' : '4', '115' : '4', '116' : '4', '117' : '4', '118' : '4', '119' : '5', '120' : '5', '121' : '5', '122' : '5', '123' : '4', '124' : '4', '125' : '4', '126' : '4', '127' : '4', '128' : '4', '129' : '4', '130' : '4', '131' : '4', '132' : '4', '133' : '4', '134' : '4', '135' : '4', '136' : '4', '137' : '4', '138' : '4', '139' : '6', '140' : '6', '141' : '6', '142' : '6', '143' : '12', '144' : '12', '145' : '12', '146' : '13', '147' : '12', '148' : '13', '149' : '12', '150' : '12', '151' : '12', '152' : '12', '153' : '12', '154' : '12', '155' : '13', '156' : '12', '157' : '12', '158' : '12', '159' : '12', '160' : '13', '161' : '13', '162' : '12', '163' : '12', '164' : '12', '165' : '12', '166' : '13', '167' : '13', '168' : '12', '169' : '12', '170' : '12', '171' : '12', '172' : '12', '173' : '12', '174' : '12', '175' : '12', '176' : '12', '177' : '12', '178' : '12', '179' : '12', '180' : '12', '181' : '12', '182' : '12', '183' : '12', '184' : '12', '185' : '12', '186' : '12', '187' : '12', '188' : '12', '189' : '12', '190' : '12', '191' : '12', '192' : '12', '193' : '12', '194' : '12', '195' : '1', '196' : '2', '197' : '3', '198' : '1', '199' : '3', '200' : '1', '201' : '1', '202' : '2', '203' : '2', '204' : '3', '205' : '1', '206' : '3', '207' : '1', '208' : '1', '209' : '2', '210' : '2', '211' : '3', '212' : '1', '213' : '1', '214' : '3', '215' : '1', '216' : '2', '217' : '3', '218' : '1', '219' : '2', '220' : '3', '221' : '1', '222' : '1', '223' : '1', '224' : '1', '225' : '2', '226' : '2', '227' : '2', '228' : '2', '229' : '3', '230' : '3'}
		CIF = open(cif_name, 'r')
		#reading part
		for line in CIF:
			if '_cell_length_a' in line:
				try:
					avec = float(line.split()[1].split('(')[0])
				except:
					print 'a vector from cif file is not proper value'
			if '_cell_length_b' in line:
				try:
                                	bvec = float(line.split()[1].split('(')[0])
				except:
                                        print 'b vector from cif file is not proper value'
			if '_cell_length_c' in line:
				try:
                                	cvec = float(line.split()[1].split('(')[0])
				except:
                                        print 'c vector from cif file is not proper value'
			if '_cell_angle_alpha' in line:
                                try:
                                        alpha = int(line.split()[1])
                                except:
                                        print 'angle alpha from cif file is not proper value'
			if '_cell_angle_beta' in line:
                                try:
                                        beta = int(line.split()[1])
                                except:
                                        print 'angle beta from cif file is not proper value'
			if '_cell_angle_gamma' in line:
                                try:
                                        gamma = int(line.split()[1])
                                except:
                                        print 'angle gamma from cif file is not proper value'
			if "_symmetry_Int_Tables_number" in line:
				try:
					symnum = int(line.split()[1])
				except:
					print 'sym number from cif file is not proper value'

		CIF.close()

#class Cif(Cellinfo):
#class Msi(Cellinfo):
