def setBasis(x,y,z) :
        x = x - int(x)
        y = y - int(y)
        z = z - int(z)
        if x < 0 :
                x = x+1
        if y < 0 :
                y = y+1
        if z < 0 :
                z = z+1
        if abs(x-1) < 0.0001 :
                x = 0
        if abs(y-1) < 0.0001 :
                y = 0
        if abs(z-1) < 0.0001 :
                z = 0
        return [x,y,z]

def lat_fc(a,b,c) :
	I = [a,b,c]
	L = [[0,0,0],[0,0,0],[0,0,0]]
	M = [[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]]
	for i in range(3) :
		for j in range(3) :
			L[i][j] = M[i][0]*I[0][j]+M[i][1]*I[1][j]+M[i][2]*I[2][j]
	return [L[0],L[1],L[2]]

def atom_fc(final) :	
	M = [[-1.0,1.0,1.0],[1.0,-1.0,1.0],[1.0,1.0,-1.0]]
	n = len(final)
	i = 0
	while i < n :
		xx = final[i][0]*M[0][0]+final[i][1]*M[1][0]+final[i][2]*M[2][0]
		yy = final[i][0]*M[0][1]+final[i][1]*M[1][1]+final[i][2]*M[2][1]
		zz = final[i][0]*M[0][2]+final[i][1]*M[1][2]+final[i][2]*M[2][2]
		[xx,yy,zz] = setBasis(xx,yy,zz)
		if xx >= 0. and xx < 1. and yy >= 0. and yy < 1. and zz >= 0. and zz < 1. :
			check = 0
			for l in range(i) :
				if abs(xx-final[l][0]) < 0.001 and abs(yy-final[l][1]) < 0.001 and abs(zz-final[l][2]) < 0.001 :
					check = 1
			if check == 0 :
				final[i] = [xx,yy,zz,final[i][3],final[i][4]]
				i = i + 1
			else :
				final[i:i+1] = []
				n = n - 1
		else :
			final[i:i+1] = []
			n = n - 1
	return final

def lat_bc(a,b,c) :
        I = [a,b,c]
        L = [[0,0,0],[0,0,0],[0,0,0]]
        M = [[-0.5,0.5,0.5],[0.5,-0.5,0.5],[0.5,0.5,-0.5]]
        for i in range(3) :
                for j in range(3) :
                        L[i][j] = M[i][0]*I[0][j]+M[i][1]*I[1][j]+M[i][2]*I[2][j]
        return [L[0],L[1],L[2]]

def atom_bc(final) :
        M = [[0,1.0,1.0],[1.0,0,1.0],[1.0,1.0,0]]
        n = len(final)
        i = 0
        while i < n :
                xx = final[i][0]*M[0][0]+final[i][1]*M[1][0]+final[i][2]*M[2][0]
                yy = final[i][0]*M[0][1]+final[i][1]*M[1][1]+final[i][2]*M[2][1]
                zz = final[i][0]*M[0][2]+final[i][1]*M[1][2]+final[i][2]*M[2][2]
		[xx,yy,zz] = setBasis(xx,yy,zz)
                if xx >= 0. and xx < 1. and yy >= 0. and yy < 1. and zz >= 0. and zz < 1. :
                        check = 0
                        for l in range(i) :
                                if abs(xx-final[l][0]) < 0.001 and abs(yy-final[l][1]) < 0.001 and abs(zz-final[l][2]) < 0.001 :
                                        check = 1
                        if check == 0 :
                                final[i] = [xx,yy,zz,final[i][3],final[i][4]]
                                i = i + 1
                        else :
                                final[i:i+1] = []
                                n = n - 1
                else :
                        final[i:i+1] = []
                        n = n - 1
        return final

def lat_cc1(a,b,c) :
        I = [a,b,c]
        L = [[0,0,0],[0,0,0],[0,0,0]]
        M = [[0.5,-0.5,0],[0.5,0.5,0],[0,0,1.0]]
        for i in range(3) :
                for j in range(3) :
                        L[i][j] = M[i][0]*I[0][j]+M[i][1]*I[1][j]+M[i][2]*I[2][j]
        return [L[0],L[1],L[2]]

def atom_cc1(final) :
        M = [[1.0,1.0,0],[-1.0,1.0,0],[0,0,1.0]]
        n = len(final)
        i = 0
        while i < n :
                xx = final[i][0]*M[0][0]+final[i][1]*M[1][0]+final[i][2]*M[2][0]
                yy = final[i][0]*M[0][1]+final[i][1]*M[1][1]+final[i][2]*M[2][1]
                zz = final[i][0]*M[0][2]+final[i][1]*M[1][2]+final[i][2]*M[2][2]
		[xx,yy,zz] = setBasis(xx,yy,zz)
                if xx >= 0. and xx < 1. and yy >= 0. and yy < 1. and zz >= 0. and zz < 1. :
                        check = 0
                        for l in range(i) :
                                if abs(xx-final[l][0]) < 0.001 and abs(yy-final[l][1]) < 0.001 and abs(zz-final[l][2]) < 0.001 :
                                        check = 1
                        if check == 0 :
                                final[i] = [xx,yy,zz,final[i][3],final[i][4]]
                                i = i + 1
                        else :
                                final[i:i+1] = []
                                n = n - 1
                else :
                        final[i:i+1] = []
                        n = n - 1
        return final

def lat_cc2(a,b,c) :
        I = [a,b,c]
        L = [[0,0,0],[0,0,0],[0,0,0]]
        M = [[0.5,0.5,0],[-0.5,0.5,0],[0,0,1.0]]
        for i in range(3) :
                for j in range(3) :
                        L[i][j] = M[i][0]*I[0][j]+M[i][1]*I[1][j]+M[i][2]*I[2][j]
        return [L[0],L[1],L[2]]

def atom_cc2(final) :
        M = [[1.0,-1.0,0],[1.0,1.0,0],[0,0,1.0]]
        n = len(final)
        i = 0
        while i < n :
                xx = final[i][0]*M[0][0]+final[i][1]*M[1][0]+final[i][2]*M[2][0]
                yy = final[i][0]*M[0][1]+final[i][1]*M[1][1]+final[i][2]*M[2][1]
                zz = final[i][0]*M[0][2]+final[i][1]*M[1][2]+final[i][2]*M[2][2]
                [xx,yy,zz] = setBasis(xx,yy,zz)
                if xx >= 0. and xx < 1. and yy >= 0. and yy < 1. and zz >= 0. and zz < 1. :
                        check = 0
                        for l in range(i) :
                                if abs(xx-final[l][0]) < 0.001 and abs(yy-final[l][1]) < 0.001 and abs(zz-final[l][2]) < 0.001 :
                                        check = 1
                        if check == 0 :
                                final[i] = [xx,yy,zz,final[i][3],final[i][4]]
                                i = i + 1
                        else :
                                final[i:i+1] = []
                                n = n - 1
                else :
                        final[i:i+1] = []
                        n = n - 1
        return final

def lat_rh(a,b,c) :
        I = [a,b,c]
        L = [[0,0,0],[0,0,0],[0,0,0]]
        M = [[2.0/3,1.0/3,1.0/3],[-1.0/3,1.0/3,1.0/3],[-1.0/3,-2.0/3,1.0/3]]
        for i in range(3) :
                for j in range(3) :
                        L[i][j] = M[i][0]*I[0][j]+M[i][1]*I[1][j]+M[i][2]*I[2][j]
        return [L[0],L[1],L[2]]

def atom_rh(final) :
        M = [[1.0,-1.0,0],[0,1.0,-1.0],[1,1,1]]
        n = len(final)
        i = 0
        while i < n :
                xx = final[i][0]*M[0][0]+final[i][1]*M[1][0]+final[i][2]*M[2][0]
                yy = final[i][0]*M[0][1]+final[i][1]*M[1][1]+final[i][2]*M[2][1]
                zz = final[i][0]*M[0][2]+final[i][1]*M[1][2]+final[i][2]*M[2][2]
                [xx,yy,zz] = setBasis(xx,yy,zz)
                if xx >= 0. and xx < 1. and yy >= 0. and yy < 1. and zz >= 0. and zz < 1. :
                        check = 0
                        for l in range(i) :
                                if abs(xx-final[l][0]) < 0.001 and abs(yy-final[l][1]) < 0.001 and abs(zz-final[l][2]) < 0.001 :
                                        check = 1
                        if check == 0 :
                                final[i] = [xx,yy,zz,final[i][3],final[i][4]]
                                i = i + 1
                        else :
                                final[i:i+1] = []
                                n = n - 1
                else :
                        final[i:i+1] = []
                        n = n - 1
        return final
