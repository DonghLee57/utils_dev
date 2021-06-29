from numpy import *

def setPlane(p_list):
	n1 = array(p_list[0]); n2 = array(p_list[1]); n3 = array(p_list[2])
	v1 = n3-n1; v2 = n2-n1;

	nv = cross(v1,v2)
	d = dot(nv, n3)

	return [nv[0], nv[1], nv[2], d]

def distToPoint(point, coef_list):
	return abs(dot(point, coef_list[0:3]) - coef_list[3]) / linalg.norm(coef_list[0:3])
