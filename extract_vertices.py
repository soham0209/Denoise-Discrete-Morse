import numpy as np
import sys

if __name__ == '__main__':
	dir_name = sys.argv[1]
	delta = sys.argv[2]
	filepath = dir_name + '/' + dir_name + '_regions' + delta + '.txt'
	vertices_path = dir_name + '/' + dir_name + '_vert.txt'
	out_file = dir_name + '/' + dir_name + '_denoised.txt'
	vertices = np.loadtxt(vertices_path)
	faces = np.loadtxt(filepath).astype(int)
	ids = np.unique(faces.flatten()).astype(int)
	vertices = vertices[ids]
	np.savetxt(out_file, vertices, fmt='%.6f')
	print('Denoised point-cloud is saved to', out_file)
