import numpy as np
from ase import Atoms
from ase.io import read, write

# 구조 A와 B 읽기
structure_A = read('structure_A.lammps')
structure_B = read('structure_B.lammps')

# 두 구조의 cell 정보 가져오기
cell_A = structure_A.get_cell()
cell_B = structure_B.get_cell()

# a, b 벡터가 같은지 확인
if not np.allclose(cell_A[:2], cell_B[:2]):
    raise ValueError("A와 B의 a, b 벡터가 다릅니다.")

# 새로운 cell 생성 (z 방향으로 합친 크기)
new_cell = cell_A.copy()
new_cell[2] += cell_B[2]

# 새로운 구조 생성
combined_structure = Atoms()

# 구조 A를 25% 위치로 이동
for atom in structure_A:
    new_position = atom.position + [0, 0, 0.25 * new_cell[2][2] - cell_A[2][2]]
    combined_structure.append(atom.symbol, position=new_position)

# 구조 B를 75% 위치로 이동
for atom in structure_B:
    new_position = atom.position + [0, 0, 0.75 * new_cell[2][2] - cell_B[2][2]]
    combined_structure.append(atom.symbol, position=new_position)

# 새로운 cell 설정
combined_structure.set_cell(new_cell)

# Interface 두께 설정 (사용자 지정)
thk = 5.0  # Angstrom 단위로 설정 (예시 값)

# Interface 영역의 z 위치 계산
interface_z_min_25 = 0.25 * new_cell[2][2] - thk / 2
interface_z_max_25 = 0.25 * new_cell[2][2] + thk / 2
interface_z_min_75 = 0.75 * new_cell[2][2] - thk / 2
interface_z_max_75 = 0.75 * new_cell[2][2] + thk / 2

print(f"25% 위치 Interface 영역의 z 위치:")
print(f"하한: {interface_z_min_25:.3f}")
print(f"상한: {interface_z_max_25:.3f}")

print(f"\n75% 위치 Interface 영역의 z 위치:")
print(f"하한: {interface_z_min_75:.3f}")
print(f"상한: {interface_z_max_75:.3f}")

# 결과 저장
write('combined_structure.xyz', combined_structure)
