from pymatgen.core.structure import Structure 
from pymatgen.analysis.diffraction.xrd import XRDCalculator
import matplotlib.pyplot as plt

def plot_xrd_from_vasp(poscar_file, wavelength='CuKa', two_theta_range=(20,90)):
    """
    VASP 형식의 구조 파일(POSCAR)로부터 XRD 패턴을 계산하고 플로팅합니다.

    Args:
        poscar_file (str): POSCAR 파일의 경로.
        wavelength (str): X-선 파장 소스. 기본값은 'CuKa'입니다[2].
                         사용 가능한 소스는 XRDCalculator.AVAILABLE_RADIATION 에서
                         확인할 수 있습니다.
    """
    # 1. POSCAR 파일로부터 구조(Structure) 객체 생성
    try:
        structure = Structure.from_file(poscar_file)
        print(f"Successfully loaded structure from '{poscar_file}'.")
        print(f"Formula: {structure.composition.reduced_formula}, Space group: {structure.get_space_group_info()}")
    except Exception as e:
        print(f"Error loading structure from file: {e}")
        return

    # 2. XRDCalculator 객체 생성
    # symprec는 대칭성을 찾을 때의 정밀도를 의미합니다.
    xrd_calculator = XRDCalculator(wavelength=wavelength, symprec=0.1)

    # 3. 구조로부터 XRD 패턴 계산
    # get_pattern 메서드는 2-theta, 강도(intensity), hkl 값 등을 포함하는
    # DiffractionPattern 객체를 반환합니다[6][9].
    pattern = xrd_calculator.get_pattern(structure)

    # 4. 계산된 패턴 플로팅 (기본 플롯)
    # show_plot 메서드를 사용하면 빠르게 기본 형태의 XRD 패턴을 확인할 수 있습니다[3].
    print("\nGenerating default pymatgen XRD plot...")
    
    # get_plot 메서드는 matplotlib의 axes 객체를 반환하여 추가적인 수정이 가능합니다.
    fig, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)
    
    # 피크 데이터 가져오기
    two_theta = pattern.x
    intensity = pattern.y
    hkls = pattern.hkls
    d_hkls = pattern.d_hkls

    # xrd_calculator의 get_plot 메서드를 활용하여 plot 생성
    xrd_calculator.get_plot(structure, two_theta_range=two_theta_range, annotate_peaks=True, ax=ax)
    
    # 그래프 제목 및 라벨 설정
    ax.set_title(f"XRD Pattern for {structure.composition.reduced_formula} ({wavelength})")
    ax.set_xlabel(r"2$\theta$ (degrees)")
    ax.set_ylabel("Intensity (a.u.)")
    
    plt.show()

if __name__ == "__main__":
    poscar_filename = "POSCAR"
    plot_xrd_from_vasp(poscar_filename, two_theta_range=(20,40))
