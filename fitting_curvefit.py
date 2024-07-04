import numpy as np
from scipy.optimize import curve_fit, Bounds

# 예제 데이터 생성 (임의의 데이터)
np.random.seed(0)
a_data = np.linspace(0, 10, 100)
b_data = np.linspace(0, 10, 100)
c_data = np.linspace(0, 10, 100)
y_data = 2 * a_data + 3 * b_data**2 - c_data + np.random.normal(0, 1, 100)

# 모델 정의
def model_function(vars, a, b, c):
    x, y, z = vars
    return a * x + b * y**2 - c * z

# 데이터 피팅
p0 = [1, 1, 1]  # 초기 추정 파라미터
bounds = Bounds([0, 0, 0], [10, 10, 10])  # 파라미터 제한 조건 (하한, 상한)

# 피팅 함수
def fit_function(data, a, b, c):
    return model_function((data[0], data[1], data[2]), a, b, c)

# 피팅
params, covariance = curve_fit(fit_function, (a_data, b_data, c_data), y_data, p0=p0, bounds=bounds)

# 결과 출력
print("Fitted parameters:", params)
print("Covariance of parameters:", covariance)

# 피팅 결과 시각화
import matplotlib.pyplot as plt

y_fitted = fit_function((a_data, b_data, c_data), *params)

plt.figure(figsize=(10, 6))
plt.plot(y_data, 'b-', label='Original data')
plt.plot(y_fitted, 'r-', label='Fitted data')
plt.legend()
plt.xlabel('Sample Index')
plt.ylabel('Value')
plt.title('Data Fitting with Multivariable Regression Model')
plt.show()
