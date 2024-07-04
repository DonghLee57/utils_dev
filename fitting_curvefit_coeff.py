import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

# 예제 데이터 생성 (임의의 다변수 데이터)
np.random.seed(0)
n_samples = 100
X = np.random.rand(n_samples, 3) * 10  # 3개의 특성을 가진 데이터
y = 2 * X[:, 0] + 3 * X[:, 1]**2 - X[:, 2] + np.random.normal(0, 1, n_samples)

# 데이터 스케일링
scaler_X = StandardScaler()
scaler_y = StandardScaler()
X_scaled = scaler_X.fit_transform(X)
y_scaled = scaler_y.fit_transform(y.reshape(-1, 1)).ravel()

# 학습 데이터와 테스트 데이터 분리
X_train, X_test, y_train, y_test = train_test_split(X_scaled, y_scaled, test_size=0.2, random_state=0)

# 선형 회귀 모델 정의 및 학습
lin_reg = LinearRegression()
lin_reg.fit(X_train, y_train)

# 변수 중요성 (계수) 출력
coefficients = lin_reg.coef_
print("Standardized Coefficients:", coefficients)

# 결과 시각화
plt.figure(figsize=(10, 6))
plt.bar(['Feature 1', 'Feature 2', 'Feature 3'], coefficients)
plt.xlabel('Features')
plt.ylabel('Standardized Coefficient')
plt.title('Feature Importance (Standardized Linear Regression)')
plt.show()
