import numpy as np
import matplotlib.pyplot as plt
from sklearn.svm import SVR
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

# SVR 모델 정의 및 학습
svr_rbf = SVR(kernel='rbf', C=100, gamma=0.1, epsilon=0.1)
svr_rbf.fit(X_train, y_train)

# 예측
y_train_pred_scaled = svr_rbf.predict(X_train)
y_test_pred_scaled = svr_rbf.predict(X_test)

# 스케일링된 예측 값을 원래 스케일로 변환
y_train_pred = scaler_y.inverse_transform(y_train_pred_scaled.reshape(-1, 1)).ravel()
y_test_pred = scaler_y.inverse_transform(y_test_pred_scaled.reshape(-1, 1)).ravel()

# 원래 스케일로 변환된 실제 값
y_train = scaler_y.inverse_transform(y_train.reshape(-1, 1)).ravel()
y_test = scaler_y.inverse_transform(y_test.reshape(-1, 1)).ravel()

# 결과 출력 및 시각화
plt.figure(figsize=(10, 6))
plt.scatter(y_train, y_train_pred, color='blue', label='Train data')
plt.scatter(y_test, y_test_pred, color='red', label='Test data')
plt.plot([min(y), max(y)], [min(y), max(y)], color='black', lw=2, linestyle='--')
plt.xlabel('Actual values')
plt.ylabel('Predicted values')
plt.title('SVR with Multivariable Data')
plt.legend()
plt.show()

# 모델 성능 평가
from sklearn.metrics import mean_squared_error, r2_score

print("Train MSE:", mean_squared_error(y_train, y_train_pred))
print("Test MSE:", mean_squared_error(y_test, y_test_pred))
print("Train R2 Score:", r2_score(y_train, y_train_pred))
print("Test R2 Score:", r2_score(y_test, y_test_pred))
