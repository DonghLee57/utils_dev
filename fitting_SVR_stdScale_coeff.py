# SHAP 라이브러리 설치
#!pip install shap

import numpy as np
import matplotlib.pyplot as plt
from sklearn.svm import SVR
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.inspection import permutation_importance, PartialDependenceDisplay
import shap

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

# 1. Permutation Importance
result = permutation_importance(svr_rbf, X_test, y_test, n_repeats=30, random_state=0)
sorted_idx = result.importances_mean.argsort()

plt.figure(figsize=(10, 6))
plt.boxplot(result.importances[sorted_idx].T, vert=False, labels=np.array(['Feature 1', 'Feature 2', 'Feature 3'])[sorted_idx])
plt.xlabel('Importance')
plt.title('Permutation Importance (SVR)')
plt.show()

# 2. Partial Dependence Plot (PDP)
features = [0, 1, 2]  # 각 특성에 대한 PDP
fig, ax = plt.subplots(figsize=(10, 6))
PartialDependenceDisplay.from_estimator(svr_rbf, X_train, features, ax=ax, grid_resolution=50)
plt.show()

# 3. SHAP (SHapley Additive exPlanations)
explainer = shap.KernelExplainer(svr_rbf.predict, X_train)
shap_values = explainer.shap_values(X_test)

# SHAP summary plot
shap.summary_plot(shap_values, X_test, feature_names=['Feature 1', 'Feature 2', 'Feature 3'])
