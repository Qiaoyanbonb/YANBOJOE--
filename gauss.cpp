#include <cmath>
#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

double roundToPrecision(double value, int precision) {
    double factor = pow(10.0, precision - ceil(log10(fabs(value))));
    return round(value * factor) / factor;
}

void gaussianElimination(vector<vector<double>>& matrix, vector<double>& b) {
    int n = matrix.size();

    // Forward Elimination
    for (int k = 0; k < n-1; k++) {
        for (int i = k+1; i < n; i++) {
            double factor = roundToPrecision(matrix[i][k] / matrix[k][k], 4);
            for (int j = k; j < n; j++) {
                matrix[i][j] = roundToPrecision(matrix[i][j] - factor * matrix[k][j], 4);
            }
            b[i] = roundToPrecision(b[i] - factor * b[k], 4);
        }
    }

    // Backward Substitution
    vector<double> x(n);
    for (int i = n-1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i+1; j < n; j++) {
            x[i] = roundToPrecision(x[i] - matrix[i][j] * x[j], 4);
        }
        x[i] = roundToPrecision(x[i] / matrix[i][i], 4);
    }

    // Print the solution
    cout << "Solution using Gaussian Elimination:" << endl;
    cout << fixed << setprecision(4); // 设置输出格式为固定的小数点表示，保留4位小数
    for (int i = 0; i < n; i++) {
        cout << "x" << i << " = " << x[i] << endl;
    }
}


void gaussianEliminationWithPivoting(vector<vector<double>>& matrix, vector<double>& b) {
    int n = matrix.size();

    // Forward Elimination with Pivoting
    for (int k = 0; k < n-1; k++) {
        // Pivoting
        int maxIndex = k;
        double maxValue = abs(matrix[k][k]);
        for (int i = k+1; i < n; i++) {
            if (abs(matrix[i][k]) > maxValue) {
                maxValue = abs(matrix[i][k]);
                maxIndex = i;
            }
        }
        swap(matrix[k], matrix[maxIndex]);
        swap(b[k], b[maxIndex]);

        // Elimination
        for (int i = k+1; i < n; i++) {
            double factor = matrix[i][k] / matrix[k][k];
            for (int j = k; j < n; j++) {
                matrix[i][j] -= factor * matrix[k][j];
            }
            b[i] -= factor * b[k];
        }
    }

    // Backward Substitution
    vector<double> x(n);
    for (int i = n-1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i+1; j < n; j++) {
            x[i] -= matrix[i][j] * x[j];
        }
        x[i] /= matrix[i][i];
    }

    // Print the solution
    cout << "Solution using Gaussian Elimination with Partial Pivoting:" << endl;
    for (int i = 0; i < n; i++) {
        cout << "x" << i << " = " << x[i] << endl;
    }
}

int main() {
    // 初始化一个3x3系数矩阵
    vector<vector<double>> originalMatrix = {
            {0.001, 2.000, 3.000},
            {-1.000, 3.712, 4.632},
            {-2.000, 1.072, 5.643}
    };

    // 初始化常数项向量
    vector<double> originalB = {1.000, 2.000, 3.000};

    // 输出矩阵和向量以验证
    cout << fixed << setprecision(4); // 设置输出格式为固定的小数点表示，保留4位小数
    cout << "Matrix:" << endl;
    for (const auto &row : originalMatrix) {
        for (const auto &elem : row) {
            cout << elem << " ";
        }
        cout << endl;
    }

    cout << "Vector b:" << endl;
    for (const auto &elem : originalB) {
        cout << elem << endl;
    }

    // 为每个函数创建矩阵和向量的副本
    vector<vector<double>> matrixForGaussian = originalMatrix;
    vector<double> bForGaussian = originalB;
    vector<vector<double>> matrixForGaussianWithPivoting = originalMatrix;
    vector<double> bForGaussianWithPivoting = originalB;

    // 调用不带主元法的高斯消元函数
    gaussianElimination(matrixForGaussian, bForGaussian);

    // 调用带主元法的高斯消元函数
    gaussianEliminationWithPivoting(matrixForGaussianWithPivoting, bForGaussianWithPivoting);

    return 0;
}