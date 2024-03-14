#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

vector<double> gaussSeidel(const vector<vector<double>>& A, const vector<double>& b, int maxIterations = 100, double tolerance = 1e-5) {
    int n = A.size();
    int counter = 0;
    vector<double> x(n, 0); // 初始化解向量为0

    for (counter = 0; counter < maxIterations; ++counter) {
        vector<double> oldX = x;

        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    sum += A[i][j] * x[j];
                }
            }
            x[i] = (b[i] - sum) / A[i][i];
        }

        // 检查收敛
        double maxDiff = 0;
        for (int i = 0; i < n; ++i) {
            maxDiff = max(maxDiff, abs(x[i] - oldX[i]));
        }
        if (maxDiff < tolerance) {
            break;
        }
    }
    cout<<"迭代次数为: "<<counter<<endl;
    return x;
}

int main() {
    // 示例：解方程组 Ax = b
    vector<vector<double>> A = {{4, 1, 0}, {1, 4, 1}, {0, 1, 4}};
    vector<double> b = {1, 2, 3};

    vector<double> x = gaussSeidel(A, b);

    // 输出结果
    cout << "Solution: ";
    for (double xi : x) {
        cout << xi << " ";
    }
    cout << endl;

    return 0;
}