#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

vector<double> jacobiIteration(const vector<vector<double>>& A, const vector<double>& b, int maxIter, double tol) {
    int counter;
    int n = A.size();
    vector<double> x(n, 0); // 初始化x为0
    vector<double> x_new(n, 0);
    bool converge = false;

    for (counter = 0; counter < maxIter && !converge; ++counter) {
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    sum += A[i][j] * x[j];
                }
            }
            x_new[i] = (b[i] - sum) / A[i][i];
        }

        // 检查收敛
        double norm = 0.0;
        for (int i = 0; i < n; ++i) {
            norm += (x_new[i] - x[i]) * (x_new[i] - x[i]);
        }
        norm = sqrt(norm);

        if (norm < tol) {
            converge = true;
        }

        x = x_new;
    }
    cout <<"迭代次数: "<<counter<<endl;
    return x;
}

int main() {
    // 示例：Ax = b
    vector<vector<double>> A = {
        { 10, -1, 2, 0 },
        { -1, 11, -1, 3 },
        { 2, -1, 10, -1 },
        { 0, 3, -1, 8 }
    };
    vector<double> b = { 6, 25, -11, 15 };

    int maxIter = 100;  // 最大迭代次数
    double tol = 1e-4;  // 收敛阈值

    vector<double> x = jacobiIteration(A, b, maxIter, tol);

    cout << "Solution:" << endl;
    for (int i = 0; i < x.size(); ++i) {
        cout << "x[" << i << "] = " << x[i] << endl;
    }
    return 0;
}
