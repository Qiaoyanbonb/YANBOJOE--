#include <iostream>
#include <vector>

using namespace std;

// 函数用于执行LU分解
void LU_Decomposition(vector<vector<double>>& A, int n, vector<vector<double>>& L, vector<vector<double>>& U, vector<int>& P) {
    P.resize(n);
    iota(P.begin(), P.end(), 0); // 初始化P为0到n-1的序列

    for (int i = 0; i < n; i++) {
        // 部分主元选取
        double max = abs(A[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (abs(A[k][i]) > max) {
                max = abs(A[k][i]);
                maxRow = k;
            }
        }
        swap(P[i], P[maxRow]);

        for (int k = i; k < n; k++) {
            double sum = 0.0;
            for (int j = 0; j < i; j++) {
                sum += L[i][j] * U[j][k];
            }
            U[i][k] = A[P[i]][k] - sum; // 使用P[i]来考虑行交换
        }

        for (int k = i; k < n; k++) {
            if (i == k) L[i][i] = 1; // 对角线元素
            else {
                double sum = 0.0;
                for (int j = 0; j < i; j++) {
                    sum += L[k][j] * U[j][i];
                }
                L[k][i] = (A[P[k]][i] - sum) / U[i][i]; // 使用P[k]来考虑行交换
            }
        }
    }
}

vector<double> forwardSubstitution(const vector<vector<double>>& L, const vector<double>& b) {
    int n = b.size();
    vector<double> y(n, 0);

    for (int i = 0; i < n; i++) {
        double sum = b[i];
        for (int j = 0; j < i; j++) {
            sum -= L[i][j] * y[j];
        }
        y[i] = sum / L[i][i]; // L[i][i] 是 1，但这里保持公式完整性
    }

    return y;
}

vector<double> backwardSubstitution(const vector<vector<double>>& U, const vector<double>& y) {
    int n = y.size();
    vector<double> x(n, 0);

    for (int i = n - 1; i >= 0; i--) {
        double sum = y[i];
        for (int j = i + 1; j < n; j++) {
            sum -= U[i][j] * x[j];
        }
        x[i] = sum / U[i][i];
    }

    return x;
}

int main() {
    int n = 3; // 矩阵的大小
    vector<vector<double>> A = {
            {2, -1, -2},
            {-4, 6, 3},
            {-4, -2, 8}
    };
    vector<double> b = {-3, 9, -6};

    vector<vector<double>> L(n, vector<double>(n, 0));
    vector<vector<double>> U(n, vector<double>(n, 0));
    vector<int> P; // 置换向量

    // 执行LU分解
    LU_Decomposition(A, n, L, U, P);

    // 根据P调整b的顺序
    vector<double> b_permuted(n);
    for(int i = 0; i < n; ++i) {
        b_permuted[i] = b[P[i]];
    }

    // 使用置换后的b进行计算
    vector<double> y = forwardSubstitution(L, b_permuted);
    vector<double> x = backwardSubstitution(U, y);

    // 打印解向量x
    cout << "Solution x:" << endl;
    for (int i = 0; i < n; i++) {
        cout << x[i] << " ";
    }
    cout << endl;

    return 0;
}
