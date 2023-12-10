
#include <iostream>
#include <vector>
#include <iomanip>
#include <assert.h>
#include <cmath>


using namespace std;

void get_matrix(vector<vector<double>> matrix, int size)
{
}
void print(const vector<vector<double>>& matrix, int m, int n)

{
    for (int i = 0;i < m;i++)
    {
        for (int j = 0;j < n;j++)
        {
            cout << setw(10) << matrix[i][j] << ' ';

        }

        cout << endl;
    }
}
void print(const vector<double>& matrix, int m)

{
    for (int i = 0;i < m;i++)
    {
       
            cout << setw(10) << matrix[i] << ' ';

      
    }
}

vector<double> find_x(vector<vector<double>> matrix, vector<vector<double>>& matrix1, int m, int n)
{
    double x3 = matrix1[m - 1][n - 1] / matrix1[m - 1][n - 2];
    double x2 = (matrix1[m - 2][n - 1] - x3 * (matrix1[m - 2][n - 2])) / matrix1[m - 2][n - 3];
    double x1 = (matrix1[m - 3][n - 1] - x2 * matrix1[m - 3][n - 3] - x3 * (matrix1[m - 3][n - 2])) / matrix1[m - 3][n - 4];


  
    vector<double> res_x = { x1,x2,x3 };
    return res_x;
}

vector<double> find_b(vector<vector<double>> &matrix, vector<double> &res, int m, int n)
{
    double x4 = res[0] * matrix[m - 3][n - 4] + res[1] * matrix[m - 3][n - 3] + res[2] * matrix[m - 3][n - 2];
    double x5 = res[0] * matrix[m - 2][n - 4] + res[1] * matrix[m - 2][n - 3] + res[2] * matrix[m - 2][n - 2];
    double x6 = res[0] * matrix[m - 1][n - 4] + res[1] * matrix[m - 1][n - 3] + res[2] * matrix[m - 1][n - 2];
    cout << "B" << endl;
    vector<double> res_b = { x4,x5,x6 };
  
    return res_b;
}

vector<double> find_p(vector<vector<double>> matrix, vector<double> &res_b, int m, int n) {
    double p1 = matrix[m - 3][n - 1] - res_b[0];
    double p2 = matrix[m - 2][n - 1] - res_b[1];
    double p3 = matrix[m - 1][n - 1] - res_b[2];

   
    vector<double> res_p = { p1,p2,p3 };
    return res_p;
}


void copy(vector<vector<double>>& vect1, vector<vector<double>>& vect2)
{
    for (int i = 0; i < vect1.size(); i++)
        vect2.push_back(vect1[i]);
}


vector<vector<double>> Gauss_method(const vector<vector<double>>& matrix, vector<vector<double>>& matrix1, int m, int n, int k)
{
    for (int i = 0;i < 3;i++) {
        assert(!(matrix[i][i] == 0));
    }
    
    if (k <= m)
    {
        double max_coef = 1;
        int index = 0;
        for (int i = k;i < m;i++)
        {
            if (max_coef <= matrix1[i][k])
            {
                index = i;
                max_coef = matrix1[i][k];
            }
        }
        if (index != 0)
        {
            swap(matrix1[k], matrix1[index]);

            print(matrix1, 3, 4);
        }
        cout << endl;
        for (int i = k;i < m;i++)
        {
            double del = matrix1[i][k];
            for (int j = k;j < n;j++)
            {
                matrix1[i][j] = matrix1[i][j] / del;
            }


        }



        for (int i = m - 1;i > k;i--)
            for (int j = k;j < n;j++)
            {
                matrix1[i][j] = matrix1[i][j] - matrix1[k][j];
            }

        print(matrix1, 3, 4);
        cout << endl;
        k++;
        if (k != m)
        {
            Gauss_method(matrix, matrix1, 3, 4, k);

        }
       
        k++;
     
    }
    return matrix1;
}
double find_otn_pogr(vector<double>& r1, vector<double>& r2)
{
    double max_1= abs(r1[0]);

    
        for (int i = 0; i < 3;i++)
        {
            if (max_1 <= abs(r1[i]))
            {
                max_1 = abs(r1[i]);
            }

        }
    
    vector<double> res_v={0,0,0};
    
    for (int i = 0; i < 3;i++)
    {
        res_v[i] = r1[i] - r2[i];
    }
    double max_2 = abs(res_v[0]);
    
        for (int i = 0; i < 3;i++)
        {
            if (max_2 <= abs(res_v[i]))
            
                max_2=abs(res_v[i]);
            
        }
    double res=max_2/max_1;
    cout << res;
    return res;
}
vector<vector<double>> LDLt_factorization(const vector<vector<double>>& A) {
    int n = A.size();
    vector<vector<double>> L(n, vector<double>(n, 0.0));
    vector<vector<double>> D(n, vector<double>(n, 0.0));

for (int i = 0; i < n; i++) {
        D[i][i] = A[i][i];
        for (int k = 0; k < i; k++) {
            D[i][i] -= L[i][k] * L[i][k] * D[k][k];
        }
        if (D[i][i] == 0) {
            cout << "Матрица не может быть разложена на LDL^T." << endl;
            exit(1);
        }

        for (int j = i + 1; j < n; j++) {
            L[j][i] = A[j][i];
            for (int k = 0; k < i; k++) {
                L[j][i] -= L[j][k] * L[i][k] * D[k][k];
            }
            L[j][i] /= D[i][i];
        }
    }

    // Построение L^T
    vector<vector<double>> LT(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            LT[i][j] = L[j][i];
        }
    }

    // Вывод матриц L и D
    cout << "Матрица L:" << endl;
    print(L, n, n);
    cout << "Матрица D:" << endl;
    print(D, n, n);
    cout << "Матрица Lt:" << endl;
    print(LT, n, n);

    

    return L;
}

int main()
{
    setlocale(LC_ALL, "ru");
    vector<vector<double>> matrix
    {
       {21.547,-95.510,-96.121,-49.930},
          {10.223,-91.065, -7.343, -12.465},
        {51.218, 12.264,86.457,60.812}
         
    };
    print(matrix, 3, 4);
    vector<vector<double>> matrix1;
    vector<vector<double>> matrix2;
        copy(matrix, matrix1);
    Gauss_method(matrix,matrix1 ,3, 4, 0);

    vector<double> res_x1=find_x(matrix, matrix1, 3, 4);
    cout << setw(15) << "X1" << endl;
    print(res_x1, 3);
    cout << endl;
    vector<double> res_b= find_b(matrix, res_x1, 3, 4);
    cout << setw(15)<< "B" << endl;
    print(res_b, 3);
    cout << endl;
    vector<double> abs = find_p(matrix, res_b, 3, 4);
    cout <<setw(15) << "Абсолютная погрешность" << endl;
    print(abs, 3);
    cout << endl;
    matrix[0][3] = res_b[0];
    matrix[1][3] = res_b[1];
    matrix[2][3] = res_b[2];
    cout << endl;
    cout << endl;
    print(matrix, 3,4);
    cout << endl;

    copy(matrix, matrix2);
    cout << setw(15) << "Подставим вектор B" << endl;
    Gauss_method(matrix, matrix2, 3, 4, 0);
    cout << setw(15) << "X2" << endl;
    vector<double> res_x2 = find_x(matrix, matrix2, 3, 4);
    print(res_x2, 3);
     
    cout << endl;
    cout <<setw(15) << "Относительная пограшность" << endl;
    find_otn_pogr(res_x1, res_x2);
    cout << endl;
    vector<vector<double>> matrix0 = LDLt_factorization(matrix);
 

  
   

    return 0;
}

