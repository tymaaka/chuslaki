
#include <iostream>
#include <cmath>

using namespace std;

void gauss(double** A, double* B, int N, double* dy) //гаусс
{
  double** matrixA = new double* [N];
  for (int i = 0; i < N; i++)
  {
    matrixA[i] = new double[N + 1];
  }
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N + 1; j++)
    {
      if (j == N)
      {
        matrixA[i][j] = B[i];
      }
      else
      {
        matrixA[i][j] = A[i][j];
      }
    }
  }
  double tmp;
  for (int i = 0; i < N; i++)
  {
    double max = abs(matrixA[i][i]);
    int save = i;
    for (int j = i + 1; j < N; j++)
    {
      if (abs(matrixA[j][i]) > max)
      {
        max = abs(matrixA[i][j]);
        save = j;
      }
    }
    std::swap(matrixA[i], matrixA[save]);

    tmp = matrixA[i][i];
    for (int j = N; j >= i; j--)
    {
      matrixA[i][j] /= tmp;
    }
    for (int j = i + 1; j < N; j++)
    {
      tmp = matrixA[j][i];
      for (int k = N; k >= i; k--)
      {
        matrixA[j][k] -= tmp * matrixA[i][k];
      }
    }
  }
  dy[N - 1] = matrixA[N - 1][N];
  for (int i = N - 2; i >= 0; i--)
  {
    dy[i] = matrixA[i][N];
    for (int j = i + 1; j < N; j++)
    {
      dy[i] -= matrixA[i][j] * dy[j];
    }
  }

  for (int i = 0; i < N; i++)
  {
    delete[] matrixA[i];
  }
  delete[] matrixA;
}

double sussi(double* y, double* COEF, int N, int m, double* x, double SUSSI)
{
  double* sum = new double[m + 1];
  SUSSI = 0;
  for (int i = 0; i < m + 1; i++)
  {
    sum[i] = y[i];
    for (int j = 0; j < m + 1; j++)
    {
      sum[i] -= pow(x[i], j) * COEF[j];
    }
    SUSSI += pow(sum[i], 2);
  }


  delete[] sum;
  return sqrt(SUSSI / (N - m - 1));
}

int main()
{
  int N = 7;
  int m = 1;

  double* y = new double[N];
  double* H = new double[N] { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
  double* p = new double[N] { log(760.0), log(674.8), log(598.0), log(528.9), log(466.6), log(410.6), log(360.2)};


  double* POWERX = new double[2 * m];

  for (int i = 0; i < 2 * m; i++) //инициализация массива
  {
    POWERX[i] = 0;
  }

  for (int k = 1; k <= 2 * m; k++) //суммы х в разных степенях
  {
    for (int i = 0; i < N; i++)
    {
      POWERX[k - 1] += pow(H[i], k); //cум=х1^k+x2^k+...

    }
  }

  double** SUMX = new double* [m + 1]; //матрица коэфф
  for (int i = 0; i < m + 1; i++) {
    SUMX[i] = new double[m + 1];
  }

  int k = 0;

  cout << "Коэфиценты матрицы: " << endl;
  for (int l = 0; l < m + 1; l++) //Сформировать матрицу коэффициентов SUMX размером(m+1)*(m+1) путем выполнения операции присваивания
  {
    for (int j = 0; j < m + 1; j++)
    {
      if (l == 0 && j == 0)
      {
        SUMX[l][j] = N;
      }
      else
      {
        k = l + j - 1;
        SUMX[l][j] = POWERX[k];
      }
      cout << SUMX[l][j] << " ";
    } cout << endl;
  }
  cout << endl;

  cout << "Правая часть системы: " << endl;

  double* PRAW = new double[m + 1]; //правая часть системы

  for (int i = 0; i < m + 1; i++)
  {
    PRAW[i] = 0;
  }

  for (int l = 0; l < m + 1; l++)
  {
    for (int i = 0; i < N; i++) //сумм = y1*x1^l + y2*x2^l + ...
    {
      if (l == 0)
      {
        PRAW[l] += p[i];
      }
      else
      {
        PRAW[l] += p[i] * pow(H[i], (l));
      }
    }
    cout << PRAW[l] << " ";
  }

  double* COEF = new double[m + 1]; //коэфф
  cout << endl;

  gauss(SUMX, PRAW, m + 1, COEF);

  //cout << "коэффициенты:" << endl;

  /*for (int i = 0; i < m + 1; i++) {
    cout << COEF[i] <<" ";
  }*/

  cout << endl;

  double SUSSI = 0;

  SUSSI = sussi(p, COEF, N, m, H, SUSSI);

  cout << "SUSSI = " << SUSSI << endl << endl;

  double* newCOEF = new double[m + 1];

  newCOEF[0] = exp(COEF[0]);
  newCOEF[1] = COEF[1] / log(10);

  cout << "Коэфиценты: " << newCOEF[0] << "; " << newCOEF[1] << endl << endl;


  cout << "Новый y: ";
  for (int i = 0; i < N; i++) {
    y[i] = newCOEF[0] * pow(10, newCOEF[1] * H[i]);
    cout << y[i] << "; ";
  }
  cout << endl;

  delete[] y;
  delete[] newCOEF;
  delete[] COEF;

  for (int i = 0; i < m + 1; i++) {
    delete[] SUMX[i];
  }delete[] SUMX;
  delete[] POWERX;
  delete[] PRAW;
  delete[] H;
  delete[] p;
}
