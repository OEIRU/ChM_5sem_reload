#include <iostream>
#include <vector>
#include <iomanip> // Для форматирования вывода
#include <fstream> // Для работы с файлами
#include <cmath>   // Для sqrt
using namespace std;

// Функция для вывода матрицы
void printMatrix(const vector<vector<double>>& matrix) 
{
   int n = matrix.size();
   for (int i = 0; i < n; i++)
   {
      for (int j = 0; j < n; j++) 
         cout << setw(10) << matrix[i][j] << " ";
      cout << endl;
   }
}

// Функция для вывода вектора
void printVector(const vector<double>& vec) 
{
   for (double v : vec) 
      cout << setw(10) << v << " ";
   cout << endl;
}

// Функция для выполнения LU(sq)-разложения
bool LU_SQ_Decomposition(const vector<vector<double>>& A, vector<vector<double>>& L, vector<vector<double>>& U) 
{
   int n = A.size();
   L = vector<vector<double>>(n, vector<double>(n, 0));
   U = vector<vector<double>>(n, vector<double>(n, 0));

   for (int i = 0; i < n; i++) 
   {
      for (int j = 0; j <= i; j++) 
      {
         double sum = 0;

         // Вычисление элементов нижнетреугольной матрицы L
         if (i == j) {  // Диагональный элемент
            for (int k = 0; k < j; k++) 
            {
               sum += L[j][k] * L[j][k];
            }
            double value = A[j][j] - sum;
            if (value <= 0) 
            {
               cout << "Matrix is NOT LU(sq) decomposable!" << endl;
               return false;
            }
            L[j][j] = sqrt(value);
            U[j][j] = L[j][j]; // L и U имеют одинаковые диагональные элементы
         }
         else 
         {  // Вне диагонали
            for (int k = 0; k < j; k++) 
            {
               sum += L[i][k] * L[j][k];
            }
            L[i][j] = (A[i][j] - sum) / L[j][j];
            U[j][i] = L[i][j]; // U является транспонированной копией L
         }
      }
   }

   return true;
}

// Функция для выполнения матрично-векторного умножения F = A * x
vector<double> multiplyMatrixVector(const vector<vector<double>>& A, const vector<double>& x) {
   int n = A.size();
   vector<double> F(n, 0.0);

   for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
         F[i] += A[i][j] * x[j];
      }
   }

   return F;
}

// Функция для чтения матрицы A и вектора x из файла
bool readMatrixAndVector(const string& matrixFile, const string& vectorFile, vector<vector<double>>& A, vector<double>& x) {
   ifstream matrixStream(matrixFile);
   ifstream vectorStream(vectorFile);

   if (!matrixStream.is_open() || !vectorStream.is_open()) 
   {
      cerr << "Не удалось открыть файл." << endl;
      return false;
   }

   int n;
   matrixStream >> n; // Считываем размер матрицы
   A = vector<vector<double>>(n, vector<double>(n, 0));
   x = vector<double>(n, 0);

   // Считывание матрицы A
   for (int i = 0; i < n; i++) 
   {
      for (int j = 0; j < n; j++) 
      {
         matrixStream >> A[i][j];
      }
   }

   // Считывание вектора x
   for (int i = 0; i < n; i++) 
   {
      vectorStream >> x[i];
   }

   matrixStream.close();
   vectorStream.close();
   return true;
}

int main() {
   setlocale(LC_ALL, "Russian");

   vector<vector<double>> A;
   vector<double> x, F;
   vector<vector<double>> L, U;

   // Считывание матрицы и вектора из файлов
   if (!readMatrixAndVector("matrix.txt", "vector.txt", A, x)) 
   {
      return 1;
   }

   // 1. LU(sq)-разложение
   cout << "Матрица A:" << endl;
   printMatrix(A);
   cout << endl << "Вектор x:" << endl;
   printVector(x);

   if (LU_SQ_Decomposition(A, L, U)) 
   {
      cout << endl << "Матрица L:" << endl;
      printMatrix(L);
      cout << endl << "Матрица U:" << endl;
      printMatrix(U);

      // 2. Прямой ход: решение системы Ly = F
      vector<double> y = multiplyMatrixVector(L, x);
      cout << endl << "Вектор y (результат умножения L * x):" << endl;
      printVector(y);

      // 3. Обратный ход: решение системы Ux = y
      F = multiplyMatrixVector(U, y);
      cout << endl << "Вектор F (результат умножения U * y):" << endl;
      printVector(F);
   }
   else {
      cout << "Разложение невозможно." << endl;
   }

   // Дополнительно, умножаем A на x для получения F
   F = multiplyMatrixVector(A, x);
   cout << endl << "Вектор F (результат умножения A * x):" << endl;
   printVector(F);

   return 0;
}
