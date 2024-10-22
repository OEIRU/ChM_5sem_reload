#include <iostream>
#include <vector>
#include <iomanip> 
#include <fstream> 
#include <cmath>  
#include <string>
#include <limits>
using namespace std;

// Структура для подсчета операций
struct OperationCounter {
    long long additions = 0;
    long long multiplications = 0;
    long long divisions = 0;
    long long square_roots = 0;
    long long swaps = 0;

    void reset() {
        additions = multiplications = divisions = square_roots = 0;
    }

    void print() const {
        cout << "Количество операций:" << endl;
        cout << "Сложения: " << additions << endl;
        cout << "Умножения: " << multiplications << endl;
        cout << "Деления: " << divisions << endl;
        cout << "Извлечения корней: " << square_roots << endl;
    }
};

// Глобальный счетчик операций
OperationCounter opCount;

// Функция для печати матрицы
void printMatrix(const vector<vector<double>>& matrix) 
{
   int n = matrix.size();
   for (int i = 0; i < n; i++)
   {
      for (int j = 0; j < n; j++) 
         cout << setw(10) << fixed << setprecision(4) << matrix[i][j] << " ";
      cout << endl;
   }
}

// Функция для печати вектора
void printVector(const vector<double>& vec) 
{
   for (double v : vec) 
      cout << setw(10) << fixed << setprecision(4) << v << " ";
   cout << endl;
}

// Функция для записи матрицы в файл
bool writeMatrix(const string& matrixFile, const vector<vector<double>>& A) {
    ofstream matrixStream(matrixFile);
    if (!matrixStream.is_open()) {
        cerr << "Не удалось открыть файл для записи матрицы." << endl;
        return false;
    }

    int n = A.size();
    matrixStream << n << "\n";
    for (int i = 0; i < n; i++) {
        for(int j=0; j < n; j++) {
            matrixStream << fixed << setprecision(4) << A[i][j] << " ";
        }
        matrixStream << "\n";
    }
    matrixStream.close();
    return true;
}

// Функция для записи вектора в файл
bool writeVector(const string& vectorFile, const vector<double>& vec) {
    ofstream vectorStream(vectorFile);
    if (!vectorStream.is_open()) {
        cerr << "Не удалось открыть файл для записи вектора." << endl;
        return false;
    }

    for(double v : vec) {
        vectorStream << fixed << setprecision(4) << v << "\n";
    }
    vectorStream.close();
    return true;
}

// Функция для чтения матрицы и вектора из файлов
bool readMatrixAndVector(const string& matrixFile, const string& vectorFile, vector<vector<double>>& A, vector<double>& x) {
   ifstream matrixStream(matrixFile);
   ifstream vectorStream(vectorFile);

   if (!matrixStream.is_open() || !vectorStream.is_open()) 
   {
      cerr << "Не удалось открыть файл." << endl;
      return false;
   }

   int n;
   matrixStream >> n; 
   A = vector<vector<double>>(n, vector<double>(n, 0));
   x = vector<double>(n, 0);

   for (int i = 0; i < n; i++) 
   {
      for (int j = 0; j < n; j++) 
      {
         matrixStream >> A[i][j];
      }
   }

   for (int i = 0; i < n; i++) 
   {
      vectorStream >> x[i];
   }

   matrixStream.close();
   vectorStream.close();
   return true;
}

// Функция для преобразования матрицы в профильный тип не написана

// LU_SQ_Decomposition с подсчетом операций
bool LU_SQ_Decomposition(vector<vector<double>>& A) 
{
    int n = A.size();
    vector<vector<double>> L(n, vector<double>(n, 0));
    vector<vector<double>> U(n, vector<double>(n, 0));

    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j <= i; j++) 
        {
            double sum = 0;

            if (i == j)
            {
                for (int k = 0; k < j; k++) 
                {
                    sum += L[j][k] * L[j][k];
                    opCount.additions += 1;
                    opCount.multiplications +=1;
                }
                double value = A[j][j] - sum;
                opCount.additions +=1;
                opCount.multiplications +=1;
                if (value <= 0) 
                {
                    cout << "Matrix is NOT LU(sq) decomposable!" << endl;
                    return false;
                }
                L[j][j] = sqrt(value);
                opCount.square_roots +=1;
                U[j][j] = L[j][j]; 
            } 
            else // Вне диагонали
            {
                for (int k = 0; k < j; k++) 
                {
                    sum += L[i][k] * L[j][k];
                    opCount.additions +=1;
                    opCount.multiplications +=1;
                }
                L[i][j] = (A[i][j] - sum) / L[j][j];
                opCount.additions +=1;
                opCount.multiplications +=1;
                opCount.divisions +=1;
                U[j][i] = L[i][j]; 
            }
        }
    }

    // Копируем результаты обратно в A
    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j < n; j++) 
        {
            if (i < j) 
                A[i][j] = U[i][j]; 
            else 
                A[i][j] = L[i][j]; 
        }
    }

    return true;
}

// Функция умножения матрицы и вектора с подсчетом операций
vector<double> multiplyMatrixVector(const vector<vector<double>>& A, const vector<double>& x) {
   int n = A.size();
   vector<double> F(n, 0.0);

   for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
         F[i] += A[i][j] * x[j];
         opCount.multiplications +=1;
         opCount.additions +=1;
      }
   }

   return F;
}

// Генерация матрицы Гильберта
vector<vector<double>> generateHilbertMatrix(int n) {
    vector<vector<double>> H(n, vector<double>(n, 0.0));
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            H[i][j] = 1.0 / (i + j + 1);
        }
    }
    return H;
}

// Функция для разложения матрицы методом Гаусса с частичным выбором ведущего элемента
bool GaussianEliminationPartialPivoting(vector<vector<double>> A, vector<double>& x) {
    int n = A.size();
    // Создаем расширенную матрицу
    // Здесь предполагается, что система уже задана как A и b, но в исходном коде вектор b не задан
    // Поэтому будем считать, что разложение без определения решения
    // Для полноты добавим разложение без учёта b


    for(int i=0; i<n; i++) {
        // Поиск максимального элемента в столбце
        int maxRow = i;
        double maxElem = abs(A[i][i]);
        for(int k=i+1; k<n; k++) {
            if(abs(A[k][i]) > maxElem) {
                maxElem = abs(A[k][i]);
                maxRow = k;
            }
        }

        // Переместить строку с максимальным элементом на позицию i
        if(maxRow != i) {
            swap(A[i], A[maxRow]);
            // Если у вас есть вектор b, его тоже нужно менять
            opCount.swaps +=1;
        }

        // Проверка на нулевой элемент
        if(abs(A[i][i]) < numeric_limits<double>::epsilon()) {
            cout << "Матрица вырождена!" << endl;
            return false;
        }

        // Приведение к верхнетреугольному виду
        for(int k=i+1; k<n; k++) {
            double factor = A[k][i] / A[i][i];
            opCount.divisions +=1;
            A[k][i] = 0.0; // Обнуляем элемент
            for(int j=i+1; j<n; j++) {
                A[k][j] -= factor * A[i][j];
                opCount.multiplications +=1;
                opCount.additions +=1;
            }
        }
    }

    // В данном примере мы не решаем систему, а просто демонстрируем разложение

    return true;
}

// Функция для выполнения теста разложения LU_SQ на матрице Гильберта
void testHilbertLU() {
    cout << "=== Тест разложения LU_SQ на матрице Гильберта ===" << endl;
    int n;
    cout << "Введите размер матрицы Гильберта: ";
    cin >> n;
    vector<vector<double>> H = generateHilbertMatrix(n);
    vector<double> x(n, 1.0); // Пример вектора

    cout << "Матрица Гильберта H:" << endl;
    printMatrix(H);
    cout << endl << "Вектор x:" << endl;
    printVector(x);

    opCount.reset();
    bool decomposed = LU_SQ_Decomposition(H);
    cout << endl;
    if(decomposed) {
        cout << "Матрица L и U (вместе в H):" << endl;
        printMatrix(H);
        cout << endl;
        opCount.print();
    } else {
        cout << "Разложение не удалось." << endl;
    }
    cout << "==============================================\n" << endl;
}

// Функция для выполнения теста разложения методом Гаусса с выбором ведущего элемента
void testGaussianElimination() {
    cout << "=== Тест разложения методом Гаусса с выбором ведущего ===" << endl;
    int n;
    cout << "Введите размер плотной матрицы: ";
    cin >> n;
    vector<vector<double>> A(n, vector<double>(n, 0.0));
    cout << "Введите элементы матрицы по строкам:" << endl;
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            cin >> A[i][j];
        }
    }
    vector<double> x(n, 0.0); // Не используется в текущем контексте

    cout << "Введенная матрица A:" << endl;
    printMatrix(A);
    cout << endl;

    opCount.reset();
    bool success = GaussianEliminationPartialPivoting(A, x);
    if(success) {
        cout << "Разложение методом Гаусса выполнено успешно." << endl;
        cout << "Верхнетреугольная матрица A после разложения:" << endl;
        printMatrix(A);
        cout << endl;
        opCount.print();
    } else {
        cout << "Разложение не удалось." << endl;
    }
    cout << "==============================================\n" << endl;
}

// Основная функция
int main() {
   setlocale(LC_ALL, "Russian");

   // Первоначальный тест: чтение из файлов
   vector<vector<double>> A;
   vector<double> x, F;

   if (!readMatrixAndVector("matrix.txt", "vector.txt", A, x)) 
   {
      return 1;
   }

   cout << "Матрица A:" << endl;
   printMatrix(A);
   cout << endl << "Вектор x:" << endl;
   printVector(x);

   opCount.reset();
   if (LU_SQ_Decomposition(A)) 
   {
      cout << endl << "Матрица L и U (вместе в A):" << endl;
      printMatrix(A);

      vector<double> y = multiplyMatrixVector(A, x);
      cout << endl << "Вектор y (результат умножения A * x):" << endl;
      printVector(y);

      F = multiplyMatrixVector(A, y);
      cout << endl << "Вектор F (результат умножения A * y):" << endl;
      printVector(F);

      opCount.print();

      // Запись результатов обратно в файлы
      writeMatrix("matrix_decomposed.txt", A);
      writeVector("vector_F.txt", F);
   }
   else {
      cout << "Разложение невозможно." << endl;
   }

   cout << "\n==============================================\n" << endl;

   // Дополнительные тесты
   int testChoice;
   do {
       cout << "Выберите тест для выполнения:" << endl;
       cout << "1. Разложение LU_SQ на матрице Гильберта" << endl;
       cout << "2. Разложение методом Гаусса с выбором ведущего элемента" << endl;
       cout << "0. Выход" << endl;
       cout << "Ваш выбор: ";
       cin >> testChoice;

       switch(testChoice) {
           case 1:
               testHilbertLU();
               break;
           case 2:
               testGaussianElimination();
               break;
           case 0:
               cout << "Выход из программы." << endl;
               break;
           default:
               cout << "Неверный выбор. Попробуйте снова." << endl;
       }
   } while(testChoice != 0);

   return 0;
}
