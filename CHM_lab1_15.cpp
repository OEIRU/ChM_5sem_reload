#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <string>
#include <limits>
#include <algorithm> // Для std::max
#include <fstream>   // Для чтения файлов

using namespace std;

// Структура для подсчета операций
struct OperationCounter
{
   long long additions = 0;
   long long multiplications = 0;
   long long divisions = 0;
   long long square_roots = 0;
   long long swaps = 0;

   void reset()
   {
      additions = multiplications = divisions = square_roots = swaps = 0;
   }

   void print() const
   {
      cout << "Количество операций:" << endl;
      cout << "Сложения: " << additions << endl;
      cout << "Умножения: " << multiplications << endl;
      cout << "Деления: " << divisions << endl;
      cout << "Извлечения корней: " << square_roots << endl;
      cout << "Перестановки строк: " << swaps << endl;
   }
};

OperationCounter opCount;

// Структура профильной матрицы
struct ProfileMatrix
{
   int n;                       // Размерность
   vector<int> first_non_zero;  // Индекс первого ненулевого элемента в каждой строке
   vector<vector<double>> rows; // Ненулевые элементы каждой строки

   ProfileMatrix(int size = 0) : n(size), first_non_zero(size, 0), rows(size, vector<double>()) {}
};

// Структура для хранения тестового случая с именем и связанными файлами
struct TestCase
{
   string name;
   string matrixFile;
   string vectorFile;
};

// Преобразование плотной матрицы в профильную
ProfileMatrix denseToProfile(const vector<vector<double>> &dense)
{
   int n = dense.size();
   ProfileMatrix profile(n);
   for (int i = 0; i < n; i++)
   {
      // Найти первый ненулевой элемент в строке
      int first = 0;
      while (first < n && abs(dense[i][first]) < numeric_limits<double>::epsilon())
      {
         first++;
      }
      profile.first_non_zero[i] = first;
      // Сохранить ненулевые элементы
      for (int j = first; j < n; j++)
      {
         profile.rows[i].push_back(dense[i][j]);
      }
   }
   return profile;
}

// Преобразование профильной матрицы обратно в плотную
vector<vector<double>> profileToDense(const ProfileMatrix &profile)
{
   int n = profile.n;
   vector<vector<double>> dense(n, vector<double>(n, 0.0));
   for (int i = 0; i < n; i++)
   {
      int start = profile.first_non_zero[i];
      for (int j = start; j < n; j++)
      {
         dense[i][j] = profile.rows[i][j - start];
      }
   }
   return dense;
}

// Получение элемента из профильной матрицы
double getProfileElement(const ProfileMatrix &profile, int i, int j)
{
   if (j < profile.first_non_zero[i])
      return 0.0;
   return profile.rows[i][j - profile.first_non_zero[i]];
}

// Установка элемента в профильной матрице
void setProfileElement(ProfileMatrix &profile, int i, int j, double value)
{
   if (j < profile.first_non_zero[i])
   {
      // Нужно расширить профиль строки
      int shift = j - profile.first_non_zero[i];
      if (shift > 0)
      {
         profile.rows[i].insert(profile.rows[i].begin(), shift, 0.0);
         profile.first_non_zero[i] = j;
      }
   }
   profile.rows[i][j - profile.first_non_zero[i]] = value;
}

// Функция вывода профильной матрицы в плотном виде
void printMatrix(const ProfileMatrix &profile)
{
   vector<vector<double>> dense = profileToDense(profile);
   int n = dense.size();
   for (int i = 0; i < n; i++)
   {
      for (int j = 0; j < n; j++)
         cout << setw(10) << fixed << setprecision(4) << dense[i][j] << " ";
      cout << endl;
   }
}

// Функция вывода вектора
void printVector(const vector<double> &vec)
{
   for (double v : vec)
      cout << setw(10) << fixed << setprecision(4) << v << " ";
   cout << endl;
}

// Функция умножения матрицы на вектор (плотная матрица)
vector<double> multiplyMatrixVector(const vector<vector<double>> &A, const vector<double> &x)
{
   int n = A.size();
   vector<double> F(n, 0.0);

   for (int i = 0; i < n; i++)
   {
      for (int j = 0; j < n; j++)
      {
         F[i] += A[i][j] * x[j];
         opCount.multiplications += 1;
         opCount.additions += 1;
      }
   }

   return F;
}

// Генерация матрицы Гильберта
vector<vector<double>> generateHilbertMatrix(int n)
{
   vector<vector<double>> H(n, vector<double>(n, 0.0));
   for (int i = 0; i < n; i++)
   {
      for (int j = 0; j < n; j++)
      {
         H[i][j] = 1.0 / (i + j + 1);
      }
   }
   return H;
}

// Функция LU-разложения для профильной матрицы с подсчетом операций
bool LU_SQ_Decomposition(ProfileMatrix &profileA)
{
   int n = profileA.n;
   // Преобразуем профильную матрицу обратно в плотную
   vector<vector<double>> A = profileToDense(profileA);
   vector<vector<double>> L(n, vector<double>(n, 0.0));
   vector<vector<double>> U(n, vector<double>(n, 0.0));

   for (int i = 0; i < n; i++)
   {
      for (int j = 0; j <= i; j++)
      {
         double sum = 0.0;

         if (i == j)
         {
            for (int k = 0; k < j; k++)
            {
               sum += L[j][k] * U[k][j];
               opCount.additions += 1;
               opCount.multiplications += 1;
            }
            double value = A[j][j] - sum;
            // Проверка на отрицательное или нулевое значение перед извлечением корня
            if (value <= 0)
            {
               cout << "Matrix is NOT LU(sq) decomposable!" << endl;
               return false;
            }
            L[j][j] = sqrt(value);
            opCount.square_roots += 1;
            U[j][j] = 1.0; // Диагональ U всегда 1 в Doolittle's методе
         }
         else
         {
            for (int k = 0; k < j; k++)
            {
               sum += L[i][k] * U[k][j];
               opCount.additions += 1;
               opCount.multiplications += 1;
            }
            if (L[j][j] == 0)
            {
               cout << "Деление на ноль при LU-разложении!" << endl;
               return false;
            }
            U[j][i] = (A[j][i] - sum) / L[j][j];
            opCount.additions += 1;
            opCount.multiplications += 1;
            opCount.divisions += 1;
            L[i][j] = U[j][i];
         }
      }
   }

   // Объединяем L и U в одну матрицу: U выше диагонали, L ниже (диагональ L - корни)
   vector<vector<double>> LU(n, vector<double>(n, 0.0));
   for (int i = 0; i < n; i++)
   {
      for (int j = 0; j < n; j++)
      {
         if (i > j)
         {
            LU[i][j] = L[i][j];
         }
         else if (i == j)
         {
            LU[i][j] = L[i][j]; // Хранение диагонали L
         }
         else
         {
            LU[i][j] = U[j][i];
         }
      }
   }

   // Преобразуем обратно в профильную матрицу
   profileA = denseToProfile(LU);

   return true;
}

// Метод Гаусса с выбором ведущего элемента для профильной матрицы с подсчетом операций
bool GaussianEliminationPartialPivoting(ProfileMatrix &profileA, vector<double> &b, vector<double> &solution)
{
   int n = profileA.n;
   // Преобразуем профильную матрицу обратно в плотную
   vector<vector<double>> A = profileToDense(profileA);

   // Создаем копию вектора b для работы
   vector<double> augmented_b = b;

   // Применяем метод Гаусса с выбором ведущего элемента
   for (int i = 0; i < n; i++)
   {
      // Поиск максимального элемента для частичного выбора ведущего
      int maxRow = i;
      double maxElem = abs(A[i][i]);
      for (int k = i + 1; k < n; k++)
      {
         if (abs(A[k][i]) > maxElem)
         {
            maxElem = abs(A[k][i]);
            maxRow = k;
         }
      }

      if (maxElem < numeric_limits<double>::epsilon())
      {
         cout << "Матрица вырождена!" << endl;
         return false;
      }

      // Перестановка строк, если необходимо
      if (maxRow != i)
      {
         swap(A[i], A[maxRow]);
         swap(augmented_b[i], augmented_b[maxRow]); // Перестановка вектора b
         opCount.swaps += 1;
      }

      // Приведение к верхнетреугольному виду
      for (int k = i + 1; k < n; k++)
      {
         double factor = A[k][i] / A[i][i];
         opCount.divisions += 1;
         A[k][i] = 0.0;
         for (int j = i + 1; j < n; j++)
         {
            A[k][j] -= factor * A[i][j];
            opCount.multiplications += 1;
            opCount.additions += 1;
         }
         augmented_b[k] -= factor * augmented_b[i];
         opCount.multiplications += 1;
         opCount.additions += 1;
      }
   }

   // Обратный ход для решения системы
   solution = vector<double>(n, 0.0);
   for (int i = n - 1; i >= 0; i--)
   {
      solution[i] = augmented_b[i];
      for (int j = i + 1; j < n; j++)
      {
         solution[i] -= A[i][j] * solution[j];
         opCount.multiplications += 1;
         opCount.additions += 1;
      }
      solution[i] /= A[i][i];
      opCount.divisions += 1;
   }

   return true;
}

// Функция для генерации всех тестовых случаев
vector<TestCase> generateTestCases()
{
   vector<TestCase> tests;

   // Тест 1: Минимальный размер (1x1)
   tests.push_back(TestCase{
       "Тест 1: Минимальный размер (1x1)",
       "tests/matrix1.txt",
       "tests/vector1.txt"});

   // Тест 2: Маленькая регулярная матрица (2x2)
   tests.push_back(TestCase{
       "Тест 2: Маленькая регулярная матрица (2x2)",
       "tests/matrix2.txt",
       "tests/vector2.txt"});

   // Тест 3: Сингулярная матрица (3x3)
   tests.push_back(TestCase{
       "Тест 3: Сингулярная матрица (3x3)",
       "tests/matrix3.txt",
       "tests/vector3.txt"});

   // Тест 4: Единичная матрица (4x4)
   tests.push_back(TestCase{
       "Тест 4: Единичная матрица (4x4)",
       "tests/matrix4.txt",
       "tests/vector4.txt"});

   // Тест 5: Матрица Гильберта (5x5)
   tests.push_back(TestCase{
       "Тест 5: Матрица Гильберта (5x5)",
       "tests/matrix5.txt",
       "tests/vector5.txt"});

   // Тест 6: Матрица с нулевыми элементами (3x3)
   tests.push_back(TestCase{
       "Тест 6: Матрица с нулевыми элементами (3x3)",
       "tests/matrix6.txt",
       "tests/vector6.txt"});

   // Тест 7: Несимметричная матрица (4x4)
   tests.push_back(TestCase{
       "Тест 7: Несимметричная матрица (4x4)",
       "tests/matrix7.txt",
       "tests/vector7.txt"});

   // Тест 8: Матрица с большинством нулей (5x5)
   tests.push_back(TestCase{
       "Тест 8: Матрица с большинством нулей (5x5)",
       "tests/matrix8.txt",
       "tests/vector8.txt"});

   // Тест 9: Матрица с повторяющимися строками (3x3)
   tests.push_back(TestCase{
       "Тест 9: Матрица с повторяющимися строками (3x3)",
       "tests/matrix9.txt",
       "tests/vector9.txt"});

   // Тест 10: Большая матрица (10x10) с некоторыми нулями
   tests.push_back(TestCase{
       "Тест 10: Большая матрица (10x10) с некоторыми нулями",
       "tests/matrix10.txt",
       "tests/vector10.txt"});

   return tests;
}

// Функция для чтения матрицы из файла
bool readMatrixFromFile(const string &filename, vector<vector<double>> &matrix, int &size)
{
   ifstream matrixStream(filename);
   if (!matrixStream.is_open())
   {
      cerr << "Не удалось открыть файл матрицы: " << filename << endl;
      return false;
   }

   matrixStream >> size;
   if (size <= 0)
   {
      cerr << "Неверный размер матрицы в файле: " << filename << endl;
      return false;
   }

   matrix.assign(size, vector<double>(size, 0.0));

   for (int i = 0; i < size; i++)
   {
      for (int j = 0; j < size; j++)
      {
         if (!(matrixStream >> matrix[i][j]))
         {
            cerr << "Ошибка чтения матрицы из файла: " << filename << endl;
            return false;
         }
      }
   }

   matrixStream.close();
   return true;
}

// Функция для чтения вектора из файла
bool readVectorFromFile(const string &filename, vector<double> &vec, int size)
{
   ifstream vectorStream(filename);
   if (!vectorStream.is_open())
   {
      cerr << "Не удалось открыть файл вектора: " << filename << endl;
      return false;
   }

   vec.assign(size, 0.0);
   for (int i = 0; i < size; i++)
   {
      if (!(vectorStream >> vec[i]))
      {
         cerr << "Ошибка чтения вектора из файла: " << filename << endl;
         return false;
      }
   }

   vectorStream.close();
   return true;
}

// Функция для загрузки теста из файлов
bool loadTest(const TestCase &test, ProfileMatrix &A, vector<double> &b, int &size)
{
   vector<vector<double>> denseMatrix;
   if (!readMatrixFromFile(test.matrixFile, denseMatrix, size))
   {
      return false;
   }

   if (!readVectorFromFile(test.vectorFile, b, size))
   {
      return false;
   }

   A = denseToProfile(denseMatrix);
   return true;
}

// Функция для вывода текущего теста
void displayCurrentTest(const TestCase &test, const ProfileMatrix &A, const vector<double> &b)
{
   cout << "\n--- Текущий Загруженный Тест ---" << endl;
   cout << "Название теста: " << test.name << endl;
   cout << "Матрица A (из файла " << test.matrixFile << "):" << endl;
   printMatrix(A);
   cout << endl
        << "Вектор b (из файла " << test.vectorFile << "):" << endl;
   printVector(b);
   cout << endl;
}

// Основная функция с меню
int main()
{
   setlocale(LC_ALL, "Russian");

   // Генерация всех тестовых случаев
   vector<TestCase> tests = generateTestCases();

   // Переменные для текущего выбранного теста
   ProfileMatrix A;
   vector<double> b, solution, F;

   // Индекс текущего теста (-1 означает, что тест не выбран)
   int currentTestIndex = -1;
   int mainChoice = -1;
   int algorithmChoice = -1;

   while (true)
   {
      // Главное меню
      cout << "================== Главное Меню ==================" << endl;
      cout << "1. Выбрать тест" << endl;
      cout << "2. Выбрать алгоритм для выполнения" << endl;
      cout << "3. Проверка по Гильберту" << endl;
      cout << "4. Вывести текущий тест" << endl;
      cout << "4. Выход" << endl;
      cout << "Ваш выбор: ";
      cin >> mainChoice;

      if (mainChoice == 1)
      {
         // Меню выбора теста
         cout << "\n--- Список Тестов ---" << endl;
         for (int i = 0; i < tests.size(); i++)
         {
            cout << i + 1 << ". " << tests[i].name << endl;
         }
         cout << "Введите номер теста для загрузки: ";
         int selectedTest;
         cin >> selectedTest;

         if (selectedTest < 1 || selectedTest > tests.size())
         {
            cout << "Неверный номер теста. Попробуйте снова.\n"
                 << endl;
            continue;
         }

         // Загрузка выбранного теста
         TestCase currentTest = tests[selectedTest - 1];
         int size = 0;
         bool loaded = loadTest(currentTest, A, b, size);
         if (loaded)
         {
            currentTestIndex = selectedTest - 1;
            cout << "Тест \"" << currentTest.name << "\" успешно загружен.\n"
                 << endl;
         }
         else
         {
            cout << "Не удалось загрузить тест \"" << currentTest.name << "\".\n"
                 << endl;
         }
      }
      else if (mainChoice == 2)
      {
         // Проверка, загружен ли тест
         if (currentTestIndex == -1)
         {
            cout << "Сначала загрузите тест (опция 1).\n"
                 << endl;
            continue;
         }

         // Меню выбора алгоритма
         cout << "\n--- Выбор Алгоритма ---" << endl;
         cout << "1. LU-разложение" << endl;
         cout << "2. Метод Гаусса с выбором ведущего элемента" << endl;
         cout << "Введите номер алгоритма для выполнения: ";
         cin >> algorithmChoice;

         if (algorithmChoice == 1)
         {
            // Применение LU-разложения
            ProfileMatrix LU = A; // Копия матрицы для разложения
            cout << "\nВыполнение LU-разложения..." << endl;
            opCount.reset();
            bool decomposed = LU_SQ_Decomposition(LU);
            if (decomposed)
            {
               cout << "Разложение LU выполнено успешно." << endl;
               cout << "Матрица L и U (вместе в LU):" << endl;
               printMatrix(LU);
               cout << endl;
               opCount.print();

               // Умножение матриц L и U на вектор b
               vector<vector<double>> denseLU = profileToDense(LU);
               vector<double> y = multiplyMatrixVector(denseLU, b);
               cout << "\nВектор y (результат умножения LU * b):" << endl;
               printVector(y);

               F = multiplyMatrixVector(denseLU, y);
               cout << "\nВектор F (результат умножения LU * y):" << endl;
               printVector(F);
            }
            else
            {
               cout << "Разложение LU не удалось.\n"
                    << endl;
            }
            cout << "==============================================\n"
                 << endl;
         }
         else if (algorithmChoice == 2)
         {
            // Применение метода Гаусса
            ProfileMatrix GA = A;  // Копия матрицы для разложения
            vector<double> gb = b; // Копия вектора
            solution.clear();

            cout << "\nВыполнение метода Гаусса..." << endl;
            opCount.reset();
            bool success = GaussianEliminationPartialPivoting(GA, gb, solution);
            if (success)
            {
               cout << "Разложение методом Гаусса выполнено успешно." << endl;
               cout << "Верхнетреугольная матрица A после разложения:" << endl;
               printMatrix(GA);
               cout << endl;
               opCount.print();

               // Вывод решения системы
               cout << "\nРешение системы AX = b:" << endl;
               printVector(solution);
            }
            else
            {
               cout << "Разложение методом Гаусса не удалось.\n"
                    << endl;
            }
            cout << "==============================================\n"
                 << endl;
         }
         else
         {
            cout << "Неверный выбор алгоритма. Попробуйте снова.\n"
                 << endl;
         }
      }
      else if (mainChoice == 3)
      {
         // МЕСТО ДЛЯ ГИЛБЕРТА
         int size_matrix; 
         cout << "Введите размерность матрицы Гильберта" << endl;
         cin >> size_matrix;
         generateHilbertMatrix(size_matrix);
         // подкрутить рандом
         // прогон через LU(sq)
         // счетчики
         cout << "Функция находится в процессе разработки" << endl;

      }
      else if (mainChoice == 4)
      {
         // Вывод текущего теста
         if (currentTestIndex == -1)
         {
            cout << "Нет загруженного теста. Сначала загрузите тест (опция 1).\n"
                 << endl;
            continue;
         }

         TestCase currentTest = tests[currentTestIndex];
         cout << "\n--- Текущий Загруженный Тест ---" << endl;
         displayCurrentTest(currentTest, A, b);
      }
      else if (mainChoice == 4)
      {
         // Выход из программы
         cout << "Выход из программы. До свидания!" << endl;
         break;
      }
      else
      {
         cout << "Неверный выбор. Попробуйте снова.\n"
              << endl;
      }
   }

   return 0;
}