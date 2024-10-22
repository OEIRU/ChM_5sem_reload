#include <iostream>
#include <vector>
#include <iomanip> 
#include <fstream> 
#include <cmath>  
#include <string>
#include <limits>
#include <algorithm> // Для std::max
using namespace std;

// Структура для подсчета операций
struct OperationCounter {
    long long additions = 0;
    long long multiplications = 0;
    long long divisions = 0;
    long long square_roots = 0;
    long long swaps = 0;

    void reset() {
        additions = multiplications = divisions = square_roots = swaps = 0;
    }

    void print() const {
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
struct ProfileMatrix {
    int n; // Размерность
    vector<int> first_non_zero; // Индекс первого ненулевого элемента в каждой строке
    vector<vector<double>> rows; // Ненулевые элементы каждой строки

    ProfileMatrix(int size = 0) : n(size), first_non_zero(size, 0), rows(size, vector<double>()) {}
};

// Преобразование плотной матрицы в профильную
ProfileMatrix denseToProfile(const vector<vector<double>>& dense) {
    int n = dense.size();
    ProfileMatrix profile(n);
    for(int i = 0; i < n; i++) {
        // Найти первый ненулевой элемент в строке
        int first = 0;
        while(first < n && abs(dense[i][first]) < numeric_limits<double>::epsilon()) {
            first++;
        }
        profile.first_non_zero[i] = first;
        // Сохранить ненулевые элементы
        for(int j = first; j < n; j++) {
            profile.rows[i].push_back(dense[i][j]);
        }
    }
    return profile;
}

// Преобразование профильной матрицы обратно в плотную
vector<vector<double>> profileToDense(const ProfileMatrix& profile) {
    int n = profile.n;
    vector<vector<double>> dense(n, vector<double>(n, 0.0));
    for(int i = 0; i < n; i++) {
        int start = profile.first_non_zero[i];
        for(int j = start; j < n; j++) {
            dense[i][j] = profile.rows[i][j - start];
        }
    }
    return dense;
}

// Получение элемента из профильной матрицы
double getProfileElement(const ProfileMatrix& profile, int i, int j) {
    if(j < profile.first_non_zero[i]) return 0.0;
    return profile.rows[i][j - profile.first_non_zero[i]];
}

// Установка элемента в профильной матрице
void setProfileElement(ProfileMatrix& profile, int i, int j, double value) {
    if(j < profile.first_non_zero[i]) {
        // Нужно расширить профиль строки
        int old_nz = profile.rows[i].size();
        int new_first = j;
        int shift = new_first - profile.first_non_zero[i];
        if(shift > 0) {
            profile.rows[i].insert(profile.rows[i].begin(), shift, 0.0);
            profile.first_non_zero[i] = new_first;
        }
    }
    profile.rows[i][j - profile.first_non_zero[i]] = value;
}

// Функция вывода профильной матрицы в плотном виде
void printMatrix(const ProfileMatrix& profile) 
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
void printVector(const vector<double>& vec) 
{
    for (double v : vec) 
        cout << setw(10) << fixed << setprecision(4) << v << " ";
    cout << endl;
}

// Функция записи матрицы в файл в плотном формате
bool writeMatrix(const string& matrixFile, const ProfileMatrix& profileMatrix) {
    ofstream matrixStream(matrixFile);
    if (!matrixStream.is_open()) {
        cerr << "Не удалось открыть файл для записи матрицы." << endl;
        return false;
    }

    vector<vector<double>> A = profileToDense(profileMatrix);
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

// Функция записи вектора в файл
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

// Функция чтения матрицы и вектора из файлов и преобразования матрицы в профильный формат
bool readMatrixAndVector(const string& matrixFile, const string& vectorFile, ProfileMatrix& A, vector<double>& x) {
    ifstream matrixStream(matrixFile);
    ifstream vectorStream(vectorFile);

    if (!matrixStream.is_open() || !vectorStream.is_open()) 
    {
        cerr << "Не удалось открыть файл." << endl;
        return false;
    }

    int n;
    matrixStream >> n; 
    vector<vector<double>> denseA(n, vector<double>(n, 0.0));
    x = vector<double>(n, 0);

    // Считываем матрицу
    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j < n; j++) 
        {
            matrixStream >> denseA[i][j];
        }
    }

    // Считываем вектор
    for (int i = 0; i < n; i++) 
    {
        vectorStream >> x[i];
    }

    matrixStream.close();
    vectorStream.close();

    // Преобразуем в профильное представление
    A = denseToProfile(denseA);

    return true;
}

// Функция умножения матрицы на вектор (плотная матрица)
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

// LU-разложение для профильной матрицы с подсчетом операций
bool LU_SQ_Decomposition(ProfileMatrix& profileA) 
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
                U[j][j] = 1.0; // Диагональ U всегда 1 в Doolittle's метод (можно изменить при необходимости)
            } 
            else 
            {
                for (int k = 0; k < j; k++) 
                {
                    sum += L[i][k] * U[k][j];
                    opCount.additions +=1;
                    opCount.multiplications +=1;
                }
                if (L[j][j] == 0) {
                    cout << "Деление на ноль при LU-разложении!" << endl;
                    return false;
                }
                U[j][i] = (A[j][i] - sum) / L[j][j];
                opCount.additions +=1;
                opCount.multiplications +=1;
                opCount.divisions +=1;
                L[i][j] = U[j][i];
            }
        }
    }

    // Объединяем L и U в одну матрицу: U выше диагонали, L ниже (диагональ L - 1)
    vector<vector<double>> LU(n, vector<double>(n, 0.0));
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            if(i > j) {
                LU[i][j] = L[i][j];
            }
            else if(i == j) {
                LU[i][j] = L[i][j]; // Можно хранить единицы для U[j][j], если используется Doolittle
            }
            else {
                LU[i][j] = U[i][j];
            }
        }
    }

    // Преобразуем обратно в профильную матрицу
    profileA = denseToProfile(LU);

    return true;
}

// Метод Гаусса с выбором ведущего элемента для профильной матрицы с подсчетом операций
bool GaussianEliminationPartialPivoting(ProfileMatrix& profileA, vector<double>& b) {
    int n = profileA.n;
    // Преобразуем профильную матрицу обратно в плотную
    vector<vector<double>> A = profileToDense(profileA);
    
    // Если вектор b не задан, можно инициализировать его нулями или считать
    // В данном случае считать, что система Ax = b решается отдельно

    // Создаем расширенную матрицу [A | b]
    // Для простоты, если b не используется непосредственно, можно пропустить

    // Применяем метод Гаусса с выбором ведущего элемента
    for(int i=0; i<n; i++) {
        // Поиск максимального элемента для частичного выбора ведущего
        int maxRow = i;
        double maxElem = abs(A[i][i]);
        for(int k=i+1; k<n; k++) {
            if(abs(A[k][i]) > maxElem) {
                maxElem = abs(A[k][i]);
                maxRow = k;
            }
        }

        if(maxElem < numeric_limits<double>::epsilon()) {
            cout << "Матрица вырождена!" << endl;
            return false;
        }

        // Перестановка строк, если необходимо
        if(maxRow != i) {
            swap(A[i], A[maxRow]);
            swap(b[i], b[maxRow]); // Если есть вектор b
            opCount.swaps +=1;
        }

        // Приведение к верхнетреугольному виду
        for(int k=i+1; k<n; k++) {
            double factor = A[k][i] / A[i][i];
            opCount.divisions +=1;
            A[k][i] = 0.0;
            for(int j=i+1; j<n; j++) {
                A[k][j] -= factor * A[i][j];
                opCount.multiplications +=1;
                opCount.additions +=1;
            }
        }
    }

    // Установка верхнетреугольной матрицы обратно в профильную
    profileA = denseToProfile(A);

    return true;
}

// Тестирование LU-разложения на матрице Гильберта
void testHilbertLU() {
    cout << "=== Тест разложения LU_SQ на матрице Гильберта ===" << endl;
    int n;
    cout << "Введите размер матрицы Гильберта: ";
    cin >> n;
    vector<vector<double>> denseH = generateHilbertMatrix(n);
    ProfileMatrix H = denseToProfile(denseH);
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

// Тестирование метода Гаусса на профильной матрице
void testGaussianEliminationProfile(const ProfileMatrix& originalProfileA, const vector<double>& original_b) {
    cout << "=== Тест разложения методом Гаусса с выбором ведущего элемента ===" << endl;
    // Создаем копию матрицы и вектора для работы
    ProfileMatrix A = originalProfileA;
    vector<double> b = original_b;

    cout << "Исходная матрица A:" << endl;
    printMatrix(A);
    cout << endl << "Вектор b:" << endl;
    printVector(b);
    cout << endl;

    opCount.reset();
    bool success = GaussianEliminationPartialPivoting(A, b);
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

    ProfileMatrix A;
    vector<double> x, F;

    // Чтение матрицы и вектора из файлов
    if (!readMatrixAndVector("matrix.txt", "vector.txt", A, x)) 
    {
        return 1;
    }

    cout << "Матрица A:" << endl;
    printMatrix(A);
    cout << endl << "Вектор x:" << endl;
    printVector(x);

    // Копия исходной матрицы для повторного использования в тестах
    ProfileMatrix originalA = A;
    vector<double> original_b = x;

    // LU-разложение
    opCount.reset();
    if (LU_SQ_Decomposition(A)) 
    {
        cout << endl << "Матрица L и U (вместе в A):" << endl;
        printMatrix(A);

        // Умножение матрицы на вектор
        vector<double> y = multiplyMatrixVector(profileToDense(A), x);
        cout << endl << "Вектор y (результат умножения A * x):" << endl;
        printVector(y);

        F = multiplyMatrixVector(profileToDense(A), y);
        cout << endl << "Вектор F (результат умножения A * y):" << endl;
        printVector(F);

        opCount.print();

        // Запись результатов в файлы
        writeMatrix("matrix_decomposed.txt", A);
        writeVector("vector_F.txt", F);
    }
    else {
        cout << "Разложение невозможно." << endl;
    }

    cout << "\n==============================================\n" << endl;

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
                testGaussianEliminationProfile(originalA, original_b);
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