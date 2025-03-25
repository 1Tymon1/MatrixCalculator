#include <iostream>
#include <string.h>
#include <math.h>
#include <chrono>
#include <fstream>
#include <thread>
#include <mutex>
#include <functional>
using namespace std;

const int c = 2;
const int d = 9;
const int e = 7;
const int f = 3;
const int maxIterations = 250;
mutex mtx;

struct Data {
    int iterations;
    long double residuum;
    int time;
};

class Matrix {
public:
    long double* numbers;
    int width;
    int heigth;

    Matrix(int heigth, int width, double a1) {
        this->width = width;
        this->heigth = heigth;
        int a2 = -1;
        int a3 = -1;
        numbers = new long double[width * heigth];
        for (int i = 0; i < heigth; i++) {
            for (int j = 0; j < width; j++) {
                if (i == j) {
                    numbers[i * width + j] = a1;
                }
                else if (i + 1 == j || i == j + 1) {
                    numbers[i * width + j] = a2;
                }
                else if (i + 2 == j || i == j + 2) {
                    numbers[i * width + j] = a3;
                }
                else {
                    numbers[i * width + j] = 0;
                }
            }
        }
    }

    Matrix(int heigth, int width, string c) {
        this->width = width;
        this->heigth = heigth;
        int a2 = -1;
        int a3 = -1;
        numbers = new long double[width * heigth];
        for (int i = 0; i < heigth; i++) {
            for (int j = 0; j < width; j++) {
                numbers[i * width + j] = sin((i * width + j + 1) * (f + 1));
            }
        }
    }

    Matrix(int heigth, int width) {
        this->width = width;
        this->heigth = heigth;
        int a2 = -1;
        int a3 = -1;
        numbers = new long double[width * heigth];
        for (int i = 0; i < heigth; i++) {
            for (int j = 0; j < width; j++) {
                numbers[i * width + j] = 0;
            }
        }
    }

    Matrix(Matrix&& a) {
        width = a.width;
        heigth = a.heigth;
        numbers = new long double[width * heigth];
        memcpy(numbers, a.numbers, width * heigth * sizeof(long double));
        a.numbers = NULL;
    }

    ~Matrix() {
        if (numbers != NULL)
            delete[](numbers);
    }

    void print() {
        cout << "[";
        for (int i = 0; i < heigth; i++) {
            if (i != 0) cout << " ";
            for (int j = 0; j < width; j++) {
                cout << numbers[i * width + j];
                cout << ", ";
            }
            if (i != heigth - 1) cout << "\n";
        }
        cout << "]\n";
    }

    Matrix operator*(Matrix& mult) {
        if (width != mult.heigth) {
            throw std::runtime_error("matrixes must be of compatible size!");
        }
        Matrix returnal(heigth, mult.width);
        long double sum = 0;
        for (int i = 0; i < returnal.heigth; i++)
        {
            for (int j = 0; j < returnal.width; j++) {
                sum = 0;
                for (int k = 0; k < mult.heigth; k++)
                {
                    sum += mult.numbers[k * mult.width + j] * numbers[i * width + k];
                    //cout << mult.numbers[k * mult.width + j] << " <-x, y-> " << numbers[i * width + k] << endl;
                }
                returnal.numbers[i * returnal.width + j] = sum;
            }
        }
        return returnal;
    }

    void operator-=(Matrix& sub) {
        if (width != sub.width || heigth != sub.heigth) {
            throw std::runtime_error("matrixes must be of compatible size!");
        }
        for (int i = 0; i < heigth; i++)
        {
            for (int j = 0; j < width; j++) {
                numbers[i * width + j] -= sub.numbers[i * width + j];
            }
        }
    }

    long double& operator[](int i) {
        return numbers[i];
    }

    void operator=(Matrix&& a) {
        width = a.width;
        heigth = a.heigth;
        numbers = new long double[width * heigth];
        memcpy(numbers, a.numbers, width * heigth * sizeof(long double));
        a.numbers = NULL;
    }

    void fill(long double filler) {
        for (int i = 0; i < width * heigth; i++) numbers[i] = filler;
    }

    void copy(Matrix& a) {
        for (int i = 0; i < width * heigth; i++) numbers[i] = a[i];
    }

    void eye() {
        for (int i = 0; i < heigth; i++) {
            for (int j = 0; j < width; j++) {
                if (i == j)
                    numbers[i * width + j] = 1;
            }
        }
    }

};

long double residuum(Matrix& a, Matrix& x, Matrix& b)
{
    Matrix c = a * x;
    c -= b;
    long double result = 0;
    for (int i = 0; i < c.heigth; i++)
    {
        result += c[i] * c[i];
    }
    return sqrt(result);
}

Data jacobi(Matrix& a, Matrix& x, Matrix& b, ofstream& file) {
    Matrix fakeX(x.heigth, x.width);
    int k = 0;
    auto time = chrono::system_clock::now();
    long double res;

    for (; k < maxIterations; k++)
    {
        file << k << ",";
        for (int i = 0; i < a.heigth; i++) {
            long double Sum1 = 0;
            long double Sum2 = 0;

            for (int j = 0; j <= i - 1; j++)
            {
                Sum1 += a[i * a.width + j] * x[j];
            }
            for (int j = i + 1; j < x.heigth; j++)
            {
                Sum2 += a[i * a.width + j] * x[j];
            }

            fakeX[i] = (b[i] - Sum1 - Sum2) / a[i * a.width + i];
        }
        long double* temp;
        temp = fakeX.numbers;
        fakeX.numbers = x.numbers;
        x.numbers = temp;
        res = residuum(a, x, b);
        file << res << endl;
        if (res < 1e-9) {
            k++;
            break;
        }
    }
    Data data;
    data.iterations = k ;
    data.residuum = res;
    data.time = static_cast<unsigned int>(chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now() - time).count());
    return data;
}

Data jacobi(Matrix& a, Matrix& x, Matrix& b) {
    Matrix fakeX(x.heigth, x.width);
    int k = 0;
    auto time = chrono::system_clock::now();
    long double res;

    for (; k < maxIterations; k++)
    {
        for (int i = 0; i < a.heigth; i++) {
            long double Sum1 = 0;
            long double Sum2 = 0;

            for (int j = 0; j <= i - 1; j++)
            {
                Sum1 += a[i * a.width + j] * x[j];
            }
            for (int j = i + 1; j < x.heigth; j++)
            {
                Sum2 += a[i * a.width + j] * x[j];
            }

            fakeX[i] = (b[i] - Sum1 - Sum2) / a[i * a.width + i];
        }
        long double* temp;
        temp = fakeX.numbers;
        fakeX.numbers = x.numbers;
        x.numbers = temp;
        res = residuum(a, x, b);
        if (res < 1e-9) {
            k++;
            break;
        }
    }
    Data data;
    data.iterations = k;
    data.residuum = res;
    data.time = static_cast<unsigned int>(chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now() - time).count());
    return data;
}

Data gauss(Matrix& a, Matrix& x, Matrix& b, ofstream& file) {
    Matrix fakeX(x.heigth, x.width);
    int k = 0;
    auto time = chrono::system_clock::now();
    long double res;

    for (; k < maxIterations; k++)
    {
        file << k << ",";
        for (int i = 0; i < a.heigth; i++) {
            long double Sum1 = 0;
            long double Sum2 = 0;

            for (int j = 0; j <= i - 1; j++)
            {
                Sum1 += a[i * a.width + j] * fakeX[j];
            }
            for (int j = i + 1; j < x.heigth; j++)
            {
                Sum2 += a[i * a.width + j] * x[j];
            }

            fakeX[i] = (b[i] - Sum1 - Sum2) / a[i * a.width + i];
        }
        long double* temp;
        temp = fakeX.numbers;
        fakeX.numbers = x.numbers;
        x.numbers = temp;
        res = residuum(a, x, b);
        file << res << endl;
        if (res < 1e-9) {
            k++;
            break;
        }
    }
    Data data;
    data.iterations = k;
    data.residuum = res;
    data.time = static_cast<unsigned int>(chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now() - time).count());
    return data;
}

Data gauss(Matrix& a, Matrix& x, Matrix& b) {
    Matrix fakeX(x.heigth, x.width);
    int k = 0;
    auto time = chrono::system_clock::now();
    long double res;

    for (; k < maxIterations; k++)
    {
        for (int i = 0; i < a.heigth; i++) {
            long double Sum1 = 0;
            long double Sum2 = 0;

            for (int j = 0; j <= i - 1; j++)
            {
                Sum1 += a[i * a.width + j] * fakeX[j];
            }
            for (int j = i + 1; j < x.heigth; j++)
            {
                Sum2 += a[i * a.width + j] * x[j];
            }

            fakeX[i] = (b[i] - Sum1 - Sum2) / a[i * a.width + i];
        }
        long double* temp;
        temp = fakeX.numbers;
        fakeX.numbers = x.numbers;
        x.numbers = temp;
        res = residuum(a, x, b);
        if (res < 1e-9) {
            k++;
            break;
        }
    }
    Data data;
    data.iterations = k;
    data.residuum = res;
    data.time = static_cast<unsigned int>(chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now() - time).count());
    return data;
}

Data LU(Matrix& a, Matrix& x, Matrix& b) {
    Matrix y(x.heigth, x.width);
    y.copy(x);
    int its = 0;
    auto time = chrono::system_clock::now();
    long double res;

    Matrix U(a.heigth, a.width);
    U.copy(a);
    Matrix L(a.heigth, a.width);
    L.eye();

    //L, U subdivision
    for (int i = 1; i < L.width; i++) {
        for (int j = 0; j < i; j++) {
            L[i * L.width + j] = U[i * L.width + j] / U[j * L.width + j];
            for (int k = 0; k < L.width; k++) {
                U[i * U.width + k] = U[i * U.width + k] - L[i * L.width + j] * U[j * U.width + k];
            }
        }
    }
 
    for (int i = 0; i < a.width; i++) {
        long double sum = 0;
        for (int j = 0; j < i; j++) {
            sum += L[i * L.width + j] * y[j];
        }
        y[i] = b[i] - sum;
    }
    for (int i = a.heigth - 1; i >= 0; i--) {
        long double sum = 0;
        for (int j = a.width; j > i; j--) {
            sum += U[i * L.width + j] * x[j];
        }
        x[i] = (y[i] - sum) / U[i * L.width + i];
    }

    res = residuum(a, x, b);
    
    Data data;
    data.iterations = its;
    data.residuum = res;
    data.time = static_cast<unsigned int>(chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now() - time).count());
    return data;
}

void compute(int start, int end,int step, ofstream& dataFile) {
    for (int i = start; i <= end; i += step) {
        int totime = 0;
        Matrix A(i * 500, i * 500, 5 + e);
        Matrix b(i * 500, 1, "");
        Matrix x(i * 500, 1);

        Data data1 = jacobi(A, x, b);
        x.fill(0);

        Data data2= gauss(A, x, b);
        x.fill(0);

        Data data3 = LU(A, x, b);

        mtx.lock();
        dataFile << i * 500 << "," << data1.time << "," << data2.time << "," << data3.time << endl;
        dataFile.flush();
        totime = data1.time + data2.time + data3.time;
        cout << i << " loop took: " << totime << endl;
        mtx.unlock();

    }
}


int main()
{
    ofstream dataFile("Matrix_solving.csv");
    dataFile << "task A" << endl << "Size of A:,a1:,f:" << endl << 900 + c * 10 + d << "," << 5 + e << "," << f << endl << endl << "task B" << endl << "iteration#,residuum,jacobi" << endl;
    dataFile.flush();
    Matrix A(900 + c * 10 + d, 900 + c * 10 + d, 5 + e);
    Matrix b(900 + c * 10 + d, 1, "");
    Matrix x(900 + c * 10 + d, 1);

    Data data = jacobi(A, x, b, dataFile);
    dataFile << "total iterations,end residuum,total time[ms]" << endl << data.iterations << "," << data.residuum << "," << data.time << endl << endl << "iteration#,residuum,gauss" << endl;
    dataFile.flush();
    cout << "iterations: " << data.iterations << " time: " << data.time << " residuum: " << data.residuum << "\n";
    x.fill(0);

    data = gauss(A, x, b, dataFile);
    dataFile << "total iterations,end residuum,total time[ms]" << endl << data.iterations << "," << data.residuum << "," << data.time << endl << endl << "task C" << endl << "iteration#,residuum,jacobi" << endl;
    dataFile.flush();
    cout << "iterations: " << data.iterations << " time: " << data.time << " residuum: " << data.residuum << "\n";
    x.fill(0);

    Matrix A2(900 + c * 10 + d, 900 + c * 10 + d, 3);
    data = jacobi(A2, x, b, dataFile);
    dataFile << "total iterations,end residuum,total time[ms]" << endl << data.iterations << "," << data.residuum << "," << data.time << endl << endl << "iteration#,residuum,gauss" << endl;
    dataFile.flush();
    cout << "iterations: " << data.iterations << " time: " << data.time << " residuum: " << data.residuum << "\n";
    x.fill(0);

    data = gauss(A2, x, b, dataFile);
    dataFile << "total iterations,end residuum,total time[ms]" << endl << data.iterations << "," << data.residuum << "," << data.time << endl << endl << "task D" << endl;
    dataFile.flush();
    cout << "iterations: " << data.iterations << " time: " << data.time << " residuum: " << data.residuum << "\n"; 
    x.fill(0);

    data = LU(A2, x, b);
    dataFile << "total iterations,end residuum,total time[ms]" << endl << data.iterations << "," << data.residuum << "," << data.time << endl << endl << "task E" << endl;
    dataFile.flush();
    cout << "iterations: " << data.iterations << " time: " << data.time << " residuum: " << data.residuum << "\n";

    dataFile << "total time\n" << "size of matrix,jacobi,gauss,LU\n";

    thread t1(std::bind(compute, 1, 10, 5, std::ref(dataFile)));
    thread t2(std::bind(compute, 2, 10, 5, std::ref(dataFile)));
    thread t3(std::bind(compute, 3, 10, 5, std::ref(dataFile)));
    thread t4(std::bind(compute, 4, 10, 5, std::ref(dataFile)));
    thread t5(std::bind(compute, 5, 10, 5, std::ref(dataFile)));

    t1.join();
    t2.join();
    t3.join();
    t4.join();
    t5.join();
    //for (int i = 1; i <= 10; i++) {
    //    int totime = 0;
    //    A = Matrix(i * 500, i * 500, 5 + e);
    //    b = Matrix(i * 500, 1, "");
    //    x = Matrix(i * 500, 1);

    //    dataFile << i * 500 << ",";
    //    data = jacobi(A, x, b);
    //    dataFile <<  data.time << ",";
    //    totime += data.time;
    //    x.fill(0);

    //    data = gauss(A, x, b);
    //    dataFile << data.time << ",";
    //    totime += data.time;
    //    x.fill(0);

    //    //data = LU(A, x, b);
    //    //dataFile << data.time;
    //    //totime += data.time;
    //    dataFile << endl;
    //    dataFile.flush();

    //    cout << i << " loop took: " << totime << endl;
    //}


    dataFile.close();
}

