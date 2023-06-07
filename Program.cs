using System;
using System.Collections.Generic;
using MyMathematica;
using System.Diagnostics;

namespace LinearAlgebra_CompMathCourse
{
    partial class Program
    {
        static bool Test_LUP_Matrix()
        {
            Random rand = new Random();
            int n = rand.Next(2, 10);
            Matrix A = Matrix.GenerateRandomMatrix(n, n, -100, 100);
            MyMath.LUP_Matrix(A, out Matrix C, out Matrix P, out int sign);
            Matrix L = Matrix.Eye(n);
            for (int i = 1; i < n; ++i)
            {
                for (int j = 0; j < i; ++j)
                {
                    L[i, j] = C[i, j];
                }
            }
            Matrix U = new Matrix(n, n);
            for (int i = 0; i < n; ++i)
            {
                for (int j = i; j < n; ++j)
                {
                    U[i, j] = C[i, j];
                }
            }
            Matrix A_2 = L * U;
            return P * A == A_2;
        }
        static bool Test_SLAE()
        {
            Random rand = new Random();
            int n = rand.Next(2, 10);
            Matrix A = Matrix.GenerateRandomMatrix(n, n, -1000, 1000);
            Matrix b = Matrix.GenerateRandomMatrix(n, 1, -1000, 1000);

            MyMath.SLAE(new Matrix(A), new Matrix(b), out Vector X, out int opCount);
            for (int i = 0; i < n; ++i)
            {
                double leftPart = 0;
                for (int j = 0; j < n; ++j)
                {
                    leftPart += X[j] * A[i, j];
                }
                if (!MyMath.DoubleCompare(leftPart, b[i, 0]))
                {
                    Console.WriteLine($"{leftPart} {b[i, 0]}");
                    return false;
                }
            }
            return true;
        }
        static bool Test_DeterminantLU()
        {
            Random rand = new Random();
            int n = rand.Next(2, 3);
            Matrix A = Matrix.GenerateRandomMatrix(n, n, -1000, 1000);
            double funcRes = MyMath.DeterminantLU(A);
            double ans;
            if (n == 2)
            {
                ans = A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0];
            }
            else
            {
                ans = A[0, 0] * A[1, 1] * A[2, 2] + A[0, 1] * A[1, 2] * A[2, 0]
                    + A[0, 2] * A[1, 0] * A[2, 1] - A[0, 2] * A[1, 1] * A[2, 0]
                    - A[0, 0] * A[1, 2] * A[2, 1] - A[0, 1] * A[1, 0] * A[2, 2];
            }

            if (!MyMath.DoubleCompare(funcRes, ans))
            {
                Console.WriteLine($"{funcRes} {ans}");
                return false;
            }
            return true;
        }
        static bool Test_RevertMatrix()
        {
            Random rand = new Random();
            int n = rand.Next(2, 10);
            Matrix A = Matrix.GenerateRandomMatrix(n, n, -100, 100);
            Matrix A_REV = MyMath.RevertMatrix(A);
            Matrix E = Matrix.Eye(n);
            return (A * A_REV == E) && (A_REV * A == E);
        }
        static bool Test_Rank()
        {
            int rank = 0;
            Matrix A = new Matrix(new double[,] {
                { 2, 1, 3, 0 },
                { 1, 4, 2, 1 },
                { 3, 2, 2, 6 },
                { 8, 5, 7, 12}
            });
            rank = MyMath.MatrixRank(A, out Matrix Q, out Matrix P, out int opCount);
            if (rank != 3)
                return false;

            A = new Matrix(new double[,] {
                { 0, 4, 10, 1 },
                { 4, 8, 18, 7 },
                { 10, 18, 40, 17 },
                { 1, 7, 17, 3}
            });


            if (MyMath.MatrixRank(A, out Q, out P, out opCount) != 2)
                return false;
            //A.PrintMatrix();

            A = new Matrix(new double[,] {
                { 1, -1, 2},
                { 2, -2, 4},
                { -1, 1, -2 }
            });

            if (MyMath.MatrixRank(A, out Q, out P, out opCount) != 1)
                return false;
            //A.PrintMatrix();

            A = new Matrix(new double[,] {
                { 1, 0, 0},
                { 0, 1, 0},
                { 0, 0, 1}
            });

            if (MyMath.MatrixRank(A, out Q, out P, out opCount) != 3)
                return false;
            //A.PrintMatrix();

            Matrix b = Matrix.GenerateRandomMatrix(A.GetLength(0), 1, -1000, 1000);
            MyMath.SLAE(A, b, out Vector X, out opCount);
            for (int i = 0; i < A.GetLength(0); ++i)
            {
                double leftPart = 0;
                for (int j = 0; j < A.GetLength(0); ++j)
                {
                    leftPart += X[j] * A[i, j];
                }
                if (!MyMath.DoubleCompare(leftPart, b[i, 0]))
                {
                    Console.WriteLine($"{leftPart} {b[i, 0]}");
                    return false;
                }
            }

            //Console.WriteLine("Solution: ");
            //X.PrintVector();

            A = new Matrix(new double[,] {
                { 1, 2, 0},
                { 2, 4, 0},
                { 3, 6, 1}
            });
            b = new Matrix(new double[,] { { 3 }, { 6 }, { 10 } });

            if (MyMath.MatrixRank(new Matrix(A), out Q, out P, out opCount) != 2)
                return false;
            //A.PrintMatrix();
            MyMath.SLAE(new Matrix(A), new Matrix(b), out X, out opCount);
            for (int i = 0; i < A.GetLength(0); ++i)
            {
                double leftPart = 0;
                for (int j = 0; j < A.GetLength(0); ++j)
                {
                    leftPart += X[j] * A[i, j];
                }
                if (!MyMath.DoubleCompare(leftPart, b[i, 0]))
                {
                    Console.WriteLine($"{leftPart} {b[i, 0]}");
                    return false;
                }
            }
            //X.PrintVector();
            return true;

        }
        static void Test_ConditionNumber()
        {
            Random rand = new Random();
            int n = rand.Next(2, 10);
            Matrix A = Matrix.GenerateRandomMatrix(n, n, -100, 100);
            Console.WriteLine(MyMath.ConditionNumber(A));
        }
        static bool Test_QRdecomposition()
        {
            Random rand = new Random();
            int n = rand.Next(2, 10);
            Matrix A = Matrix.GenerateRandomMatrix(n, n, -100, 100);
            MyMath.QRdecomposition(A, out Matrix Q, out Matrix R);
            return A == Q * R;
        }
        static bool Test_SLAE_QR()
        {
            Random rand = new Random();
            int n = rand.Next(2, 10);
            Matrix A = Matrix.GenerateRandomMatrix(n, n, -1000, 1000);
            Matrix b = Matrix.GenerateRandomMatrix(n, 1, -1000, 1000);

            MyMath.SLAE_QR(new Matrix(A), new Matrix(b), out Vector X);
            for (int i = 0; i < n; ++i)
            {
                double leftPart = 0;
                for (int j = 0; j < n; ++j)
                {
                    leftPart += X[j] * A[i, j];
                }
                if (!MyMath.DoubleCompare(leftPart, b[i, 0]))
                {
                    Console.WriteLine($"{leftPart} {b[i, 0]}");
                    return false;
                }
            }
            return true;
        }
        static bool Test_SLAE_Seidel_DiagDominMat()
        {
            Random rand = new Random();
            int n = rand.Next(2, 10);
            Matrix A = Matrix.GenerateRandMatWithDiagDomin(n);
            Matrix b = Matrix.GenerateRandomMatrix(n, 1, -1000, 1000);

            MyMath.SLAE_QR(new Matrix(A), new Matrix(b), out Vector ans);
            MyMath.SeidelMethod(new Matrix(A), new Matrix(b), out Vector X, out int iterCount, out double aPrioriEst);
            Console.WriteLine($"iteration: {iterCount} a priori estimate: {aPrioriEst} ");
            for (int i = 0; i < n; ++i)
            {
                if (!MyMath.DoubleCompare(ans[i], X[i]))
                {
                    Console.WriteLine($"{ans[i]} {X[i]}");
                    return false;
                }
            }
            return true;
        }
        static bool Test_SLAE_Seidel_PositDefMat(Matrix A)
        {
            int n = A.GetLength(0);
            Matrix b = Matrix.GenerateRandomMatrix(n, 1, -1000, 1000);

            MyMath.SLAE_QR(new Matrix(A), new Matrix(b), out Vector ans);
            MyMath.SeidelMethod(new Matrix(A), new Matrix(b), out Vector X, out int iterCount, out double aPrioriEst);
            Console.WriteLine($"iteration: {iterCount}   a priori estimate: {aPrioriEst} ");
            for (int i = 0; i < n; ++i)
            {
                if (!MyMath.DoubleCompare(ans[i], X[i]))
                {
                    Console.WriteLine($"{ans[i]} {X[i]}");
                    return false;
                }
            }
            return true;
        }
        static bool Test_SLAE_Jacobi_DiagDominMat()
        {
            Random rand = new Random();
            int n = rand.Next(2, 10);
            Matrix A = Matrix.GenerateRandMatWithDiagDomin(n);
            Matrix b = Matrix.GenerateRandomMatrix(n, 1, -1000, 1000);

            MyMath.SLAE_QR(new Matrix(A), new Matrix(b), out Vector ans);
            MyMath.JacobiMethod(new Matrix(A), new Matrix(b), out Vector X, out int iterCount, out double aPrioriEst);
            Console.WriteLine($"iteration: {iterCount} a priori estimate: {aPrioriEst} ");
            for (int i = 0; i < n; ++i)
            {
                if (!MyMath.DoubleCompare(ans[i], X[i]))
                {
                    Console.WriteLine($"{ans[i]} {X[i]}");
                    return false;
                }
            }
            return true;
        }
        static bool Test_SLAE_Jacobi_PositDefMat(Matrix A)
        {
            int n = A.GetLength(0);
            Matrix b = Matrix.GenerateRandomMatrix(n, 1, -1000, 1000);

            MyMath.SLAE_QR(new Matrix(A), new Matrix(b), out Vector ans);
            MyMath.JacobiMethod(new Matrix(A), new Matrix(b), out Vector X, out int iterCount, out double aPrioriEst);
            Console.WriteLine($"iteration: {iterCount}   a priori estimate: {aPrioriEst} ");


            for (int i = 0; i < n; ++i)
            {
                if (!MyMath.DoubleCompare(ans[i], X[i]))
                {
                    Console.WriteLine($"{ans[i]} {X[i]}");
                    return false;
                }
            }
            return true;
        }

        static public double answer = 10.83954510946909397740794566485262705081;
        static double f(double x)
        {
            return 4 * Math.Cos(0.5 * x) * Math.Exp((-5.0 * x) / 4.0)
                + 2 * Math.Sin(4.5 * x) * Math.Exp(x / 8.0) + 2;
        }
        static double p(double x, double b)
        {
            return  Math.Pow(b - x, -5.0/6.0);
        }
        static double function(double x, double b = 2.2)
        {
            return f(x) * p(x,b);
        }
        static double integ_function(int n = 500, double a = 1.3, double b = 2.2)
        {
            double h = (b - a) / n;
            double res = 0;
            for(int i = 1; i <= n; ++i)
            {
                res += function(a + (i - 0.5)*h,b)*h;
            }
            return res;
        }
        static public double M(double a = 0, double b = 0.9, int k = 0)
        {
            return (Math.Pow(b, (k + 1.0 / 6.0)) - Math.Pow(a, (k + 1.0 / 6.0)))/(k + 1.0/6.0);
        }
        static public double NewtonKots(double a = 1.3, double b = 2.2, int n = 3, double s = 1.3)
        {
            Matrix X = new Matrix(n);
            for(int i = 0; i < n; ++i)
                for(int j = 0; j < n; ++j)
                    X[i, j] = Math.Pow(((a - s) + (b - a) * j / (n - 1)), i);

            Matrix B = new Matrix(n,1);
            for (int i = 0; i < n; ++i)
                B[i,0] = M(a - s, b - s, i);

            MyMath.SLAE_QR(X, B, out Vector A);

            double res = 0;
            for(int i = 0; i < n; ++i)
            {
                res += A[i] * (f(2.2 - X[1,i]));
            }
            return res;

        }
        static public double NewtonKots_modified(int steps, double a = 1.3, double b = 2.2)
        {
            double res = 0;
            for(int i = 0; i < steps; ++i)
            {
                res += NewtonKots(a + i * (b-a)/steps, a + (i + 1) * (b - a) / steps, 3, a);
            }
            return res;
        }
        static public double Etkin(List<double> E, double L = 2)
        {
            return -Math.Log((E[E.Count-1] - E[E.Count - 2]) / (E[E.Count - 2] - E[E.Count - 3])) / Math.Log(L);
        }
        static public List<double> Runge(List<double> E, double h = 0.45, double L = 2)
        {
            double m = Etkin(E);
            return new List<double>(new double[]{
                (E[E.Count-1] - E[E.Count-2])/(1-Math.Pow(L,-m)),
                (E[E.Count-1] - E[E.Count-2])/(Math.Pow(L,m)-1)});
        }
        static public double Richardson(List<double> E, double h = 0.45, double L = 2)
        {
            double m = Etkin(E);
            int n = E.Count;
            Matrix A = new Matrix(n);
            for(int i = 0; i < n; ++i)
            {
                for(int j = 0; j < n-1; ++j)
                {
                    A[i, j] = Math.Pow(h / Math.Pow(L, i),m+j);
                }
                A[i, n - 1] = 1;
            }
            Matrix b = new Matrix(n, 1);
            for(int i = 0; i < n; ++i)
            {
                b[i, 0] = E[i];
            }
            MyMath.SLAE_QR(A, b, out Vector X);
            //X.PrintVector();
            return X[n - 1];
        }
        static void Main(string[] args)
        {
            //Console.WriteLine("----------Test_LUP_Matrix-----------");
            //for (int i = 0; i < 10; ++i)
            //{
            //    if (Test_LUP_Matrix())
            //    {
            //        Console.ForegroundColor = ConsoleColor.Green;
            //        Console.WriteLine("TEST_PASSED");
            //    }
            //    else
            //    {
            //        Console.ForegroundColor = ConsoleColor.Red;
            //        Console.WriteLine("TEST_FAILED");
            //    }
            //    Console.ForegroundColor = ConsoleColor.White;
            //}

            //Console.WriteLine("----------Test_SLAE-----------");
            //for (int i = 0; i < 10; ++i)
            //{
            //    if (Test_SLAE())
            //    {
            //        Console.ForegroundColor = ConsoleColor.Green;
            //        Console.WriteLine("TEST_PASSED");
            //    }
            //    else
            //    {
            //        Console.ForegroundColor = ConsoleColor.Red;
            //        Console.WriteLine("TEST_FAILED");
            //    }
            //    Console.ForegroundColor = ConsoleColor.White;
            //}

            //Console.WriteLine("----------Test_DeterminantLU-----------");
            //for (int i = 0; i < 10; ++i)
            //{
            //    if (Test_DeterminantLU())
            //    {
            //        Console.ForegroundColor = ConsoleColor.Green;
            //        Console.WriteLine("TEST_PASSED");
            //    }
            //    else
            //    {
            //        Console.ForegroundColor = ConsoleColor.Red;
            //        Console.WriteLine("TEST_FAILED");
            //    }
            //    Console.ForegroundColor = ConsoleColor.White;
            //}

            //Console.WriteLine("----------Test_RevertMatrix-----------");
            //for (int i = 0; i < 10; ++i)
            //{
            //    if (Test_RevertMatrix())
            //    {
            //        Console.ForegroundColor = ConsoleColor.Green;
            //        Console.WriteLine("TEST_PASSED");
            //    }
            //    else
            //    {
            //        Console.ForegroundColor = ConsoleColor.Red;
            //        Console.WriteLine("TEST_FAILED");
            //    }
            //    Console.ForegroundColor = ConsoleColor.White;
            //}

            //Console.WriteLine("----------Test_Rank-----------");
            //if (Test_Rank())
            //{
            //    Console.ForegroundColor = ConsoleColor.Green;
            //    Console.WriteLine("TEST_PASSED");
            //}
            //else
            //{
            //    Console.ForegroundColor = ConsoleColor.Red;
            //    Console.WriteLine("TEST_FAILED");
            //}
            //Console.ForegroundColor = ConsoleColor.White;

            //Console.WriteLine("----------Test_ConditionNumber-----------");
            //for (int i = 0; i < 10; ++i)
            //{
            //    Test_ConditionNumber();
            //}

            //Console.WriteLine("----------Test_QRdecoposition-----------");
            //for (int i = 0; i < 10; ++i)
            //{
            //    if (Test_QRdecomposition())
            //    {
            //        Console.ForegroundColor = ConsoleColor.Green;
            //        Console.WriteLine("TEST_PASSED");
            //    }
            //    else
            //    {
            //        Console.ForegroundColor = ConsoleColor.Red;
            //        Console.WriteLine("TEST_FAILED");
            //    }
            //    Console.ForegroundColor = ConsoleColor.White;
            //}

            //Console.WriteLine("----------Test_SLAE_QR-----------");
            //for (int i = 0; i < 10; ++i)
            //{
            //    if (Test_SLAE_QR())
            //    {
            //        Console.ForegroundColor = ConsoleColor.Green;
            //        Console.WriteLine("TEST_PASSED");
            //    }
            //    else
            //    {
            //        Console.ForegroundColor = ConsoleColor.Red;
            //        Console.WriteLine("TEST_FAILED");
            //    }
            //    Console.ForegroundColor = ConsoleColor.White;
            //}

            //Console.WriteLine("----------Test_SLAE_Seidel_DiagDomin-----------");
            //for (int i = 0; i < 10; ++i)
            //{
            //    if (Test_SLAE_Seidel_DiagDominMat())
            //    {
            //        Console.ForegroundColor = ConsoleColor.Green;
            //        Console.WriteLine("TEST_PASSED");
            //    }
            //    else
            //    {
            //        Console.ForegroundColor = ConsoleColor.Red;
            //        Console.WriteLine("TEST_FAILED");
            //    }
            //    Console.ForegroundColor = ConsoleColor.White;
            //}

            //Console.WriteLine("----------Test_SLAE_Seidel_PositDefMat-----------");
            //Matrix A = new Matrix(new double[,] {
            //    {2, -1, 0 },
            //    {-1, 2, -1 },
            //    {0, -1, 2 }
            //});
            //if (Test_SLAE_Seidel_PositDefMat(A))
            //{
            //    Console.ForegroundColor = ConsoleColor.Green;
            //    Console.WriteLine("TEST_PASSED");
            //}
            //else
            //{
            //    Console.ForegroundColor = ConsoleColor.Red;
            //    Console.WriteLine("TEST_FAILED");
            //}
            //Console.ForegroundColor = ConsoleColor.White;

            //A = new Matrix(new double[,] {
            //    {5, 2, -4 },
            //    {2, 1, -2},
            //    {-4, -2, 5}
            //});
            //if (Test_SLAE_Seidel_PositDefMat(A))
            //{
            //    Console.ForegroundColor = ConsoleColor.Green;
            //    Console.WriteLine("TEST_PASSED");
            //}
            //else
            //{
            //    Console.ForegroundColor = ConsoleColor.Red;
            //    Console.WriteLine("TEST_FAILED");
            //}
            //Console.ForegroundColor = ConsoleColor.White;

            //A = new Matrix(new double[,] {
            //    {1, 2, 0},
            //    {2, 5, -1},
            //    {0, -1, 3}
            //});
            //if (Test_SLAE_Seidel_PositDefMat(A))
            //{
            //    Console.ForegroundColor = ConsoleColor.Green;
            //    Console.WriteLine("TEST_PASSED");
            //}
            //else
            //{
            //    Console.ForegroundColor = ConsoleColor.Red;
            //    Console.WriteLine("TEST_FAILED");
            //}
            //Console.ForegroundColor = ConsoleColor.White;

            //Console.WriteLine("----------Test_SLAE_Jacobi_DiagDomin-----------");
            //for (int i = 0; i < 10; ++i)
            //{
            //    if (Test_SLAE_Jacobi_DiagDominMat())
            //    {
            //        Console.ForegroundColor = ConsoleColor.Green;
            //        Console.WriteLine("TEST_PASSED");
            //    }
            //    else
            //    {
            //        Console.ForegroundColor = ConsoleColor.Red;
            //        Console.WriteLine("TEST_FAILED");
            //    }
            //    Console.ForegroundColor = ConsoleColor.White;
            //}

            //Console.WriteLine("----------Test_SLAE_Jacobi_PositDefMat-----------");
            //A = new Matrix(new double[,] {
            //    {2, -1, 0 },
            //    {-1, 2, -1 },
            //    {0, -1, 2 }
            //});
            //if (Test_SLAE_Jacobi_PositDefMat(A))
            //{
            //    Console.ForegroundColor = ConsoleColor.Green;
            //    Console.WriteLine("TEST_PASSED");
            //}
            //else
            //{
            //    Console.ForegroundColor = ConsoleColor.Red;
            //    Console.WriteLine("TEST_FAILED");
            //}
            //Console.ForegroundColor = ConsoleColor.White;

            //A = new Matrix(new double[,] {
            //    {1, 2, 0},
            //    {2, 5, -1},
            //    {0, -1, 3}
            //});
            //if (Test_SLAE_Jacobi_PositDefMat(A))
            //{
            //    Console.ForegroundColor = ConsoleColor.Green;
            //    Console.WriteLine("TEST_PASSED");
            //}
            //else
            //{
            //    Console.ForegroundColor = ConsoleColor.Red;
            //    Console.WriteLine("TEST_FAILED");
            //}
            //Console.ForegroundColor = ConsoleColor.White;
            Console.WriteLine($"middle rectangle formula with {5000} iterations: {integ_function(5000)}");
            Console.WriteLine($"NewtonKots method: \t{NewtonKots_modified(3)}");
            Console.WriteLine($"Correct answer: \t{answer}");

            double a = 1.3;
            double b = 2.2;
            double eps = 1e-6;
            double m;
            Console.WriteLine();
            Console.WriteLine("---------------СКФ НьютонаКотса с точностью 1е-6-----------------");
            List<double> E = new List<double>(new double[] {
                NewtonKots_modified(1), NewtonKots_modified(2)});
            List<double> R = new List<double>(new double[]{
                1,1});
            int step = 1;
            while(Math.Abs(R[step]) > eps)
            {
                ++step;
                E.Add(NewtonKots_modified((int)Math.Pow(2, step)));
                m = Etkin(E);
                R.Add(Richardson(E, b - a) - E[step]);
                //R.Add(Runge(E,b-a,2)[1]);
                Console.WriteLine($"Step: {Math.Pow(2,step)} \tResult: {E[step]} \tError: {R[step]} \tReal error: {answer - E[step]} \tRich error: {answer - Richardson(E, b - a)} \tEtkin: {m}");
            }

            Console.WriteLine();
            Console.WriteLine("---------------СКФ НьютонаКотса с выбором оптимального шага-----------------");
            
            E = new List<double>(new double[] {
                NewtonKots_modified(1), NewtonKots_modified(2), NewtonKots_modified(4)});
            m = Etkin(E);
            R = new List<double>(new double[] {
                1.0, Runge(E, (b-a)/2)[0], Runge(E, (b-a)/2)[1]});
            int optimalStepsNum =
                (int)Math.Ceiling((b - a) / (Math.Pow(eps / Math.Abs(R[R.Count - 1]), 1.0 / m) * 0.95));
                //(int)(((b - a) / (Math.Pow(eps / Math.Abs(R[R.Count - 1]), 1.0 / m) * 0.95)))+1;
            Console.WriteLine($"Optimal steps num: {optimalStepsNum}");
            E = new List<double>(new double[] {
                NewtonKots_modified(optimalStepsNum),
                NewtonKots_modified(optimalStepsNum*2)});
            optimalStepsNum *= 4;
            R = new List<double>(new double[] { 1 });
            step = 2;
            while(Math.Abs(R[R.Count-1]) > eps)
            {
                E.Add(NewtonKots_modified(optimalStepsNum));
                m = Etkin(E);
                R.Add(Richardson(E, b - a) - E[E.Count-1]);
                optimalStepsNum *= 2;
                ++step;
                Console.WriteLine($"Step: {step} \tResult: {E[E.Count-1]} \tError: {R[R.Count - 1]} \tEtkin: {m}");
            }
        }

        
       

    }
}

