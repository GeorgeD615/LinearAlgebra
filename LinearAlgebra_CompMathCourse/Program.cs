using System;
using MyMathematica;

namespace CompMath_Ex2
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

            MyMath.SLAE(new Matrix(A), new Matrix(b), out Vector X);
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
            rank = MyMath.MatrixRank(A, out Matrix Q, out Matrix P);
            if (rank != 3)
                return false;

            A = new Matrix(new double[,] {
                { 0, 4, 10, 1 },
                { 4, 8, 18, 7 },
                { 10, 18, 40, 17 },
                { 1, 7, 17, 3}
            });


            if (MyMath.MatrixRank(A, out Q, out P) != 2)
                return false;
            //A.PrintMatrix();

            A = new Matrix(new double[,] {
                { 1, -1, 2},
                { 2, -2, 4},
                { -1, 1, -2 }
            });

            if (MyMath.MatrixRank(A, out Q, out P) != 1)
                return false;
            //A.PrintMatrix();

            A = new Matrix(new double[,] {
                { 1, 0, 0},
                { 0, 1, 0},
                { 0, 0, 1}
            });

            if (MyMath.MatrixRank(A, out Q, out P) != 3)
                return false;
            //A.PrintMatrix();

            Matrix b = Matrix.GenerateRandomMatrix(A.GetLength(0), 1, -1000, 1000);
            MyMath.SLAE(A, b, out Vector X);
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

            if (MyMath.MatrixRank(new Matrix(A), out Q, out P) != 2)
                return false;
            //A.PrintMatrix();
            MyMath.SLAE(new Matrix(A), new Matrix(b), out X);
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

        static bool Test_SLAE_Seidel()
        {
            Random rand = new Random();
            int n = rand.Next(2, 10);
            Matrix A = Matrix.GenerateRandomMatrix(n, n, -1000, 1000);
            Matrix b = Matrix.GenerateRandomMatrix(n, 1, -1000, 1000);

            MyMath.SeidelMethod(new Matrix(A), new Matrix(b), out Vector X);
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
        static void Main(string[] args)
        {
            Console.WriteLine("----------Test_LUP_Matrix-----------");
            for (int i = 0; i < 10; ++i)
            {
                if (Test_LUP_Matrix())
                {
                    Console.ForegroundColor = ConsoleColor.Green;
                    Console.WriteLine("TEST_PASSED");
                }
                else
                {
                    Console.ForegroundColor = ConsoleColor.Red;
                    Console.WriteLine("TEST_FAILED");
                }
                Console.ForegroundColor = ConsoleColor.White;
            }

            Console.WriteLine("----------Test_SLAE-----------");
            for (int i = 0; i < 10; ++i)
            {
                if (Test_SLAE())
                {
                    Console.ForegroundColor = ConsoleColor.Green;
                    Console.WriteLine("TEST_PASSED");
                }
                else
                {
                    Console.ForegroundColor = ConsoleColor.Red;
                    Console.WriteLine("TEST_FAILED");
                }
                Console.ForegroundColor = ConsoleColor.White;
            }

            Console.WriteLine("----------Test_DeterminantLU-----------");
            for (int i = 0; i < 10; ++i)
            {
                if (Test_DeterminantLU())
                {
                    Console.ForegroundColor = ConsoleColor.Green;
                    Console.WriteLine("TEST_PASSED");
                }
                else
                {
                    Console.ForegroundColor = ConsoleColor.Red;
                    Console.WriteLine("TEST_FAILED");
                }
                Console.ForegroundColor = ConsoleColor.White;
            }

            Console.WriteLine("----------Test_RevertMatrix-----------");
            for (int i = 0; i < 10; ++i)
            {
                if (Test_RevertMatrix())
                {
                    Console.ForegroundColor = ConsoleColor.Green;
                    Console.WriteLine("TEST_PASSED");
                }
                else
                {
                    Console.ForegroundColor = ConsoleColor.Red;
                    Console.WriteLine("TEST_FAILED");
                }
                Console.ForegroundColor = ConsoleColor.White;
            }

            Console.WriteLine("----------Test_Rank-----------");
            if (Test_Rank())
            {
                Console.ForegroundColor = ConsoleColor.Green;
                Console.WriteLine("TEST_PASSED");
            }
            else
            {
                Console.ForegroundColor = ConsoleColor.Red;
                Console.WriteLine("TEST_FAILED");
            }
            Console.ForegroundColor = ConsoleColor.White;

            Console.WriteLine("----------Test_ConditionNumber-----------");
            for (int i = 0; i < 10; ++i)
            {
                Test_ConditionNumber();
            }

            Console.WriteLine("----------Test_QRdecoposition-----------");
            for (int i = 0; i < 10; ++i)
            {
                if (Test_QRdecomposition())
                {
                    Console.ForegroundColor = ConsoleColor.Green;
                    Console.WriteLine("TEST_PASSED");
                }
                else
                {
                    Console.ForegroundColor = ConsoleColor.Red;
                    Console.WriteLine("TEST_FAILED");
                }
                Console.ForegroundColor = ConsoleColor.White;
            }

            Console.WriteLine("----------Test_SLAE_QR-----------");
            for (int i = 0; i < 10; ++i)
            {
                if (Test_SLAE_QR())
                {
                    Console.ForegroundColor = ConsoleColor.Green;
                    Console.WriteLine("TEST_PASSED");
                }
                else
                {
                    Console.ForegroundColor = ConsoleColor.Red;
                    Console.WriteLine("TEST_FAILED");
                }
                Console.ForegroundColor = ConsoleColor.White;
            }

            //Console.WriteLine("----------Test_SLAE_Seidel-----------");
            //for (int i = 0; i < 10; ++i)
            //{
            //    if (Test_SLAE_Seidel())
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


        }

    }
}

