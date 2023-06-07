using System;

namespace MyMathematica
{
    public class MyMath
    {
        public static double eps = 1e-10;
        public static bool DoubleCompare(double a, double b, double diff = 1e-10)
        {
            if (a == 0 && b == 0) return true;
            return Math.Abs(a - b) <= diff * (1 + Math.Max(Math.Abs(a), Math.Abs(b)));
        }
        public static void LUP_Matrix(Matrix A, out Matrix C, out Matrix P, out int sign)
        {
            if (A.GetLength(0) != A.GetLength(1))
            {
                C = null;
                throw new Exception("Error: Wrong Size in LUPMatrix function");
            }
            int n = A.GetLength(0);
            sign = 1;
            int opCount = 0;
            C = new Matrix(A);
            P = Matrix.Eye(n);

            for (int i = 0; i < n; i++)
            {
                double maxValue = 0;
                int row = -1;
                for (int k = i; k < n; k++)
                {
                    if (Math.Abs(C[k, i]) > maxValue)
                    {
                        maxValue = Math.Abs(C[k, i]);
                        row = k;
                    }
                }
                if (row != i)
                {
                    P.SwapRows(row, i);
                    C.SwapRows(row, i);
                    sign *= -1;
                }
                if (maxValue >= 1e-14)
                {
                    for (int j = i + 1; j < n; j++)
                    {
                        C[j, i] /= C[i, i];
                        for (int k = i + 1; k < n; k++)
                            C[j, k] -= C[j, i] * C[i, k];
                        opCount += 1 + (i + 1 - n) * 2;
                    }
                }
            }
            Console.WriteLine($"LU operations count: {opCount}");
        }
        public static void SLAE(Matrix A, Matrix b, out Vector X, out int opCount)
        {

            if (A.GetLength(0) != A.GetLength(1) || A.GetLength(0) != b.GetLength(0) || b.GetLength(1) != 1) throw new Exception("Error: Wrong Size");
            int n = A.GetLength(0);


            int rank = MatrixRank(A, out Matrix Q, out Matrix P, out opCount);
            //Console.WriteLine(rank);
            //A.PrintMatrix();
            Construct_C2LU(A, out Matrix L, out Matrix U);


            double[] y = new double[n];
            double[] x = new double[n];
            Matrix z = new Matrix(n, 1);
            b = P * b;
            y[0] = b[0, 0];
            for (int i = 1; i < n; ++i)
            {
                y[i] = b[i, 0];
                for (int j = 0; j < i; ++j)
                {
                    y[i] -= L[i, j] * y[j];
                    ++opCount;
                }
            }

            bool check = true;

            for (int i = rank; i < n; ++i)
            {
                if (!MyMath.DoubleCompare(y[i], 0))
                {
                    check = false;
                    break;
                }
            }
            if (check)
            {
                z[0, 0] = y[0];
                for (int i = rank - 1; i >= 0; --i)
                {
                    z[i, 0] = y[i];
                    for (int j = i + 1; j < rank; ++j)
                    {
                        z[i, 0] -= U[i, j] * z[j, 0];
                        ++opCount;
                    }
                    z[i, 0] /= U[i, i];
                    ++opCount;
                }
                z = Q * z;
                for (int i = 0; i < n; ++i)
                {
                    x[i] = z[i, 0];
                    ++opCount;
                }
            }
            else
            {
                //x[rank - 1] = y[rank - 1] / U[rank - 1, rank - 1];
                //for (int i = rank - 2; i >= 0; --i)
                //{
                //    x[i] = y[i];
                //    for (int j = n - 1; j > i; --j)
                //    {
                //        x[i] -= U[i, j] * x[j];
                //    }
                //    x[i] /= U[i, i];
                //}
                throw new Exception("Error");
            }


            X = new Vector(x);
        }
        public static void Construct_C2LU(Matrix C, out Matrix L, out Matrix U)
        {
            if (C.GetLength(0) != C.GetLength(1)) throw new Exception("Error: Wrong size");
            int n = C.GetLength(0);
            L = Matrix.Eye(n);
            for (int i = 1; i < n; ++i)
            {
                for (int j = 0; j < i; ++j)
                {
                    L[i, j] = C[i, j];
                }
            }
            U = new Matrix(n, n);
            for (int i = 0; i < n; ++i)
            {
                for (int j = i; j < n; ++j)
                {
                    U[i, j] = C[i, j];
                }
            }
        }
        public static double DeterminantLU(Matrix A)
        {
            LUP_Matrix(A, out Matrix C, out Matrix P, out int sign);
            Construct_C2LU(C, out Matrix L, out Matrix U);
            double res = sign;
            for (int i = 0; i < U.GetLength(0); ++i)
            {
                res *= U[i, i];
            }
            return res;
        }
        public static Matrix RevertMatrix(Matrix A)
        {
            if (A.GetLength(0) != A.GetLength(1)) throw new Exception("Error: Wrong size");
            int n = A.GetLength(0);
            Matrix reverted = new Matrix(n, n);
            Matrix b = new Matrix(n, 1);
            b[0, 0] = 1;
            for (int i = 1; i <= n; ++i)
            {
                SLAE_QR(A, b, out Vector X);
                for (int j = 0; j < n; ++j)
                {
                    reverted[j, i - 1] = X[j];
                }
                if (i == n) break;
                b[i - 1, 0] = 0;
                b[i, 0] = 1;
            }

            return reverted;

        }
        public static double MatrixNorm(Matrix A)
        {
            double res = 0;
            for (int i = 0; i < A.GetLength(0); ++i)
            {
                double sum = 0;
                for (int j = 0; j < A.GetLength(1); ++j)
                {
                    sum += Math.Abs(A[i, j]);
                }
                if (sum > res)
                    res = sum;
            }
            return res;
        }
        public static double ConditionNumber(Matrix A)
        {
            return MatrixNorm(A) * MatrixNorm(RevertMatrix(A));
        }
        public static int MatrixRank(Matrix A, out Matrix Q, out Matrix P, out int opCount)
        {
            int rank = A.GetLength(1);
            P = Matrix.Eye(A.GetLength(0));
            Q = Matrix.Eye(A.GetLength(0));
            opCount = 0;
            for (int row = 0; row < rank; row++)
            {
                double maxVal = Math.Abs(A[row, row]);
                int maxInd = row;
                for (int i = row + 1; i < A.GetLength(0); i++)
                {
                    if (!MyMath.DoubleCompare(A[i, row], 0))
                    {
                        if (Math.Abs(maxVal) < Math.Abs(A[i, row]))
                        {
                            maxVal = A[i, row];
                            maxInd = i;
                        }
                    }
                }

                if (!MyMath.DoubleCompare(maxVal, 0))
                {
                    A.SwapRows(maxInd, row);
                    P.SwapRows(maxInd, row);
                    for (int j = row + 1; j < A.GetLength(0); j++)
                    {
                        A[j, row] /= A[row, row];
                        for (int k = row + 1; k < rank; k++)
                            A[j, k] -= A[j, row] * A[row, k];
                        opCount += 1 + 2 * (rank - row - 1);
                    }
                }
                else
                {
                    rank--;
                    A.SwapColumns(row, rank);
                    Q.SwapColumns(row, rank);
                    row--;
                }
            }

            return rank;

        }
        public static Matrix TransposeMatrix(Matrix A)
        {
            int n = A.GetLength(0);
            int m = A.GetLength(1);
            Matrix res = new Matrix(m, n);
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    res[i, j] = A[j, i];
                }
            }
            return res;
        }
        public static double VectorNorm(Vector A)
        {
            double res = 0;
            for (int i = 0; i < A.Lenght; ++i)
            {
                res += Math.Pow(A[i], 2);
            }
            res = Math.Sqrt(res);
            return res;
        }
        public static void QRdecomposition(Matrix A, out Matrix Q, out Matrix R)
        {
            int n = A.GetLength(0);
            Q = new Matrix(n);
            R = new Matrix(n);
            Vector currQCol = new Vector(n);
            Vector a = new Vector(n);
            Vector prevQCol = new Vector(n);

            for (int j = 0; j < n; ++j)
            {
                a[j] = A[j, 0];
            }
            currQCol = a / VectorNorm(a);
            for (int l = 0; l < n; ++l)
            {
                Q[l, 0] = currQCol[l];
            }

            for (int i = 1; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    currQCol[j] = A[j, i];
                }
                for (int k = 0; k < i; ++k)
                {
                    for (int l = 0; l < n; ++l)
                    {
                        prevQCol[l] = Q[l, k];
                    }
                    currQCol = currQCol - (currQCol * prevQCol) * prevQCol;
                }
                currQCol = currQCol / VectorNorm(currQCol);
                for (int j = 0; j < n; ++j)
                {
                    Q[j, i] = currQCol[j];
                }
            }

            R = TransposeMatrix(Q) * A;
        }
        public static void SLAE_QR(Matrix A, Matrix b, out Vector X)
        {
            int n = A.GetLength(0);
            MyMath.QRdecomposition(A, out Matrix Q, out Matrix R);
            Matrix y = MyMath.TransposeMatrix(Q) * b;
            X = new Vector(n);
            X[n - 1] = y[n - 1, 0] / R[n - 1, n - 1];
            for (int i = n - 2; i >= 0; --i)
            {
                X[i] = y[i, 0];
                for (int j = n - 1; j > i; --j)
                {
                    X[i] -= R[i, j] * X[j];
                }
                X[i] /= R[i, i];
            }
        }
        private static bool SeidelTerminationCondition(Vector Xcur, Vector Xprev)
        {
            int n = Xcur.Lenght;
            for (int i = 0; i < n; ++i)
            {
                if (!DoubleCompare(Xcur[i], Xprev[i]))
                    return false;
            }
            return true;
        }
        public static bool IsMatrixDiagDomin(Matrix A)
        {
            int n = A.GetLength(0);
            double sum = 0;
            for (int i = 0; i < n; ++i)
            {
                sum = 0;
                for(int j = 0; j < n; ++j)
                    if (i != j)
                        sum += Math.Abs(A[i, j]);
                if(sum < A[i, i])
                    return false;
            }
            return true;
        }
        public static void SeidelMethod(Matrix A, Matrix b, out Vector X, out int iterCount, out double aPrioriEst)
        {
            int n = A.GetLength(0);
            X = new Vector(n);
            Vector Xprev = new Vector(n);
            for(int i = 0; i < n; ++i)
            {
                X[i] = 1;
            }
            Matrix Q = new Matrix(n);
            for(int i = 0; i < n; ++i)
            {
                for(int j = 0; j < i; ++j)
                    Q[i, j] = A[i, j] / A[i, i];
                for(int j = i+1; j < n-1; ++j)
                    Q[i, j - 1] = A[i, j] / A[i, i];
            }
            for(int i = 0; i < n; ++i)
            {
                Q[i, n - 2] = b[i, 0] / A[i, i];
            }

            double qNorm = MatrixNorm(Q);

            aPrioriEst = 1;
            iterCount = 0;
            do
            {
                if(iterCount == 1)
                {
                    aPrioriEst = (int)Math.Log(eps * (1 - qNorm) / VectorNorm(X - Xprev), qNorm);
                }
                for(int i = 0; i < n; ++i)
                {
                    Xprev[i] = X[i];
                }
                for(int i = 0; i < n; ++i)
                {
                    double a = 0;
                    for(int j = 0; j < i; ++j)
                        a += (A[i, j] * X[j]);
                    for(int j = i+1; j < n; ++j)
                        a += (A[i, j] * Xprev[j]);
                    X[i] = (b[i,0] - a) / A[i, i];
                }
                ++iterCount;
            } while (!SeidelTerminationCondition(X, Xprev));
        }
        public static void JacobiMethod(Matrix A, Matrix b, out Vector X, out int iterCount, out double aPrioriEst)
        {
            int n = A.GetLength(0);
            double norm = 0;
            X = new Vector(n);
            for (int i = 0; i < n; ++i)
            {
                X[i] = 1;
            }
            Vector Xprev = new Vector(n);
            Matrix Q = new Matrix(n);
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < i; ++j)
                    Q[i, j] = A[i, j] / A[i, i];
                for (int j = i + 1; j < n - 1; ++j)
                    Q[i, j - 1] = A[i, j] / A[i, i];
            }
            for (int i = 0; i < n; ++i)
            {
                Q[i, n - 2] = b[i, 0] / A[i, i];
            }

            double qNorm = MatrixNorm(Q);
            aPrioriEst = 0;
            iterCount = 0;
            do
            {
                for(int i = 0; i < n; ++i)
                {
                    Xprev[i] = b[i,0];
                    for(int j = 0; j < n; ++j)
                        if(i != j)
                            Xprev[i] -= A[i, j] * X[j];
                    Xprev[i] /= A[i, i];
                }
                norm = 0;
                if (iterCount == 0)
                {
                    aPrioriEst = (int)Math.Log(eps * (1 - qNorm) / VectorNorm(X - Xprev), qNorm);
                }
                for (int i = 0; i < n; i++)
                {
                    norm = Math.Max(norm, Math.Abs(X[i] - Xprev[i]));
                    X[i] = Xprev[i];
                }
                ++iterCount;
            } while (norm > eps);
        }

    }
}
