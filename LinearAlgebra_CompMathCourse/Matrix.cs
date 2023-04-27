using System;
using System.Text;

namespace MyMathematica
{
    public class Matrix
    {
        private double[,] matrix;
        public Matrix()
        {
            matrix = new double[0, 0];
        }
        public Matrix(double[,] m)
        {
            matrix = new double[m.GetLength(0), m.GetLength(1)];
            for (int i = 0; i < m.GetLength(0); ++i)
            {
                for (int j = 0; j < m.GetLength(1); ++j)
                {
                    matrix[i, j] = m[i, j];
                }
            }
        }
        public Matrix(int row, int col)
        {
            if (row < 0 || col < 0)
            {
                throw new Exception("Error: Wrong size.");
            }
            matrix = new double[row, col];
        }
        public Matrix(int n)
        {
            if (n < 0)
            {
                throw new Exception("Error: Wrong size.");
            }
            matrix = new double[n, n];
        }
        public Matrix(Matrix matrixToCopy)
        {
            int n = matrixToCopy.GetLength(0);
            int m = matrixToCopy.GetLength(1);
            matrix = new double[n, m];
            for(int i = 0; i < n; ++i)
            {
                for (int j = 0; j < m; ++j)
                {
                    matrix[i, j] = matrixToCopy[i, j];
                }
            }
        }
        public static Matrix Eye(int n)
        {
            Matrix matrix = new Matrix(n);

            for (int i = 0; i < n; ++i)
            {
                matrix[i, i] = 1;
            }

            return matrix;
        }
        public static Matrix Eye(int row, int col)
        {
            Matrix matrix = new Matrix(row, col);
            for (int i = 0; i < Math.Min(row, col); ++i)
            {
                matrix[i, i] = 1;
            }

            return matrix;
        }
        public static Matrix Ones(int n)
        {
            Matrix matrix = new Matrix(n, n);
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    matrix[i, j] = 1;
                }
            }
            return matrix;
        }
        public static Matrix Ones(int row, int col)
        {
            Matrix matrix = new Matrix(row, col);
            for (int i = 0; i < row; ++i)
            {
                for (int j = 0; j < col; ++j)
                {
                    matrix[i, j] = 1;
                }
            }
            return matrix;
        }
        public void Zeros()
        {
            for(int i = 0; i < matrix.GetLength(0); ++i)
            {
                for (int j = 0; j < matrix.GetLength(0); ++j)
                {
                    matrix[i, j] = 0;
                }
            }
        }
        public void PrintMatrix()
        {
            for (int i = 0; i < matrix.GetLength(0); ++i)
            {
                for (int j = 0; j < matrix.GetLength(1); ++j)
                {
                    Console.Write(matrix[i, j] + "\t");
                }
                Console.WriteLine();
            }
            Console.WriteLine();
        }

        public override string ToString()
        {
            int n = GetLength(0);
            int m = GetLength(1);
            StringBuilder str = new StringBuilder(n * m);
            for(int i = 0; i < n; ++i)
            {
                for(int j = 0; j < m; ++j)
                {
                    str.Append($"{matrix[i, j]} ");
                }
                str.Append('\n');
            }
            return str.ToString();
        }
        public int GetLength(int n) // 0 при n != 1 || n != 0
        {
            if (n != 1 && n != 0) return 0;
            else return matrix.GetLength(n);
        }
        public double this[int i, int j] { 
            get => matrix[i, j];
            set => matrix[i, j] = value;
        }
        public static Matrix operator*(Matrix A, Matrix B)
        {
            if(A.GetLength(1) != B.GetLength(0))
            {
                throw new Exception("Error: Wrong Size in Matrix multiplication ");
            }
            int n = A.GetLength(0);
            int m = B.GetLength(1);
            int l = A.GetLength(1);
            Matrix result = new Matrix(n,m);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    result[i, j] = 0;

                    for (int k = 0; k < l; k++)
                    {
                        result[i, j] += A[i, k] * B[k, j];
                    }
                }
            }

            return result;
        }

        public override bool Equals(object obj)
        {
            return base.Equals(obj);
        }
        public override int GetHashCode()
        {
            return base.GetHashCode();
        }
        public static bool operator ==(Matrix A, Matrix B)
        {
            if (A.GetLength(0) != B.GetLength(0) || A.GetLength(1) != B.GetLength(1)) 
                return false;
            int n = A.GetLength(0);
            int m = A.GetLength(1);
            for(int i = 0; i < n; ++i)
            {
                for(int j = 0; j < m; ++j)
                {
                    if (MyMath.DoubleCompare(A[i, j], B[i, j]) != true)
                    {
                        Console.WriteLine($"{A[i, j]} {B[i, j]}");
                        return false;
                    }
                }
            }

            return true;

        }
        public static bool operator !=(Matrix A, Matrix B)
        {
            return !(A == B);
        }
        public static Matrix GenerateRandomMatrix(int n, int m, int min = int.MinValue, int max = int.MaxValue)
        {
            Matrix result = new Matrix(n, m);
            Random rand = new Random();
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < m; ++j)
                {
                    result[i, j] = rand.Next(min, max);
                }
            }

            return result;
        }
        public void SwapRows(int i, int j)
        {
            int n = GetLength(1);
            double tmp;
            for(int k = 0; k < n; ++k)
            {
                tmp = matrix[i, k];
                matrix[i, k] = matrix[j, k];
                matrix[j, k] = tmp;
            }
        }
        public void SwapColumns(int i, int j)
        {
            int n = GetLength(0);
            double tmp;
            for (int k = 0; k < n; ++k)
            {
                tmp = matrix[k, i];
                matrix[k, i] = matrix[k, j];
                matrix[k, j] = tmp;
            }
        }
        public void MakeRowsPermutations(Vector P)
        {
            if (P.Lenght != GetLength(0))
            {
                throw new Exception("Error: Wrong size.");
            }
            int n = GetLength(0);
            
            for(int i = 0; i < n; ++i)
            {
                while(P[i] != i)
                {
                    SwapRows(i, (int)P[i]);
                    P.SwapElements(i, (int)P[i]);
                }
            }
        }
        public void DeleteRow(int row)
        {
            if (row < 0 && row >= GetLength(0))
                throw new Exception("Error. This row doesn't exist.");
            double[,] newMatrix = new double[GetLength(0) - 1, GetLength(1)];
            for (int i = 0; i < row; ++i)
            {
                for (int j = 0; j < GetLength(1); ++j)
                {
                    newMatrix[i, j] = matrix[i, j];
                }
            }

            for (int i = row + 1; i < GetLength(0); ++i)
            {
                for (int j = 0; j < GetLength(1); ++j)
                {
                    newMatrix[i-1, j] = matrix[i, j];
                }
            }

            matrix = newMatrix;
        }
    }
}
