using System;

namespace MyMathematica
{
    public class Vector
    {
        private double[] vector;
        public Vector()
        {
            vector = new double[0];
        }
        public Vector(int n)
        {
            vector = new double[n];
        }
        public Vector(double [] v)
        {
            vector = new double[v.Length];
            for(int i = 0; i < v.Length; ++i)
            {
                vector[i] = v[i];
            }
        }
        public Vector(Vector v)
        {
            vector = new double[v.Lenght];
            for (int i = 0; i < v.Lenght; ++i)
            {
                vector[i] = v[i];
            }
        }
        public void PrintVector()
        {
            for (int i = 0; i < vector.GetLength(0); ++i)
            {
                Console.WriteLine(vector[i]);
            }
        }
        public double this[int index]
        {
            get => vector[index];
            set => vector[index] = value;
        }
        public void SwapElements(int i, int j)
        {
            double tmp = vector[i];
            vector[i] = vector[j];
            vector[j] = tmp;
        }
        public int Lenght { get { return vector.Length; } }

        public static Vector operator +(Vector a, Vector b)
        {
            int n = a.Lenght;
            Vector res = new Vector(n);
            for (int i = 0; i < n; ++i)
            {
                res[i] = a[i] + b[i];
            }
            return res;
        }
        public static Vector operator -(Vector a, Vector b)
        {
            int n = a.Lenght;
            Vector res = new Vector(n);
            for (int i = 0; i < n; ++i)
            {
                res[i] = a[i] - b[i];
            }
            return res;
        }
        public static double operator*(Vector a, Vector b)
        {
            double res = 0;
            int n = a.Lenght;
            for(int i = 0; i < n; ++i)
            {
                res += a[i] * b[i];
            }
            return res;
        }

        public static Vector operator *(double num, Vector vec)
        {
            int n = vec.Lenght;
            Vector res = new Vector(vec);
            for(int i = 0; i < n; ++i)
            {
                res[i] *= num;
            }
            return res;
        }

        public static Vector operator *(Vector vec, double num)
        {
            int n = vec.Lenght;
            Vector res = new Vector(vec);
            for (int i = 0; i < n; ++i)
            {
                res[i] *= num;
            }
            return res;
        }

        public static Vector operator /(double num, Vector vec)
        {
            int n = vec.Lenght;
            Vector res = new Vector(vec);
            for (int i = 0; i < n; ++i)
            {
                res[i] /= num;
            }
            return res;
        }

        public static Vector operator /(Vector vec, double num)
        {
            int n = vec.Lenght;
            Vector res = new Vector(vec);
            for (int i = 0; i < n; ++i)
            {
                res[i] /= num;
            }
            return res;
        }



    }
}
