using System;
using System.IO;

namespace FEM
{
    /// <summary>
    /// метод GMRES
    /// </summary>
    [Serializable] class GMRES : ISolver
    {
        /// <summary>
        /// ссылка на решаемую СЛАУ
        /// </summary>
        public ISLAE SLAE { set; get; }
        /// <summary>
        /// предобусловливатель
        /// </summary>
        public IPreconditioner Preconditioner { set; get; }
        /// <summary>
        /// порог останова по невязке
        /// </summary>
        public double Eps { set; get; }
        /// <summary>
        /// максимальное число итераций
        /// </summary>
        public int Max_Iter { set; get; }
        /// <summary>
        /// текущая итерация
        /// </summary>
        public int Iter { set; get; }
        /// <summary>
        /// текущая норма невязки
        /// </summary>
        public double Norma_Residual { set; get; }
        /// <summary>
        /// флаг останова алгоритмов (по умолчанию 0)
        /// </summary>
        public int Flag { set; get; }
        /// <summary>
        /// результат
        /// </summary>
        public double[] Res { set; get; }

        /// <summary>
        /// глубина метода GMRES
        /// </summary>
        public int m { set; get; }

        /// <summary>
        /// вспомогательные матрицы метода
        /// </summary>
        Matrix H, H_, V, W;

        /// <summary>
        /// конструктор метода GMRES
        /// </summary>
        /// <param name="M - глубина"></param>
        public GMRES (int M = 5) { m = M; }

        /// <summary>
        /// старт алгоритма решателя GMRES
        /// </summary>
        public void Start_Method()
        {
            //размер системы
            int n = SLAE.n;

            //вспомогательные матрицы метода
            H  = new Matrix(m + 1, m);
            H_ = new Matrix(m + 1, m);
            V  = new Matrix(m + 1, n);
            W  = new Matrix(m, n);

            //вспомогательные векторы
            var r = new Vector(n);
            var Help = new Vector(n);
            var Minimization_Result = new Vector(m + 1);
            
            //параметры метода: флаг останова и глубина GMRES
            int flag_m = 0, old_m = m;

            //норма невязкаи
            double Norma_r;

            //вычисляем вектор невязки
            SLAE.MultMV(Res, Help.Elem);
            for (int i = 0; i < n; i++)
            {
                Help.Elem[i] = SLAE.f[i] - Help.Elem[i];
            }

            //процедура предобусловливания r = M^(-1) * (F - A * x)
            Preconditioner.Start_Preconditioner(Help, r);

            //норма вектора невязки
            Norma_r = r.Norma();

            //итерации метода GMRES
            while (Flag == 0 && Iter < Max_Iter)
            {
                //восстановление параметров
                m = old_m;
                flag_m = 0;

                //ЗДЕСЬ БЫЛО ЗАНУЛЕНИЕ МАТРИЦ H V W
                H.Clear_Matrix();
                V.Clear_Matrix();
                W.Clear_Matrix();
                H_.Clear_Matrix();
                

                //строим первый вектор подпространства Крылова v1 = r / ||r||
                for (int i = 0; i < n; i++) V.Elem[0][i] = r.Elem[i] / Norma_r;

                //строим оставшиеся базисные векторы
                for (int j = 0; j < m && flag_m == 0; j++)
                {
                    //Help = A * Vj
                    SLAE.MultMV(V.Elem[j], Help.Elem);
                    
                    //w = Help = M^(-1) * Help
                    Preconditioner.Start_Preconditioner(Help, Help);

                    //процедура ортогонализации Арнольди
                    for (int i = 0; i <= j; i++)
                    {
                        //скалярное произведение Hij = Wj * Vi, где Wj = Help
                        H.Elem[i][j] = 0.0;
                        for (int l = 0; l < n; l++)
                        {
                            W.Elem[j][l] = Help.Elem[l];
                            H.Elem[i][j] += W.Elem[j][l] * V.Elem[i][l];
                        }
                        
                        //в H_ хранится не модифицированная матрица коэффициентов ортогонализации
                        H_.Elem[i][j] = H.Elem[i][j];
                        
                        //Wj = Wj - Hij * Vi
                        for (int l = 0; l < n; l++)
                        {
                            W.Elem[j][l] -= H.Elem[i][j] * V.Elem[i][l];
                        }
                    }

                    //расширение матрицы H и H_
                    H.Elem[j + 1][j] = Operation.Norma(W.Elem[j]);
                    H_.Elem[j + 1][j] = H.Elem[j + 1][j];

                    if (Math.Abs(H.Elem[j + 1][j]) < CONST.EPS)
                    {
                        m = j;
                        flag_m = 1;
                    }
                    else
                    {
                        for (int l = 0; l < n; l++) V.Elem[j + 1][l] = W.Elem[j][l] / H.Elem[j + 1][j];
                    }
                }
                
                /*
                //проверка на ортогональность
                for (int I = 0; I < V.M; I++)
                {
                    for (int J = I; J < V.M; J++)
                    {
                        Console.WriteLine("V{0} * V{1} = {2}", I+1, J+1, Operation.ScalMult(V.Elem[I], V.Elem[J]));
                        Console.ReadLine();
                    }
                }
                */

                //решение задачи минимизации (МНК) и норма невязки
                Norma_Residual = Minimization_Problem(Norma_r, Minimization_Result.Elem);

                if (Norma_Residual < Eps) Flag = 1;
                else
                {
                    //вычисляем произведение Vy = V(t) * y
                    Mult_Vy(Minimization_Result.Elem, Help.Elem);

                    //вычисление нового результата
                    for (int l = 0; l < n; l++) { Res[l] += Help.Elem[l]; }

                    /*
                    //вычисление вектора невязки r = V * (||r||e1 - (H_)m * y)
                    for (int i = 0; i < m; i++)
                    {
                        r.Elem[i] = 0.0;
                        for (int j = 0; j < m; j++) r.Elem[i] -= H_.Elem[i][j] * Minimization_Result.Elem[j];
                    }
                    r.Elem[0] += Norma_r;

                    for (int i = 0; i < n; i++)
                    {
                        Help.Elem[i] = 0.0;
                        for (int j = 0; j < m; j++)
                        {
                            Help.Elem[i] += V.Elem[j][i] * r.Elem[j];
                        }
                    }
                    */
                    
                    SLAE.MultMV(Res, Help.Elem);
                    for (int i = 0; i < n; i++)
                    {
                        Help.Elem[i] = SLAE.f[i] - Help.Elem[i];
                    }
                    

                    //процедура предобусловливания r = M^(-1) * Help, где Help = f - A * x
                    Preconditioner.Start_Preconditioner(Help, r);

                    Norma_r = r.Norma();

                    Iter++;
                } 
            }
        }

        /// <summary>
        /// процедура минимизации
        /// возвращает норму невязки
        /// </summary>
        /// <param name="Norma_r - текущая норма невязки"></param>
        /// <param name="f - результат минимизации"></param>
        double Minimization_Problem(double Norma_r, double[] f)
        {
            //вектор правой части для решения СЛАУ с треуг.матрицей f = ||r|| * e1
            for (int i = 1; i < m + 1; i++) f[i] = 0.0;
            f[0] = Norma_r;
            
            //вспомогательные переменные
            double help1, help2;
           
            //косинус, синус, тангенс
            double c, s, t;

            //преобразования Гивенса: приведение матрицы Хессенберга H к верхней треугольной форме
            for (int i = 0; i < m; i++)
            {
                t = H.Elem[i + 1][i] / H.Elem[i][i];
                c = 1.0 / Math.Sqrt(1 + t * t);
                s = t * c;

                //H_new = Gt * H (минус у матрицы вращения внизу)
                for (int k = i; k < H.N; k++)
                {
                    help1 = c * H.Elem[i][k] + s * H.Elem[i + 1][k];
                    help2 = c * H.Elem[i + 1][k] - s * H.Elem[i][k];
                    H.Elem[i][k] = help1;
                    H.Elem[i + 1][k] = help2;
                }
                
                //перемножаем слева вектор правой части на трансп.матрицу преобразования
                help1 = c * f[i] + s * f[i + 1];
                help2 = c * f[i + 1] - s * f[i];
                f[i] = help1;
                f[i + 1] = help2;
            }

            //решается система m X m с верхней треугольной матрицей
            for (int i = m - 1; i >= 0; i--)
            {
                f[i] /= H.Elem[i][i];
                for (int j = i - 1; j >= 0; j--)
                    f[j] -= H.Elem[j][i] * f[i];
            }

            return Math.Abs(f[m]);
        }

        /// <summary>
        /// вычисление произведения Vy = V(t) * y
        /// </summary>
        /// <param name="y - вектор-множитель"></param>
        /// <param name="Vy - результат умножения"></param>
        void Mult_Vy(double[] y, double[] Vy)
        {
            //размер системы
            int n = SLAE.n;
            //вспомогательная переменная
            double s;
            for (int i = 0; i < n; i++)
            {
                s = 0;
                for (int j = 0; j < m; j++)
                    s += V.Elem[j][i] * y[j];
                Vy[i] = s;
            }
        }
    }
}
