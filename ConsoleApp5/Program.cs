using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ConsoleApp5
{
    internal class Program
    {
        static void Main(string[] args)
        {
            void VectorLength(double px,double py,double pz, out double pl)
            {
                pl = Math.Sqrt(px * px + py * py + pz * pz);
            }
            void AngleVect(double px1, double py1, double pz1,double px2,double py2, double pz2, out double pa)
            {
                double kinter;
                double ar;
                kinter = (px1 * px2 + py1 * py2 + pz1 * pz2) / (Math.Sqrt(px1 * px1 + py1 * py1 + pz1 * pz1) * Math.Sqrt(px2 * px2 + py2 * py2 + pz2 * pz2));
                ar = Math.Acos(kinter);
                pa = (180 / Math.PI) * ar;
            }
            void SumVect(double px1, double py1, double pz1, double px2, double py2, double pz2, out double px, out double py, out double pz)            
            {
                    px = px1 + px2;
                    py = py1 + py2;
                    pz = pz1 + pz2;
            }
            void SubstructVect(double px1, double py1, double pz1, double px2, double py2, double pz2, out double px, out double py, out double pz)
            {
                px = px1 - px2;
                py = py1 - py2;
                pz = pz1 - pz2;
            }
            void MultVectOnNum(double px1, double py1, double pz1, double pk, out double px, out double py, out double pz)            
            {
                px = pk * px1;
                py = pk * py1;
                pz = pk * pz1;
            }
            void OpposVect(double px1, double py1, double pz1, out double px, out double py, out double pz)
            {
                px = -px1;
                py = -py1;
                pz = -pz1;
            }

            void ScalarMultiply(double px1, double py1, double pz1, double px2, double py2, double pz2, out double psc)
            {
                psc = px1 * px2 + py1 * py2 + pz1 * pz2;
            }
            void VectorMultiply(double px1, double py1, double pz1, double px2, double py2, double pz2, out double px, out double py, out double pz)
            {
                px = py1 * pz2 - py2 * pz1;
                py = px2 * pz1 - px1 * pz2;
                pz = px1 * py2 - px2 * py1;
            }
            void MixedMultVect(double px1, double py1, double pz1, double px2, double py2, double pz2,  double px3, double py3,  double pz3, out double pm)
            {
                pm = px1 * py2 * pz3 + px2 * py3 * pz1 + px3 * py1 * pz2 - px3 * py2 * pz1 - px1 * py3 * pz2 - px2 * py1 * pz3;
            }

            string st;

            double ax; // вектор a
            double ay;
            double az;

            double bx; // вектор b
            double by;
            double bz;

            double cx; // вектор c
            double cy;
            double cz;

            double dx; // вектор d
            double dy;
            double dz;

            double rx; // векторный результат
            double ry;
            double rz;

            double rx1; // векторный результат левой части равенства
            double ry1;
            double rz1;

            double rx2; // векторный результат правой части равенства
            double ry2;
            double rz2;

            double r;  // скалярный результат
            double r1; // скалярный результат левой части равенства
            double r2; // скалярный результат правой части равенства

            double vx1; //  дополнительный вектор 1
            double vy1;
            double vz1;

            double vx2; //  дополнительный вектор 2
            double vy2;
            double vz2;

            double k; // скалярные коэффициенты
            double k1;
            double k2;
            double k3;

            int sel;  // выбор пункта главного меню
            int sel1; // выбор пункта подменю главного меню

            while (true)
            {
                Console.WriteLine("         Главное меню");
                Console.WriteLine("1 - Общие свойства векторов");
                Console.WriteLine("2 - Свойства скалярного произведения векторов");
                Console.WriteLine("3 - Свойства векторного произведения векторов");
                Console.WriteLine("4 - Свойства смешанного произведения векторов");
                Console.WriteLine("5 - Вычисления с векторами");
                Console.WriteLine("6 - Выход из программы");
                Console.Write("Введите пункт меню: ");
                sel = int.Parse(Console.ReadLine());
                Console.WriteLine();
                switch (sel)
                {
                    /* общие свойства векторов */
                    case 1:
                        Console.WriteLine("      Меню общие свойства векторов");
                        Console.WriteLine("1 - Коммутативность сложения векторов");
                        Console.WriteLine("2 - Ассоциативность сложения векторов");
                        Console.WriteLine("3 - Ассоциативность умножения вектора на число");
                        Console.WriteLine("4 - Дистрибутивность умножения вектора на число относительно сложения векторов");
                        Console.WriteLine("5 - Дистрибутивность умножения вектора на число относительно сложения чисел");
                        Console.Write("Введите пункт меню: ");
                        sel1 = int.Parse(Console.ReadLine());
                        Console.WriteLine();
                        switch (sel1)                        {
                            
                            case 1:
                                /* va+vb=vb+va */
                                Console.WriteLine("Проверка коммутативности сложения векторов");
                                Console.WriteLine("va+vb=vb+va");
                                Console.Write("ax: ");
                                st = Console.ReadLine();
                                ax = double.Parse(st);
                                Console.Write("ay: ");
                                st = Console.ReadLine();
                                ay = double.Parse(st);
                                Console.Write("az: ");
                                st = Console.ReadLine();
                                az = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("bx: ");
                                st = Console.ReadLine();
                                bx = double.Parse(st);
                                Console.Write("by: ");
                                st = Console.ReadLine();
                                by = double.Parse(st);
                                Console.Write("bz: ");
                                st = Console.ReadLine();
                                bz = double.Parse(st);
                                Console.WriteLine();
                                SumVect(ax, ay, az, bx, by, bz, out rx1, out ry1, out rz1);
                                SumVect(bx, by, bz, ax, ay, az, out rx2, out ry2, out rz2);
                                Console.WriteLine($" Левая часть равенства: rx1={rx1}   ry1={ry1}   rz1={rz1}");
                                Console.WriteLine($"Правая часть равенства: rx2={rx2}   ry2={ry2}   rz2={rz2}");
                                break;
                            case 2:
                                /* (va+vb)+vc=va+(vb+vc) */
                                Console.WriteLine("Проверка ассоциативности сложения векторов");
                                Console.WriteLine("(va+vb)+vc=va+(vb+vc)");
                                Console.Write("ax: ");
                                st = Console.ReadLine();
                                ax = double.Parse(st);
                                Console.Write("ay: ");
                                st = Console.ReadLine();
                                ay = double.Parse(st);
                                Console.Write("az: ");
                                st = Console.ReadLine();
                                az = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("bx: ");
                                st = Console.ReadLine();
                                bx = double.Parse(st);
                                Console.Write("by: ");
                                st = Console.ReadLine();
                                by = double.Parse(st);
                                Console.Write("bz: ");
                                st = Console.ReadLine();
                                bz = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("cx: ");
                                st = Console.ReadLine();
                                cx = double.Parse(st);
                                Console.Write("cy: ");
                                st = Console.ReadLine();
                                cy = double.Parse(st);
                                Console.Write("cz: ");
                                st = Console.ReadLine();
                                cz = double.Parse(st);
                                Console.WriteLine();
                                SumVect(ax, ay, az, bx, by, bz, out vx1, out vy1, out vz1);
                                SumVect(vx1, vy1, vz1, cx, cy, cz, out rx1, out ry1, out rz1);
                                SumVect(bx, by, bz, cx, cy, cz, out vx2, out vy2, out vz2);
                                SumVect(ax, ay, az, vx2, vy2, vz2, out rx2, out ry2, out rz2);
                                Console.WriteLine($" Левая часть равенства: rx1={rx1}   ry1={ry1}   rz1={rz1}");
                                Console.WriteLine($"Правая часть равенства: rx2={rx2}   ry2={ry2}   rz2={rz2}");
                                break;
                            case 3:
                                /*  (k1*k2)*va=k1*(k2*va) */
                                Console.WriteLine("Проверка ассоциативности умножения вектора на число");
                                Console.WriteLine("(k1*k2)*va=k1*(k2*va)");
                                Console.Write("k1: ");
                                st = Console.ReadLine();
                                k1 = double.Parse(st);
                                Console.Write("k2: ");
                                st = Console.ReadLine();
                                k2 = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("ax: ");
                                st = Console.ReadLine();
                                ax = double.Parse(st);
                                Console.Write("ay: ");
                                st = Console.ReadLine();
                                ay = double.Parse(st);
                                Console.Write("az: ");
                                st = Console.ReadLine();
                                az = double.Parse(st);
                                Console.WriteLine();
                                k3 = k1 * k2;
                                MultVectOnNum(ax, ay, az, k3, out rx1, out ry1, out rz1);
                                MultVectOnNum(ax, ay, az, k2, out vx1, out vy1, out vz1);
                                MultVectOnNum(vx1, vy1, vz1, k1, out rx2, out ry2, out rz2);
                                Console.WriteLine($" Левая часть равенства: rx1={rx1}   ry1={ry1}   rz1={rz1}");
                                Console.WriteLine($"Правая часть равенства: rx2={rx2}   ry2={ry2}   rz2={rz2}");
                                break;
                            case 4:
                                /* k1*(va+vb)=k1*va+k1*vb */
                                Console.WriteLine("Проверка дистрибутивности умножения вектора на число относительно сложения векторов");
                                Console.WriteLine("k1*(va+vb)=k1*va+k1*vb");
                                Console.Write("k1: ");
                                st = Console.ReadLine();
                                k1 = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("ax: ");
                                st = Console.ReadLine();
                                ax = double.Parse(st);
                                Console.Write("ay: ");
                                st = Console.ReadLine();
                                ay = double.Parse(st);
                                Console.Write("az: ");
                                st = Console.ReadLine();
                                az = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("bx: ");
                                st = Console.ReadLine();
                                bx = double.Parse(st);
                                Console.Write("by: ");
                                st = Console.ReadLine();
                                by = double.Parse(st);
                                Console.Write("bz: ");
                                st = Console.ReadLine();
                                bz = double.Parse(st);
                                Console.WriteLine();
                                SumVect(ax, ay, az, bx, by, bz, out vx1, out vy1, out vz1);
                                MultVectOnNum(vx1, vy1, vz1, k1, out rx1, out ry1, out rz1);
                                MultVectOnNum(ax, ay, az, k1, out vx1, out vy1, out vz1);
                                MultVectOnNum(bx, by, bz, k1, out vx2, out vy2, out vz2);
                                SumVect(vx1, vy1, vz1, vx2, vy2, vz2, out rx2, out ry2, out rz2);
                                Console.WriteLine($" Левая часть равенства: rx1={rx1}   ry1={ry1}   rz1={rz1}");
                                Console.WriteLine($"Правая часть равенства: rx2={rx2}   ry2={ry2}   rz2={rz2}");
                                break;
                            case 5:
                                /* (k1+k2)*va=k1*va+k2*va */
                                Console.WriteLine("Проверка дистрибутивности умножения вектора на число относительно сложения чисел");
                                Console.WriteLine("(k1+k2)*va=k1*va+k2*va");
                                Console.Write("k1: ");
                                st = Console.ReadLine();
                                k1 = double.Parse(st);
                                Console.Write("k2: ");
                                st = Console.ReadLine();
                                k2 = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("ax: ");
                                st = Console.ReadLine();
                                ax = double.Parse(st);
                                Console.Write("ay: ");
                                st = Console.ReadLine();
                                ay = double.Parse(st);
                                Console.Write("az: ");
                                st = Console.ReadLine();
                                az = double.Parse(st);
                                Console.WriteLine();
                                k3 = k1 + k2;
                                MultVectOnNum(ax, ay, az, k3, out rx1, out ry1, out rz1);
                                MultVectOnNum(ax, ay, az, k1, out vx1, out vy1, out vz1);
                                MultVectOnNum(ax, ay, az, k2, out vx2, out vy2, out vz2);
                                SumVect(vx1, vy1, vz1, vx2, vy2, vz2, out rx2, out ry2, out rz2);
                                Console.WriteLine($" Левая часть равенства: rx1={rx1}   ry1={ry1}   rz1={rz1}");
                                Console.WriteLine($"Правая часть равенства: rx2={rx2}   ry2={ry2}   rz2={rz2}");
                                break;
                            default:
                                Console.WriteLine("Введён неверный пункт меню!");
                                break;
                        }
                        break;
                        /* скалярное произведение */
                    case 2:
                        Console.WriteLine("     Меню свойства скалярного произведения векторов");
                        Console.WriteLine("1 - Коммутативность скалярного произведения векторов");
                        Console.WriteLine("2 - Ассоциативность скалярного произведения векторов совместно с умножением вектора на число");
                        Console.WriteLine("3 - Дистрибутивность скалярного произведения векторов относительно сложения векторов");
                        Console.Write("Введите пункт меню: ");
                        sel1 = int.Parse(Console.ReadLine());
                        Console.WriteLine();
                        switch (sel1)
                        {
                            case 1:
                                /*  (va*vb)=(vb*va) */
                                Console.WriteLine("Проверка коммутативности скалярного умножения векторов");
                                Console.WriteLine("(va*vb)=(vb*va)");
                                Console.Write("ax: ");
                                st = Console.ReadLine();
                                ax = double.Parse(st);
                                Console.Write("ay: ");
                                st = Console.ReadLine();
                                ay = double.Parse(st);
                                Console.Write("az: ");
                                st = Console.ReadLine();
                                az = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("bx: ");
                                st = Console.ReadLine();
                                bx = double.Parse(st);
                                Console.Write("by: ");
                                st = Console.ReadLine();
                                by = double.Parse(st);
                                Console.Write("bz: ");
                                st = Console.ReadLine();
                                bz = double.Parse(st);
                                Console.WriteLine();
                                ScalarMultiply(ax, ay, az, bx, by, bz, out r1);
                                ScalarMultiply(bx, by, bz, ax, ay, az, out r2);
                                Console.WriteLine($" Левая часть равенства: r1={r1}");
                                Console.WriteLine($"Правая часть равенства: r2={r2}");
                                break;
                            case 2:
                                /* (k1*va)*vb=k1*(va*vb) */
                                Console.WriteLine("Проверка ассоциативности скалярного произведения векторов совместно с умножением вектора на число");
                                Console.WriteLine("(k1*va)*vb=k1*(va*vb)");
                                Console.Write("k1: ");
                                st = Console.ReadLine();
                                k1 = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("ax: ");
                                st = Console.ReadLine();
                                ax = double.Parse(st);
                                Console.Write("ay: ");
                                st = Console.ReadLine();
                                ay = double.Parse(st);
                                Console.Write("az: ");
                                st = Console.ReadLine();
                                az = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("bx: ");
                                st = Console.ReadLine();
                                bx = double.Parse(st);
                                Console.Write("by: ");
                                st = Console.ReadLine();
                                by = double.Parse(st);
                                Console.Write("bz: ");
                                st = Console.ReadLine();
                                bz = double.Parse(st);
                                Console.WriteLine();
                                MultVectOnNum(ax, ay, az, k1, out vx1, out vy1, out vz1);
                                ScalarMultiply(vx1, vy1, vz1, bx, by, bz, out r1);
                                ScalarMultiply(ax, ay, az, bx, by, bz, out k2);
                                r2 = k1 * k2;
                                Console.WriteLine($" Левая часть равенства: r1={r1}");
                                Console.WriteLine($"Правая часть равенства: r2={r2}");
                                break;
                            case 3:
                                /*  (va+vb)*vc=(va*vc)+(vb*vc) */
                                Console.WriteLine("Проверка дистрибутивности скалярного произведения векторов относительно сложения векторов");
                                Console.WriteLine("(va+vb)*vc=(va*vc)+(vb*vc)");
                                Console.Write("ax: ");
                                st = Console.ReadLine();
                                ax = double.Parse(st);
                                Console.Write("ay: ");
                                st = Console.ReadLine();
                                ay = double.Parse(st);
                                Console.Write("az: ");
                                st = Console.ReadLine();
                                az = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("bx: ");
                                st = Console.ReadLine();
                                bx = double.Parse(st);
                                Console.Write("by: ");
                                st = Console.ReadLine();
                                by = double.Parse(st);
                                Console.Write("bz: ");
                                st = Console.ReadLine();
                                bz = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("cx: ");
                                st = Console.ReadLine();
                                cx = double.Parse(st);
                                Console.Write("cy: ");
                                st = Console.ReadLine();
                                cy = double.Parse(st);
                                Console.Write("cz: ");
                                st = Console.ReadLine();
                                cz = double.Parse(st);
                                Console.WriteLine();
                                SumVect(ax, ay, az, bx, by, bz, out vx1, out vy1, out vz1);
                                ScalarMultiply(vx1, vy1, vz1, cx, cy, cz, out r1);
                                ScalarMultiply(ax, ay, az, cx, cy, cz, out k1);
                                ScalarMultiply(bx, by, bz, cx, cy, cz, out k2);
                                r2 = k1 + k2;
                                Console.WriteLine($" Левая часть равенства: r1={r1}");
                                Console.WriteLine($"Правая часть равенства: r2={r2}");
                                break;
                            default:
                                Console.WriteLine("Введён неверный пункт меню!");
                                break;
                        }
                        break;
                        /* векторное произведение */
                    case 3:
                        Console.WriteLine("     Меню свойства векторного произведения");
                        Console.WriteLine("1 - Антикоммутативность векторного произведения");
                        Console.WriteLine("2 - Ассоциативность векторного произведения совместно с умножением вектора на число");
                        Console.WriteLine("3 - Дистрибутивность векторного произведения относительно сложения векторов");
                        Console.WriteLine("4 - Свойство двойного векторного произведения");
                        Console.Write("Введите пункт меню: ");
                        sel1 = int.Parse(Console.ReadLine());
                        Console.WriteLine();
                        switch (sel1)
                        {
                            case 1:
                                /* [va*vb]=-[vb*va] */
                                Console.WriteLine("Проверка антикоммутативности векторного произведения");
                                Console.WriteLine("[va*vb]=-[vb*va]");
                                Console.Write("ax: ");
                                st = Console.ReadLine();
                                ax = double.Parse(st);
                                Console.Write("ay: ");
                                st = Console.ReadLine();
                                ay = double.Parse(st);
                                Console.Write("az: ");
                                st = Console.ReadLine();
                                az = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("bx: ");
                                st = Console.ReadLine();
                                bx = double.Parse(st);
                                Console.Write("by: ");
                                st = Console.ReadLine();
                                by = double.Parse(st);
                                Console.Write("bz: ");
                                st = Console.ReadLine();
                                bz = double.Parse(st);
                                Console.WriteLine();
                                VectorMultiply(ax, ay, az, bx, by, bz, out rx1, out ry1, out rz1);
                                VectorMultiply(bx, by, bz, ax, ay, az, out vx1, out vy1, out vz1);
                                OpposVect(vx1, vy1, vz1, out rx2, out ry2, out rz2);
                                Console.WriteLine($" Левая часть равенства: rx1={rx1}   ry1={ry1}   rz1={rz1}");
                                Console.WriteLine($"Правая часть равенства: rx2={rx2}   ry2={ry2}   rz2={rz2}");
                                break;
                            case 2:
                                /* [(k1*va)*vb]=k1*[va*vb] */
                                Console.WriteLine("Проверка ассоциативности векторного произведения совместно с умножением вектора на число");
                                Console.WriteLine("[(k1*va)*vb]=k1*[va*vb]");
                                Console.Write("k1: ");
                                st = Console.ReadLine();
                                k1 = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("ax: ");
                                st = Console.ReadLine();
                                ax = double.Parse(st);
                                Console.Write("ay: ");
                                st = Console.ReadLine();
                                ay = double.Parse(st);
                                Console.Write("az: ");
                                st = Console.ReadLine();
                                az = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("bx: ");
                                st = Console.ReadLine();
                                bx = double.Parse(st);
                                Console.Write("by: ");
                                st = Console.ReadLine();
                                by = double.Parse(st);
                                Console.Write("bz: ");
                                st = Console.ReadLine();
                                bz = double.Parse(st);
                                Console.WriteLine();
                                MultVectOnNum(ax, ay, az, k1, out vx1, out vy1, out vz1);
                                VectorMultiply(vx1, vy1, vz1, bx, by, bz, out rx1, out ry1, out rz1);
                                VectorMultiply(ax, ay, az, bx, by, bz, out vx1, out vy1, out vz1);
                                MultVectOnNum(vx1, vy1, vz1, k1, out rx2, out ry2, out rz2);
                                Console.WriteLine($" Левая часть равенства: rx1={rx1}   ry1={ry1}   rz1={rz1}");
                                Console.WriteLine($"Правая часть равенства: rx2={rx2}   ry2={ry2}   rz2={rz2}");
                                break;
                            case 3:
                                /* [(va+vb)*vc]=[va*vc]+[vb*vc] */
                                Console.WriteLine("Проверка дистрибутивности векторного произведения относительно сложения векторов");
                                Console.WriteLine("[(va+vb)*vc]=[va*vc]+[vb*vc]");
                                Console.Write("ax: ");
                                st = Console.ReadLine();
                                ax = double.Parse(st);
                                Console.Write("ay: ");
                                st = Console.ReadLine();
                                ay = double.Parse(st);
                                Console.Write("az: ");
                                st = Console.ReadLine();
                                az = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("bx: ");
                                st = Console.ReadLine();
                                bx = double.Parse(st);
                                Console.Write("by: ");
                                st = Console.ReadLine();
                                by = double.Parse(st);
                                Console.Write("bz: ");
                                st = Console.ReadLine();
                                bz = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("cx: ");
                                st = Console.ReadLine();
                                cx = double.Parse(st);
                                Console.Write("cy: ");
                                st = Console.ReadLine();
                                cy = double.Parse(st);
                                Console.Write("cz: ");
                                st = Console.ReadLine();
                                cz = double.Parse(st);
                                Console.WriteLine();
                                SumVect(ax, ay, az, bx, by, bz, out vx1, out vy1, out vz1);
                                VectorMultiply(vx1, vy1, vz1, cx, cy, cz, out rx1, out ry1, out rz1);
                                VectorMultiply(ax, ay, az, cx, cy, cz, out vx1, out vy1, out vz1);
                                VectorMultiply(bx, by, bz, cx, cy, cz, out vx2, out vy2, out vz2);
                                SumVect(vx1, vy1, vz1, vx2, vy2, vz2, out rx2, out ry2, out rz2);
                                Console.WriteLine($" Левая часть равенства: rx1={rx1}   ry1={ry1}   rz1={rz1}");
                                Console.WriteLine($"Правая часть равенства: rx2={rx2}   ry2={ry2}   rz2={rz2}");
                                break;
                            case 4:
                                /* [va*[vb*vc]]=(va*vc)*vb-(va*vb)*vc */
                                Console.WriteLine("Проверка свойства двойного векторного произведения");
                                Console.WriteLine("[va*[vb*vc]]=(va*vc)*vb-(va*vb)*vc");
                                Console.Write("ax: ");
                                st = Console.ReadLine();
                                ax = double.Parse(st);
                                Console.Write("ay: ");
                                st = Console.ReadLine();
                                ay = double.Parse(st);
                                Console.Write("az: ");
                                st = Console.ReadLine();
                                az = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("bx: ");
                                st = Console.ReadLine();
                                bx = double.Parse(st);
                                Console.Write("by: ");
                                st = Console.ReadLine();
                                by = double.Parse(st);
                                Console.Write("bz: ");
                                st = Console.ReadLine();
                                bz = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("cx: ");
                                st = Console.ReadLine();
                                cx = double.Parse(st);
                                Console.Write("cy: ");
                                st = Console.ReadLine();
                                cy = double.Parse(st);
                                Console.Write("cz: ");
                                st = Console.ReadLine();
                                cz = double.Parse(st);
                                Console.WriteLine();
                                VectorMultiply(bx, by, bz, cx, cy, cz, out vx1, out vy1, out vz1);
                                VectorMultiply(ax, ay, az, vx1, vy1, vz1, out rx1, out ry1, out rz1);
                                ScalarMultiply(ax, ay, az, cx, cy, cz, out k1);
                                ScalarMultiply(ax, ay, az, bx, by, bz, out k2);
                                MultVectOnNum(bx, by, bz, k1, out vx1, out vy1, out vz1);
                                MultVectOnNum(cx, cy, cz, k2, out vx2, out vy2, out vz2);
                                SubstructVect(vx1, vy1, vz1, vx2, vy2, vz2, out rx2, out ry2, out rz2);
                                Console.WriteLine($" Левая часть равенства: rx1={rx1}   ry1={ry1}   rz1={rz1}");
                                Console.WriteLine($"Правая часть равенства: rx2={rx2}   ry2={ry2}   rz2={rz2}");
                                break;
                            default:
                                Console.WriteLine("Введён неверный пункт меню!");
                                break;
                        }
                        break;
                        /* смешанное произведение */
                    case 4:
                        Console.WriteLine("                    Меню свойства смешанного произведения векторов");
                        Console.WriteLine("1 - Ассоциативность смешанного произведения совместно с умножением вектора на число");
                        Console.WriteLine("2 - Дистрибутивность смешанного произведения относительно сложения векторов");
                        Console.WriteLine("3 - Правило циклической перестановки для смешанного произведения векторов");
                        Console.WriteLine("4 - Правило перестановки двух множителей для смешанного произведения векторов");
                        Console.Write("Введите пункт меню: ");
                        sel1 = int.Parse(Console.ReadLine());
                        Console.WriteLine();
                        switch (sel1)
                        {
                            case 1:
                                /* (k1*va,vb,vc)=k1*(va,vb,vc) */
                                Console.WriteLine("Проверка ассоциативности смешанного произведения совместно с умножением вектора на число");
                                Console.WriteLine("(k1*va,vb,vc)=k1*(va,vb,vc)");
                                Console.Write("k1: ");
                                st = Console.ReadLine();
                                k1 = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("ax: ");
                                st = Console.ReadLine();
                                ax = double.Parse(st);
                                Console.Write("ay: ");
                                st = Console.ReadLine();
                                ay = double.Parse(st);
                                Console.Write("az: ");
                                st = Console.ReadLine();
                                az = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("bx: ");
                                st = Console.ReadLine();
                                bx = double.Parse(st);
                                Console.Write("by: ");
                                st = Console.ReadLine();
                                by = double.Parse(st);
                                Console.Write("bz: ");
                                st = Console.ReadLine();
                                bz = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("cx: ");
                                st = Console.ReadLine();
                                cx = double.Parse(st);
                                Console.Write("cy: ");
                                st = Console.ReadLine();
                                cy = double.Parse(st);
                                Console.Write("cz: ");
                                st = Console.ReadLine();
                                cz = double.Parse(st);
                                Console.WriteLine();
                                MultVectOnNum(ax, ay, az, k1, out vx1, out vy1, out vz1);
                                MixedMultVect(vx1, vy1, vz1, bx, by, bz, cx, cy, cz, out r1);
                                MixedMultVect(ax, ay, az, bx, by, bz, cx, cy, cz, out k2);
                                r2 = k1 * k2;
                                Console.WriteLine($" Левая часть равенства: r1={r1}");
                                Console.WriteLine($"Правая часть равенства: r2={r2}");
                                break;
                            case 2:
                                /* (va+vb,vc,vd)=(va,vc,vd)+(vb,vc,vd) */
                                Console.WriteLine("Проверка дистрибутивности смешанного произведения относительно сложения векторов");
                                Console.WriteLine("(va+vb,vc,vd)=(va,vc,vd)+(vb,vc,vd)");
                                Console.Write("ax: ");
                                st = Console.ReadLine();
                                ax = double.Parse(st);
                                Console.Write("ay: ");
                                st = Console.ReadLine();
                                ay = double.Parse(st);
                                Console.Write("az: ");
                                st = Console.ReadLine();
                                az = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("bx: ");
                                st = Console.ReadLine();
                                bx = double.Parse(st);
                                Console.Write("by: ");
                                st = Console.ReadLine();
                                by = double.Parse(st);
                                Console.Write("bz: ");
                                st = Console.ReadLine();
                                bz = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("cx: ");
                                st = Console.ReadLine();
                                cx = double.Parse(st);
                                Console.Write("cy: ");
                                st = Console.ReadLine();
                                cy = double.Parse(st);
                                Console.Write("cz: ");
                                st = Console.ReadLine();
                                cz = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("dx: ");
                                st = Console.ReadLine();
                                dx = double.Parse(st);
                                Console.Write("dy: ");
                                st = Console.ReadLine();
                                dy = double.Parse(st);
                                Console.Write("dz: ");
                                st = Console.ReadLine();
                                dz = double.Parse(st);
                                Console.WriteLine();
                                SumVect(ax, ay, az, bx, by, bz, out vx1, out vy1, out vz1);
                                MixedMultVect(vx1, vy1, vz1, cx, cy, cz, dx, dy, dz, out r1);
                                MixedMultVect(ax, ay, az, cx, cy, cz, dx, dy, dz, out k1);
                                MixedMultVect(bx, by, bz, cx, cy, cz, dx, dy, dz, out k2);
                                r2 = k1 + k2;
                                Console.WriteLine($" Левая часть равенства: r1={r1}");
                                Console.WriteLine($"Правая часть равенства: r2={r2}");
                                break;
                            case 3:
                                /* (va,vb,vc)=(vb,vc,va) */
                                Console.WriteLine("Проверка правила циклической перестановки для смешанного произведения векторов");
                                Console.WriteLine("(va,vb,vc)=(vb,vc,va)");
                                Console.Write("ax: ");
                                st = Console.ReadLine();
                                ax = double.Parse(st);
                                Console.Write("ay: ");
                                st = Console.ReadLine();
                                ay = double.Parse(st);
                                Console.Write("az: ");
                                st = Console.ReadLine();
                                az = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("bx: ");
                                st = Console.ReadLine();
                                bx = double.Parse(st);
                                Console.Write("by: ");
                                st = Console.ReadLine();
                                by = double.Parse(st);
                                Console.Write("bz: ");
                                st = Console.ReadLine();
                                bz = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("cx: ");
                                st = Console.ReadLine();
                                cx = double.Parse(st);
                                Console.Write("cy: ");
                                st = Console.ReadLine();
                                cy = double.Parse(st);
                                Console.Write("cz: ");
                                st = Console.ReadLine();
                                cz = double.Parse(st);
                                Console.WriteLine();
                                MixedMultVect(ax, ay, az, bx, by, bz, cx, cy, cz, out r1);
                                MixedMultVect(bx, by, bz, cx, cy, cz, ax, ay, az, out r2);
                                Console.WriteLine($" Левая часть равенства: r1={r1}");
                                Console.WriteLine($"Правая часть равенства: r2={r2}");
                                break;
                            case 4:
                                /* (va,vb,vc)=-(va,vc,vb) */
                                Console.WriteLine("Проверка правила перестановки двух множителей для смешанного произведения векторов");
                                Console.WriteLine("(va,vb,vc)=-(va,vc,vb)");
                                Console.Write("ax: ");
                                st = Console.ReadLine();
                                ax = double.Parse(st);
                                Console.Write("ay: ");
                                st = Console.ReadLine();
                                ay = double.Parse(st);
                                Console.Write("az: ");
                                st = Console.ReadLine();
                                az = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("bx: ");
                                st = Console.ReadLine();
                                bx = double.Parse(st);
                                Console.Write("by: ");
                                st = Console.ReadLine();
                                by = double.Parse(st);
                                Console.Write("bz: ");
                                st = Console.ReadLine();
                                bz = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("cx: ");
                                st = Console.ReadLine();
                                cx = double.Parse(st);
                                Console.Write("cy: ");
                                st = Console.ReadLine();
                                cy = double.Parse(st);
                                Console.Write("cz: ");
                                st = Console.ReadLine();
                                cz = double.Parse(st);
                                Console.WriteLine();
                                MixedMultVect(ax, ay, az, bx, by, bz, cx, cy, cz, out r1);
                                MixedMultVect(ax, ay, az, cx, cy, cz, bx, by, bz, out k1);
                                r2 = -k1;
                                Console.WriteLine($" Левая часть равенства: r1={r1}");
                                Console.WriteLine($"Правая часть равенства: r2={r2}");
                                break;
                            default:
                                Console.WriteLine("Введён неверный пункт меню!");
                                break;
                        }
                        break;
                        /* вычисления с векторами */
                    case 5:
                        Console.WriteLine("     Меню вычисления с векторами");
                        Console.WriteLine("1 - Вычисление суммы двух векторов");
                        Console.WriteLine("2 - Вычисление разности двух векторов");
                        Console.WriteLine("3 - Вычисление противоположного вектора");
                        Console.WriteLine("4 - Вычисление произведения вектора на число");
                        Console.WriteLine("5 - Вычисление скалярного произведения векторов");
                        Console.WriteLine("6 - Вычисление векторного произведения векторов");
                        Console.WriteLine("7 - Вычисление смешанного произведения векторов");
                        Console.WriteLine("8 - Вычисление длины вектора");
                        Console.WriteLine("9 - Вычисление угла между двумя векторами");
                        Console.Write("Введите пункт меню: ");
                        sel1 = int.Parse(Console.ReadLine());
                        Console.WriteLine();
                        switch (sel1)
                        {
                            case 1:
                                /* Вычисление суммы двух векторов */
                                Console.WriteLine("Введите координаты двух векторов");
                                Console.Write("ax: ");
                                st = Console.ReadLine();
                                ax = double.Parse(st);
                                Console.Write("ay: ");
                                st = Console.ReadLine();
                                ay = double.Parse(st);
                                Console.Write("az: ");
                                st = Console.ReadLine();
                                az = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("bx: ");
                                st = Console.ReadLine();
                                bx = double.Parse(st);
                                Console.Write("by: ");
                                st = Console.ReadLine();
                                by = double.Parse(st);
                                Console.Write("bz: ");
                                st = Console.ReadLine();
                                bz = double.Parse(st);
                                Console.WriteLine();
                                SumVect(ax, ay, az, bx, by, bz, out rx, out ry, out rz);
                                Console.WriteLine($"Результат: rx={rx}   ry={ry}   rz={rz}");
                                break;
                            case 2:
                                /* Вычисление разности двух векторов */
                                Console.WriteLine("Введите координаты двух векторов");
                                Console.Write("ax: ");
                                st = Console.ReadLine();
                                ax = double.Parse(st);
                                Console.Write("ay: ");
                                st = Console.ReadLine();
                                ay = double.Parse(st);
                                Console.Write("az: ");
                                st = Console.ReadLine();
                                az = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("bx: ");
                                st = Console.ReadLine();
                                bx = double.Parse(st);
                                Console.Write("by: ");
                                st = Console.ReadLine();
                                by = double.Parse(st);
                                Console.Write("bz: ");
                                st = Console.ReadLine();
                                bz = double.Parse(st);
                                Console.WriteLine();
                                SubstructVect(ax, ay, az, bx, by, bz, out rx, out ry, out rz);
                                Console.WriteLine($"Результат: rx={rx}   ry={ry}   rz={rz}");
                                break;
                            case 3:
                                /* Вычисление противоположного вектора */
                                Console.WriteLine("Введите координаты вектора");
                                Console.Write("ax: ");
                                st = Console.ReadLine();
                                ax = double.Parse(st);
                                Console.Write("ay: ");
                                st = Console.ReadLine();
                                ay = double.Parse(st);
                                Console.Write("az: ");
                                st = Console.ReadLine();
                                az = double.Parse(st);
                                Console.WriteLine();
                                OpposVect(ax, ay, az, out rx, out ry, out rz);
                                Console.WriteLine($"Результат: rx={rx}   ry={ry}   rz={rz}");
                                break;
                            case 4:
                                /* Вычисление произведения вектора на число */
                                Console.WriteLine("Введите число k и координаты вектора a");
                                Console.Write("k: ");
                                st = Console.ReadLine();
                                k = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("ax: ");
                                st = Console.ReadLine();
                                ax = double.Parse(st);
                                Console.Write("ay: ");
                                st = Console.ReadLine();
                                ay = double.Parse(st);
                                Console.Write("az: ");
                                st = Console.ReadLine();
                                az = double.Parse(st);
                                Console.WriteLine();
                                MultVectOnNum(ax, ay, az, k,out rx, out ry, out rz);
                                Console.WriteLine($"Результат: rx={rx}   ry={ry}   rz={rz}");
                                break;
                            case 5:
                                /* Вычисление скалярного произведения векторов */
                                Console.WriteLine("Введите координаты двух векторов");
                                Console.Write("ax: ");
                                st = Console.ReadLine();
                                ax = double.Parse(st);
                                Console.Write("ay: ");
                                st = Console.ReadLine();
                                ay = double.Parse(st);
                                Console.Write("az: ");
                                st = Console.ReadLine();
                                az = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("bx: ");
                                st = Console.ReadLine();
                                bx = double.Parse(st);
                                Console.Write("by: ");
                                st = Console.ReadLine();
                                by = double.Parse(st);
                                Console.Write("bz: ");
                                st = Console.ReadLine();
                                bz = double.Parse(st);
                                Console.WriteLine();
                                ScalarMultiply(ax, ay, az, bx, by, bz, out r);
                                Console.WriteLine($"Результат: r={r}");
                                break;
                            case 6:
                                /* Вычисление векторного произведения векторов */
                                Console.WriteLine("Введите координаты двух векторов");
                                Console.Write("ax: ");
                                st = Console.ReadLine();
                                ax = double.Parse(st);
                                Console.Write("ay: ");
                                st = Console.ReadLine();
                                ay = double.Parse(st);
                                Console.Write("az: ");
                                st = Console.ReadLine();
                                az = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("bx: ");
                                st = Console.ReadLine();
                                bx = double.Parse(st);
                                Console.Write("by: ");
                                st = Console.ReadLine();
                                by = double.Parse(st);
                                Console.Write("bz: ");
                                st = Console.ReadLine();
                                bz = double.Parse(st);
                                Console.WriteLine();
                                VectorMultiply(ax, ay, az, bx, by, bz, out rx, out ry, out rz);
                                Console.WriteLine($"Результат: rx={rx}   ry={ry}   rz={rz}");
                                break;
                            case 7:
                                /* Вычисление смешанного произведения векторов */
                                Console.WriteLine("Введите координаты трёх векторов");
                                Console.Write("ax: ");
                                st = Console.ReadLine();
                                ax = double.Parse(st);
                                Console.Write("ay: ");
                                st = Console.ReadLine();
                                ay = double.Parse(st);
                                Console.Write("az: ");
                                st = Console.ReadLine();
                                az = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("bx: ");
                                st = Console.ReadLine();
                                bx = double.Parse(st);
                                Console.Write("by: ");
                                st = Console.ReadLine();
                                by = double.Parse(st);
                                Console.Write("bz: ");
                                st = Console.ReadLine();
                                bz = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("cx: ");
                                st = Console.ReadLine();
                                cx = double.Parse(st);
                                Console.Write("cy: ");
                                st = Console.ReadLine();
                                cy = double.Parse(st);
                                Console.Write("cz: ");
                                st = Console.ReadLine();
                                cz = double.Parse(st);
                                Console.WriteLine();
                                MixedMultVect(ax, ay, az, bx, by, bz, cx, cy, cz, out r);
                                Console.WriteLine($"Результат: r={r}");
                                break;
                            case 8:
                                /* Вычисление длины вектора */
                                Console.WriteLine("Введите координаты вектора");
                                Console.Write("ax: ");
                                st = Console.ReadLine();
                                ax = double.Parse(st);
                                Console.Write("ay: ");
                                st = Console.ReadLine();
                                ay = double.Parse(st);
                                Console.Write("az: ");
                                st = Console.ReadLine();
                                az = double.Parse(st);
                                Console.WriteLine();
                                VectorLength(ax, ay, az, out r);
                                Console.WriteLine($"Результат: r={r}");
                                break;
                            case 9:
                                /* Вычисление угла между двумя векторами */
                                Console.WriteLine("Введите координаты двух векторов");
                                Console.Write("ax: ");
                                st = Console.ReadLine();
                                ax = double.Parse(st);
                                Console.Write("ay: ");
                                st = Console.ReadLine();
                                ay = double.Parse(st);
                                Console.Write("az: ");
                                st = Console.ReadLine();
                                az = double.Parse(st);
                                Console.WriteLine();
                                Console.Write("bx: ");
                                st = Console.ReadLine();
                                bx = double.Parse(st);
                                Console.Write("by: ");
                                st = Console.ReadLine();
                                by = double.Parse(st);
                                Console.Write("bz: ");
                                st = Console.ReadLine();
                                bz = double.Parse(st);
                                Console.WriteLine();
                                AngleVect(ax, ay, az, bx, by, bz, out r);
                                Console.WriteLine($"Результат: r={r}");
                                break;
                            default:
                                Console.WriteLine("Введён неверный пункт меню!");
                                break;
                        }
                        break;
                    case 6:
                        return;
                    default:
                        Console.WriteLine("Введён неверный пункт меню!");
                        break;
                }
                Console.WriteLine();
            }
        }
    }
}
