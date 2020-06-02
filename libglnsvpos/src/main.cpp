#include <iostream>
#include <fstream>
#include <include\libglnsvpos\glnsvpos.h>
#include <include\libglnsvpos\rungekutta.h>
#include <ctime>

using namespace std;

int main()
{
    //    14 20  2 10 13 45  0.0  .479789450765E-04  .000000000000E+00  .495000000000E+05
    //         .122349252930E+05 -.683955192566E+00  .186264514923E-08  .000000000000E+00
    //        -.154207099609E+05  .211544322968E+01  .931322574615E-09 -.700000000000E+01
    //         .162340239258E+05  .252534008026E+01 -.279396772385E-08  .000000000000E+00
    time_t start, end;
    double delt = 0.1;
    int n = (int) 12*3600/delt;
    double **koord = new double * [n];
    double *koord_matlab = new double[3];
    double max_del = 0;
    int i_max = 0;
    for (int i = 0; i < n; i++)
    {
        koord[i] = new double [6];
    }
    ofstream out;
    out.open("D:\\res_cpp.txt"); // путь сохранения результатов расчетов на С++
    time(&start);
    ifstream in("D:\\res_mat.txt"); // путь для чтения результатов матлаба
    if (!in)
    {
        cout << "File not open!" << endl;
    } else {
        cout << "File open!" << endl;
    }
    time(&start);
    glns_coord(koord);
    time(&end);
    for (int i = 0; i < n; i++)
    {
        string koord_str1 = to_string(koord[i][0]);
        string koord_str2 = to_string(koord[i][1]);
        string koord_str3 = to_string(koord[i][2]);
        out << koord_str1 << "\t" << koord_str2 << "\t" << koord_str3 << endl;
        in >> koord_matlab[0] >> koord_matlab[1] >> koord_matlab[2];
        for (int j = 0; j < 3; j++)
        {
            if (abs(koord[i][j]-koord_matlab[j]) > max_del)
            {
                max_del = abs(koord[i][j]-koord_matlab[j]);
                i_max = i;
            }
        }
        delete [] koord[i];
        koord[i] = nullptr;
    }
    in.close();
    out.close();
    delete[] koord;
    koord = nullptr;
    delete[] koord_matlab;
    koord_matlab = nullptr;
    double seconds = difftime(end, start);
    string seconds1 = to_string(seconds*1000000/n);
    cout << "Vremya raboty funkcii " << seconds1 << " s"<< endl;
    string max_del1 = to_string(max_del);
    cout << "maksimalnaya raznitca koordinat " << max_del1 << " m" << endl;
    string imax = to_string(i_max);
    cout << "Nomer otscheta s maksimalnoy raznitcey " << imax << endl;

}
