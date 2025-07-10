#ifndef DOU
#define DOU
#include<iostream>
#include<memory>
#include<vector>
#include"base.hpp"
class matrix_double : public matrix<double>{
    protected:
    std::unique_ptr<double[]> x; // 用来存方程组的特解
    std::vector<std::vector<double>> X; // 用来存基础解系
    std::unique_ptr<double[]> eigen; // 用来存特征值
    double cond, spectral_norm; // cond是条件数
    std::unique_ptr<int[]> arr; // 用来存对应基础解系的系数
    int only, vary, e; // only,vary,e用来检测对应函数是否运行,以便在程序结束时正确释放空间
    //matrix_double **characteristic_vector;//存特征向量
    public:
    matrix_double():matrix<double>(),only(0),vary(0),e(0){}
    matrix_double(const matrix_double &m):matrix<double>(m),only(0),vary(0),e(0){}
    matrix_double(matrix_double &&m):matrix<double>(std::move(m)),only(0),vary(0),e(0){}
    matrix_double(int r,int c):matrix<double>(r,c),only(0),vary(0),e(0){}
    ~matrix_double() = default; // 智能指针会自动管理内存
    matrix_double & operator = (const matrix_double &m);
    friend matrix_double operator + (matrix_double &m,matrix_double &n);
    friend matrix_double operator - (matrix_double &m,matrix_double &n);
    friend matrix_double operator * (matrix_double &m,matrix_double &n);//通过交换循环顺序减少对矩阵m和n的访问次数来提高效率
    friend matrix_double operator * (double t,matrix_double &m);
    friend matrix_double operator * (matrix_double &m,double t);
    matrix_double ARF(int i,int j);//返回该矩阵的第i行j列的余子式
    matrix_double accompany();//返回该矩阵的伴随矩阵,如果没有伴随矩阵返回0矩阵
    matrix_double inverse();//返回该矩阵的逆,如果没有逆返回0矩阵
    void equation();//解方程组  ps:用来解方程组的矩阵必须是增广矩阵
    void equation_only();//求特解
    void varEq();//求基础解系
    void showx();//展示特解
    double getxi(int i){return x ? x[i] : 0;}
    double getX(int i,int j){return (i < X.size() && j < X[i].size()) ? X[i][j] : 0;}
    double gete(){return e;}
    double get_eigen(int i){return eigen ? eigen[i] : 0;}
    double getcond(){return cond;}
    double getspectral_norm(){return spectral_norm;}
    bool test();
    void copyx(int i,matrix_double &m);
    void showX();
    void solution();//展示通解,但是必须先运行对象的equation函数
    double inner_product(std::vector<double> &a,std::vector<double> &b,int l);//返回向量a和b的内积,l是向量维数
    double coefficient(int i,int j);//施密特正交化向量的系数
    matrix_double Gram_Schmidt();//返回以列向量为基的正交化矩阵
    void matrix_standardise();//矩阵列向量单位化
    matrix_double QRfactorization_Q();//返回QR分解的Q矩阵
    matrix_double QRfactorization_R(matrix_double &Q);//返回QR分解的R矩阵
    void eigen_QR();//QR分解法求特征值
    void showEigen();//展示特征值
    void Cond();//求条件数
    void Spectral_norm();//求谱范数
};
#endif