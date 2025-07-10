#ifndef MATRIX
#define MATRIX
#include<iostream>
#include<cmath>
#include<vector>
#include<memory>
template<class T>
class matrix{
    protected:
    int row,col,Rank;//ps:第0行就是最上的那行,第0列就是最左的那列
    std::vector<std::vector<T>> p;//p是存放数据的容器;detvalue是矩阵的行列式的值,非方阵都是0;factors是高斯消元产生的因子,是求解行列式的中间量
    T detvalue,factors;
    T getfactors(){return factors;}
    public:
    matrix():factors(1),detvalue(0),Rank(0){}
    matrix(const matrix&);
    matrix(matrix&&);
    matrix(int r,int c,std::vector<std::vector<T>>&);
    matrix(int r,int c);
    virtual ~matrix() = default;
    int getrow()const{return row;}
    int getcol()const{return col;}
    std::vector<std::vector<T>>& getp(){return p;}
    void rank();
    int getrank(){if(Rank==0)rank();return Rank;}
    void setrow(int r){row=r;}
    void setcol(int c){col=c;}
    void input();
    virtual void show();
    void simpleline();//化简为阶梯矩阵
    void make_matrix(int r,int c);
    std::vector<T>& operator [] (int i);
    const std::vector<T>& operator [] (int i) const;
    bool operator == (const matrix &m)const;
    matrix operator + (const matrix&);
    matrix &operator = (const matrix&);
    matrix operator - (const matrix&);
    matrix operator * (const matrix&);
    friend matrix<T> operator * (T t, const matrix<T>& m) {
        matrix<T> temp(m);
        for(int i=0; i<m.row; ++i)
            for(int j=0; j<m.col; ++j)
                temp[i][j] *= t;
        return temp;
    }
    friend matrix<T> operator * (const matrix<T>& m, T t) {
        return t * m;
    }
    void transposition();//矩阵转置
    T detcalculate_1(int n, const std::vector<std::vector<T>>& p);
    void detcalculate_1();//返回行列式的值到detvalue的函数,用的是行列式按行展开的递归算法,不推荐
    void setvalue(T v){detvalue=v;}
    T getvalue(){detcalculate_1();return detvalue;}
    void detcalculate_2();//返回行列式的值到detvalue的函数,用的是高斯消元将方阵化为上三角
    void first_prim_matrix(int i,int j);//交换矩阵的第i和j行
    void second_prim_matrix(int i,T k);//矩阵的第i行乘k
    void third_prim_matrix(int i,int j,T k);//矩阵的第i行加第j行乘k
    double row_vector_norm(int r,double p);//矩阵第r行向量的p范数
    double col_vector_norm(int c,double p);//矩阵第c列向量的p范数
    double Frobenius_norm();//矩阵的Frobenius范数
    double matrix_row_sum_norm();//矩阵行和范数
    double matrix_col_sum_norm();//矩阵列和范数
};
#endif