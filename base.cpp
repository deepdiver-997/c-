#pragma once
#include"base.hpp"
#include<iostream>
#include<iomanip>
using namespace std;
template<class T>
matrix<T>::matrix(int r,int c):row(r),col(c),factors(1),detvalue(0),Rank(0)
{
    p.resize(row, std::vector<T>(col));
}
template<class T>
matrix<T>::matrix(const matrix &m):row(m.getrow()),col(m.getcol()),detvalue(0),factors(1),Rank(0){
    p.resize(row, std::vector<T>(col));
    for(int i=0;i<row;++i){
        for(int j=0;j<col;++j)
            p[i][j]=m[i][j];
    }
}
template<class T>
matrix<T>::matrix(matrix &&m)
:row(m.row),col(m.col),detvalue(0),factors(1),Rank(0),p(std::move(m.p))
{
    // p已经通过移动构造获取了m.p的所有权
}
template<class T>
matrix<T>::matrix(int r,int c,std::vector<std::vector<T>> &vec)
:row(r),col(c),p(std::move(vec)),factors(1),detvalue(0),Rank(0)
{
    // vec已经通过移动构造获取了所有权
}
template<class T>
void matrix<T>::input(){
        for(int i=0;i<row;++i){
        for(int j=0;j<col;++j)
        std::cin>>p[i][j];
    }
}
template<class T>
void matrix<T>::show(){
    for(int i=0;i<row;++i){
        for(int j=0;j<col;++j){
            //if(p[i][j]>=0)
            //std::cout<<fixed<<setprecision(4)<<abs(p[i][j])<<" ";
            //else
            std::cout<<fixed<<setprecision(4)<<(p[i][j])<<" ";
        }
        std::cout<<std::endl;
    }
}
template<class T>
void matrix<T>::rank(){
    Rank=0;matrix<T> q(*this);
    q.simpleline();
    for(int i=0;i<row;++i){
        for(int j=0;j<col;++j)
        if(q[i][j]!=0){
            ++Rank;
            break;
        }
    }
}
template<class T>
void matrix<T>::simpleline(){
    for(int j=0;j<getcol();++j){
        for(int i=j;i<getrow();++i){
            if((*this)[i][j]!=0){
                if(i!=j)
                first_prim_matrix(i,j);
                break;
            }
        }
    }
    for(int j=0;j<getcol()-1;++j){
        for(int i=j+1;i<getrow();++i){
            if((*this)[i][j]!=0){
                T a=(*this)[i][j],b=(*this)[j][j];
                second_prim_matrix(i,b);
                third_prim_matrix(i,j,-a);
            }
        }
    }
}
template<class T>
void matrix<T>::make_matrix(int r,int c){
    row=r,col=c;
    p.resize(row, std::vector<T>(col));
}
template<class T>
std::vector<T>& matrix<T>::operator [] (int i){
    return p[i];
}
template<class T>
const std::vector<T>& matrix<T>::operator [] (int i) const{
    return p[i];
}
template<class T>
bool matrix<T>::operator == (const matrix &m)const{
    if(row!=m.getrow()||col!=m.getcol())
    return false;
    else{
        bool b=true;
        for(int i=0;i<row;++i){
            for(int j=0;j<col;++j)
            if(p[i][j]!=m[i][j]){
                b=false;
                goto S;
            }
        }
        S:return b;
    }
}
template<class T>
matrix<T> matrix<T>::operator + (const matrix<T> &m)
{
    if(row == m.row && col == m.col)
    {
        matrix<T> temp(*this);
        for(int i=0;i<m.row;++i)
        for(int j=0;j<m.col;++j)
        temp[i][j] += m[i][j];
        return temp;
    }
    else
    {
        throw std::runtime_error("Matrix dimensions do not match for addition");
    }
}
template<class T>
matrix<T> matrix<T>::operator - (const matrix<T> &m)
{
    if(row == m.row && col == m.col)
    {
        matrix<T> temp(*this);
        for(int i=0;i<m.row;++i)
        for(int j=0;j<m.col;++j)
        temp[i][j] -= m[i][j];
        return temp;
    }
    else
    {
        throw std::runtime_error("Matrix dimensions do not match for subtraction");
    }
}
template<class T>
matrix<T> matrix<T>::operator * (const matrix<T> &m)
{
    if(col != m.row)
    {
        throw std::runtime_error("Matrix dimensions do not match for multiplication");
    }
    else
    {
        matrix<T> temp(row,m.getcol());
        for(int i=0;i<row;++i)
        {
        for(int j=0;j<m.getcol();++j)
            for(int k=0;k<col;++k)
            temp[i][j] += this->p[i][k] * m[k][j];
        }
        return temp;
    }
}
// 友元函数已在头文件中内联实现
template<class T>
void matrix<T>::transposition(){
    std::vector<std::vector<T>> temp(col, std::vector<T>(row));
    for(int i=0;i<row;++i)
        for(int j=0;j<col;++j)
            temp[j][i]=p[i][j];
    
    // 交换行列数
    std::swap(row, col);
    p = std::move(temp);
}
template<class T>
void matrix<T>::first_prim_matrix(int i,int j){
    std::swap(p[i], p[j]);
    factors*=-1;
}
template<class T>
void matrix<T>::second_prim_matrix(int i,T k){
    for(int j=0;j<col;++j)
        p[i][j]*=k;
    factors*=k;
}
template<class T>
void matrix<T>::third_prim_matrix(int i,int j,T k){
    for(int t=0;t<col;++t)
        p[i][t]+=p[j][t]*k;
}
template<class T>
T matrix<T>::detcalculate_1(int n, const std::vector<std::vector<T>>& mat){
    if(n==2)
        return mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0];
    else{
        std::vector<std::vector<T>> submat(n-1, std::vector<T>(n-1));
        T sum=0;
        
        for(int i=0;i<n;++i){
            for(int j=1;j<n;++j){
                for(int k=0;k<n-1;++k){
                    if(i>k)
                        submat[j-1][k]=mat[j][k];
                    else
                        submat[j-1][k]=mat[j][k+1];
                }
            }
            sum+=pow(-1,i)*mat[0][i]*detcalculate_1(n-1,submat);
        }
        return sum;
    }
}
template<class T>
void matrix<T>::detcalculate_1(){
    detvalue=detcalculate_1(row,p);
}
template<class T>
void matrix<T>::detcalculate_2(){
    if(row==col){
        T val=1;
        matrix t(*this);
        t.simpleline();
        for(int i=0;i<row;++i)
            val*=t.p[i][i];
        val=val/t.factors;
        detvalue=val;
    }
}
template<class T>
double matrix<T>::row_vector_norm(int r,double P){
    double d=0;
    for(int i=0;i<col;++i)
        d+=pow(p[r][i],P);
    d=pow(d,1/P);
    if(d>LLONG_MAX){
        cout<<"This is an infinite norm\n";
        double max=abs(p[r][0]);
        for(int i=1;i<col;++i)
            if(max<abs(p[r][i]))
                max=abs(p[r][i]);
        return max;
    }
    return d;
}
template<class T>
double matrix<T>::col_vector_norm(int c,double P){
    double d=0;
    for(int i=0;i<row;++i)  // 修正：应该是row而不是col
        d+=pow(p[i][c],P);
    d=pow(d,1/P);
    if(d>LLONG_MAX){
        cout<<"This is an infinite norm\n";
        double max=abs(p[0][c]);
        for(int i=1;i<row;++i)  // 修正：应该是row而不是col
            if(max<abs(p[i][c]))
                max=abs(p[i][c]);
        return max;
    }
    return d;
}
template<class T>
double matrix<T>::Frobenius_norm(){
    double f=0;
    for(int i=0;i<row;++i)
        for(int j=0;j<col;++j)
            f+=pow(p[i][j],2);
    f=pow(f,0.5);
    return f;
}
template<class T>
double matrix<T>::matrix_row_sum_norm(){
    double max=abs(row_vector_norm(0,1));
    for(int i=1;i<row;++i)
        if(max<abs(row_vector_norm(i,1)))
            max=abs(row_vector_norm(i,1));
    return max;
}
template<class T>
double matrix<T>::matrix_col_sum_norm(){
    double max=abs(col_vector_norm(0,1));
    for(int i=1;i<col;++i)
        if(max<abs(col_vector_norm(i,1)))
            max=abs(col_vector_norm(i,1));
    return max;
}