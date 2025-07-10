#pragma once
#include<iostream>
#include<memory>
#include<vector>
#include"double.hpp"
#include"base.cpp"
#include"base.hpp"
matrix_double & matrix_double::operator = (const matrix_double &m){
    row = m.getrow();
    col = m.getcol();
    p.resize(row, std::vector<double>(col));
    for(int i=0;i<getrow();++i){
        for(int j=0;j<getcol();++j)
        p[i][j]=m[i][j];
    }
    return *this;
}
matrix_double operator + (matrix_double &m,matrix_double &n){
    matrix_double M(m);
    if(m.getrow()==n.getrow()&&m.getcol()==n.getcol())
    for(int i=0;i<m.getrow();++i){
        for(int j=0;j<m.getcol();++j)
        M[i][j]+=n[i][j];
    }
    else{
        throw std::runtime_error("Matrix dimensions do not match for addition");
    }
    return M;
}
matrix_double operator - (matrix_double &m,matrix_double &n){
    matrix_double M(m);
    if(m.getrow()==n.getrow()&&m.getcol()==n.getcol())
    for(int i=0;i<m.getrow();++i){
        for(int j=0;j<m.getcol();++j)
        M[i][j]-=n[i][j]; // 修正：应该是减法而不是加法
    }
    else{
        throw std::runtime_error("Matrix dimensions do not match for subtraction");
    }
    return M;
}
matrix_double operator * (matrix_double &m,matrix_double &n){
    matrix_double q(m.getrow(),n.getcol());
    if(m.getcol()==n.getrow())
    for(int i=0;i<m.getrow();++i){
        for(int j=0;j<n.getcol();++j){
            for(int k=0;k<m.getcol();++k)
            q[i][j]+=m[i][k]*n[k][j];
        }
    }
    else{
        throw std::runtime_error("Matrix dimensions do not match for multiplication");
    }
    return q;
}
matrix_double operator * (matrix_double &m,double t){
    matrix_double q(m);
    for(int i=0;i<m.getrow();++i)
    for(int j=0;j<m.getcol();++j)
    q[i][j] *= t;
    return q;
}
matrix_double operator * (double t,matrix_double &m){
    matrix_double q(m);
    for(int i=0;i<m.getrow();++i)
    for(int j=0;j<m.getcol();++j)
    q[i][j] *= t;
    return q;
}
matrix_double matrix_double::ARF(int i,int j){
    matrix_double q(getrow()-1,getcol()-1);
    for(int u=0;u<getrow()-1;++u){
        for(int v=0;v<getcol()-1;++v){
            if(u<i){
                if(v<j)
                q[u][v]=(*this)[u][v];
                if(v>=j)
                q[u][v]=(*this)[u][v+1];
            }
            else{
                if(v<j)
                q[u][v]=(*this)[u+1][v];
                if(v>=j)
                q[u][v]=(*this)[u+1][v+1];
            }
        }
    }
    return std::move(q);
}
matrix_double matrix_double::accompany(){
    matrix_double q(getrow(),getcol());
    if(getrank()==getrow()&&getrow()==getcol())
    for(int i=0;i<getrow();++i){
        for(int j=0;j<getcol();++j){
            q[i][j]=pow(-1,i+j)*ARF(i,j).getvalue();
        }
    }
    else std::cout<<"Conditions Unsatisfied\n";
    return std::move(q);
}
matrix_double matrix_double::inverse(){
    matrix_double q(accompany());
    if(getrank()==0){
        std::cout<<"This matrix has no inverse matrix\n";
        return std::move(q);
    }
    double d=getvalue();
    for(int i=0;i<getrow();++i)
    q.second_prim_matrix(i,1/d);
    return std::move(q);
}
void matrix_double::equation_only(){
    x=std::make_unique<double[]>(getcol()-1);
    only=1;
    simpleline();
    for(int i=getcol()-2;i>=0;--i){
        double d=0.0;
        for(int j=getcol()-2;j>i;--j)
        d=d+x[j]*(*this)[i][j];
        x[i]=((*this)[i][getcol()-1]-d)/(*this)[i][i];
    }
}
void matrix_double::showx(){
    for(int i=0;i<getcol()-1;++i)
    cout<<x[i]<<"  ";
}
bool matrix_double::test(){
    int a=0;
    for(int i=0;i<getcol()-1;++i)
    if((*this)[getrank()][i]!=0)
    ++a;
    if(a>1)
    return false;
    else return true;
}
void matrix_double::copyx(int i,matrix_double &m){
    X[i]=std::vector<double>(getcol()-1);
    for(int j=0;j<getcol()-1;++j){
        X[i][j]=m.getxi(j);
    }
}
void matrix_double::varEq(){
    simpleline();
    vary=getcol()-getrank()-1;
    matrix_double temp(getcol()-1,getcol());
    int j;
    arr=std::make_unique<int[]>(vary);
    X=std::vector<std::vector<double>>(vary);
    for(int i=0;i<getrank();++i){
        for(j=0;j<getcol();++j){
            if((*this)[i][j]!=0)
            break;
        }
        for(int k=j;k<getcol();++k){
            temp[j][k]=(*this)[i][k];
        }
    }
    j=0;
    for(int i=0;i<getcol()-1;++i){
        if(temp[i][i]==0){
            temp[i][i]=1;
            arr[j]=i;X[j]=std::vector<double>(getcol()-getrank()-1);
            ++j;
        }
    }
    temp.equation_only();
    // temp.showx();
    x = std::make_unique<double[]>(getcol()-1);
    for(int i=0;i<getcol()-1;++i)
    x[i]=temp.x[i];
    only=1;
    temp[arr[0]][getcol()-1]=1;
    temp.equation_only();
    copyx(0,temp);
    for(int i=1;i<getcol()-getrank()-1;++i){
        temp[arr[i-1]][getcol()-1]=0;
        temp[arr[i]][getcol()-1]=1;
        temp.equation_only();
        copyx(i,temp);
    }
}
void matrix_double::showX(){
    for(int i=0;i<X.size();++i){
    for(int j=0;j<X[i].size();++j)
    cout<<X[i][j]<<" ";
    cout<<endl;}
}
void matrix_double::solution(){
    bool b=false;
    for(int i=0;i<getcol()-1;++i)
    if(x[i]==0)
    continue;
    else{
        b=true;
        break;
    }
    if(b){
        cout<<'(';
        for(int i=0;i<getcol()-2;++i)
        cout<<x[i]<<',';
        cout<<x[getcol()-2]<<')';
    }
    if(getcol()-1>getrank()&&X.size()>0){
        if(b){
            cout<<"+k"<<arr[0]<<'(';
            for(int j=0;j<getcol()-2;++j)
            cout<<X[0][j]<<',';
            cout<<X[0][getcol()-2]<<')';
        }
        else{
            cout<<"k"<<arr[0]<<'(';
            for(int j=0;j<getcol()-2;++j)
            cout<<X[0][j]<<',';
            cout<<X[0][getcol()-2]<<')';
        }
        for(int i=1;i<getcol()-getrank()-1;++i){
            cout<<"+k"<<arr[i]<<'(';
            for(int j=0;j<getcol()-2;++j)
            cout<<X[i][j]<<',';
            cout<<X[i][getcol()-2]<<')';
        }
    }
    //cout<<endl<<"x:"<<x<<"\nX:"<<X<<endl;
}
void matrix_double::equation(){
    matrix_double q(getrow(),getcol()-1);
    for(int i=0;i<getrow();++i){
        for(int j=0;j<getcol()-1;++j)
        q[i][j]=(*this)[i][j];
    }
    if(getrank()==q.getrank()&&getrank()==getcol()-1){
        equation_only();
    }
    if(getrank()==q.getrank()&&getrank()<getcol()-1){
        varEq();
    }
    if(getrank()!=q.getrank()){
        throw std::runtime_error("No solutions to this equation system");
    }
}
double matrix_double::inner_product(std::vector<double> &a,std::vector<double> &b,int l){
    double r=0;
    for(int i=0;i<l;++i)
    r+=a[i]*b[i];
    return r;
}
double matrix_double::coefficient(int i,int j){
    return inner_product((*this)[i],(*this)[j],getcol())/pow(row_vector_norm(j,2),2);
}
matrix_double matrix_double::Gram_Schmidt(){
    matrix_double a(*this);
    double d;
    for(int i=0;i<a.getrow();++i){
        matrix_double minus(getrow(),getcol());
        for(int j=0;j<i;++j){
            d=a.coefficient(i,j);//cout<<i<<j<<": "<<d<<endl;
            for(int k=0;k<minus.getcol();++k){
                minus[i][k] += d * a[j][k];
            }
            //minus.show();cout<<endl;
        }
        //minus.show();cout<<endl<<"a:\n";
        for(int j=0;j<getcol();++j)
        a[i][j]-=minus[i][j];
    }
    return std::move(a);
}
void matrix_double::matrix_standardise(){
    double d;
    for(int i=0;i<getrow();++i){
        d=row_vector_norm(i,2);
        second_prim_matrix(i,1/d);
    }
}
matrix_double matrix_double::QRfactorization_Q(){
    if(getrow()>=getcol()){
        transposition();
        matrix_double Q(Gram_Schmidt());
        Q.matrix_standardise();
        Q.transposition(),transposition();
        return std::move(Q);
    }
    else{
        matrix_double Q(getrow(),getrow());
        for(int i=0;i<getrow();++i)
        for(int j=0;j<getrow();++j)
        Q[i][j]=(*this)[i][j];
        Q.transposition();
        Q=Q.Gram_Schmidt();
        Q.transposition();
        return std::move(Q);
    }
}
matrix_double matrix_double::QRfactorization_R(matrix_double &Q){
    if(Q.getrow()>=Q.getcol()){
        matrix_double R(getcol(),getcol());
        Q.transposition(),transposition();
        for(int i=0;i<Q.getrow();++i)
        for(int j=i;j<getcol();++j)
        R[i][j]=inner_product(Q[i],(*this)[j],getcol());
        return std::move(R);
    }
    else{
        matrix_double R(getrow(),getcol());
        Q.transposition(),transposition();
        for(int i=0;i<Q.getrow();++i)
        for(int j=i;j<getcol();++j)
        R[i][j]=inner_product(Q[i],(*this)[j],getcol());
        return std::move(R);
    }
}
void matrix_double::eigen_QR(){
    matrix_double A(*this),B(*this);
    int i=0,b;
        do
        {
            matrix_double Q(A.QRfactorization_Q()),R(A.QRfactorization_R(Q));
            A=B;
            Q.transposition();
            B=R*Q;++i;
        }while (i<50);A.show();
    //for(int i=0;i<(getrow()>getcol()?getcol():getrow());++i)
    //cout<<A[i][i]<<endl;
    e=A.getrank();
    eigen=std::make_unique<double[]>(e);
    for(i=0;i<e;++i){
        eigen[i]=A[i][i];
        /*matrix_double temp(C.getrow(),C.getcol()+1);
        for(int j=0;j<C.getrow();++j)
        for(int k=0;k<C.getcol();++k)
        if(j==k)
        temp[j][k]=C[j][k]-eigen[i];
        else
        temp[j][k]=C[j][k];
        temp.show();
        //(*this)=temp;
        temp.equation();temp.showX();cout<<endl;
        characteristic_vector[i]=new matrix_double(C.getcol(),vary);
        for(int u=0;u<vary;++u)
        for(int v=0;v<C.getcol();++v)
        (*characteristic_vector[i])[v][u]=X[u][v];
        matrix_double::add(characteristic_vector[i]);
        //temp.solution();cout<<endl;*/
    }
}
    
void matrix_double::showEigen(){
    for(int i=0;i<e;++i)
    cout<<eigen[i]<<" ";
    cout<<endl;
}
void matrix_double::Cond(){
    if(getrank()==0)
    cout<<"Cond doesn't exist\n";
    else{
        double a=matrix_col_sum_norm(),b=inverse().matrix_col_sum_norm();
        cond=a*b;
    }
}
void matrix_double::Spectral_norm(){
    matrix_double A2(*this);A2.transposition();
    matrix_double A1(A2*(*this));
    A1.eigen_QR();
    double max=0;
    for(int i=0;i<A1.gete();++i)
    if(max<A1.get_eigen(i))
    max=A1.get_eigen(i);
    spectral_norm=pow(max,0.5);
}
// int main(){
//     matrix_double m(5,6);m.input();m.equation();m.showx();//cout<<m.getcond();
//     matrix_double::reduce();
//     return 0;
// }
/*
3 1 -1
1 3 -1
-1 -1 3

3 1 -1
1 3 -1
-1 -1 3
*/