#pragma once
#include<fstream>
template<class T>
class file
{
    file(){}
    public:
    std::fstream File;
    ~file() = default; // 确保文件流在对象销毁时关闭
    file(const file&) = delete; // 禁止拷贝构造
    file& operator=(const file&) = delete; // 禁止拷贝赋值
    file(file&&) = delete; // 禁止移动构造
    file& operator=(file&&) = delete; // 禁止移动赋值
    // 静态成员函数用于创建单例对象
    static file makeObj();
    void write(char*,T&);
    void read(char*,T&);
};
template<class T>
file<T> file<T>::makeObj()
{
    static file f;
    return f;
}
template<class T>
void file<T>::write(char *name,T &obj)
{
    File.open(name,std::ios::binary|std::ios::out);
    File.write(reinterpret_cast<char*>(&obj),sizeof(obj));
    File.close();
}
template<class T>
void file<T>::read(char *name,T &obj)
{
    File.open(name,std::ios::binary|std::ios::in);
    File.read(reinterpret_cast<char*>(&obj),sizeof(obj));
    File.close();
}