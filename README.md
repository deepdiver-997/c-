# C++ 矩阵计算库

这是一个功能强大的C++矩阵计算库，提供了丰富的矩阵操作和数值计算功能。该库包含基础矩阵模板类和专门用于双精度浮点数的矩阵类，支持各种矩阵运算、线性方程组求解、矩阵分解等高级数值计算功能。

## 项目结构

```
├── base.cpp    - 基础矩阵模板类实现
├── base.hpp    - 基础矩阵模板类定义
├── double.cpp  - 双精度浮点数矩阵类实现
├── double.hpp  - 双精度浮点数矩阵类定义
└── file.hpp    - 矩阵文件I/O操作类
```

## 功能特点

### 基础矩阵操作 (`matrix<T>` 类)

- 矩阵创建与初始化
- 基本矩阵运算：加法、减法、乘法
- 矩阵转置
- 行列式计算
- 矩阵秩的计算
- 初等变换：
  - 交换行
  - 行乘以常数
  - 行加上另一行的倍数
- 矩阵范数计算：
  - 行向量范数
  - 列向量范数
  - 矩阵列和范数
  - Frobenius范数

### 高级数值计算 (`matrix_double` 类)

- 伴随矩阵计算
- 矩阵求逆
- 线性方程组求解：
  - 唯一解求解
  - 通解求解（含参数解）
- Gram-Schmidt正交化
- QR分解
- 特征值计算
- 条件数计算
- 谱范数计算

### 文件操作 (`file<T>` 类)

- 矩阵对象的二进制文件读写
- 使用单例模式实现文件操作

## 使用示例

### 基本矩阵操作

```cpp
#include "base.hpp"

int main() {
    // 创建3x3矩阵
    matrix<double> m(3, 3);
    
    // 输入矩阵元素
    m.input();
    
    // 显示矩阵
    m.show();
    
    // 计算行列式
    double det = m.getvalue();
    std::cout << "Determinant: " << det << std::endl;
    
    // 计算矩阵秩
    int rank = m.getrank();
    std::cout << "Rank: " << rank << std::endl;
    
    // 矩阵转置
    m.transposition();
    std::cout << "Transposed matrix:" << std::endl;
    m.show();
    
    return 0;
}
```

### 线性方程组求解

```cpp
#include "double.hpp"

int main() {
    // 创建增广矩阵 (系数矩阵 | 常数向量)
    matrix_double m(3, 4);
    
    // 输入增广矩阵
    m.input();
    
    // 求解线性方程组
    try {
        m.equation();
        
        // 显示解
        std::cout << "Solution: ";
        m.solution();
        std::cout << std::endl;
    } catch (const std::runtime_error& e) {
        std::cout << e.what() << std::endl;
    }
    
    return 0;
}
```

### 矩阵分解和特征值计算

```cpp
#include "double.hpp"

int main() {
    // 创建矩阵
    matrix_double m(3, 3);
    m.input();
    
    // QR分解
    matrix_double Q = m.QRfactorization_Q();
    matrix_double R = m.QRfactorization_R(Q);
    
    std::cout << "Q matrix:" << std::endl;
    Q.show();
    std::cout << "R matrix:" << std::endl;
    R.show();
    
    // 计算特征值
    m.eigen_QR();
    std::cout << "Eigenvalues: ";
    m.showEigen();
    
    return 0;
}
```

### 文件操作

```cpp
#include "file.hpp"
#include "double.hpp"

int main() {
    // 创建矩阵
    matrix_double m(3, 3);
    m.input();
    
    // 创建文件操作对象
    auto fileObj = file<matrix_double>::makeObj();
    
    // 写入文件
    char filename[] = "matrix.dat";
    fileObj.write(filename, m);
    
    // 读取文件
    matrix_double m2(3, 3);
    fileObj.read(filename, m2);
    
    // 显示读取的矩阵
    m2.show();
    
    return 0;
}
```

## 编译和使用

该库使用标准C++编写，无需额外依赖。可以使用任何支持C++11或更高版本的编译器进行编译。

```bash
# 使用g++编译
g++ -std=c++11 your_program.cpp base.cpp double.cpp -o your_program

# 使用clang++编译
clang++ -std=c++11 your_program.cpp base.cpp double.cpp -o your_program
```

## 注意事项

- 该库使用模板编程，支持多种数据类型的矩阵操作
- 高级数值计算功能仅在 `matrix_double` 类中提供
- 文件I/O操作使用二进制模式，确保跨平台兼容性
- 矩阵索引从0开始

## 未来改进方向

- 添加更多矩阵分解方法（如LU分解、SVD分解）
- 优化大型矩阵的计算性能
- 添加稀疏矩阵支持
- 提供更多数值稳定性的算法实现
- 增加并行计算支持
