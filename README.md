## Mesh Simplify

#### 使用

- 编译采用了Makefile形式，在Makefile中修改成指定编译器并键入命令`make`即可编译
- 用法：
  - `simplify <input> <output> <ratio> (<t>)`
  - 其中`input`为输入的OBJ文件的路径，`output`为输出的OBJ文件的路径
  - `ratio`为简化比，即最后面片数对原数量的比
  - 最后一个参数`t`为可选项（默认为0.01），表示在选择Pair过程中两点距离的上限，和论文略有出入，在报告中有说明
- 本地环境
  - macOS 10.14.5
  - Homebrew GCC 9.1.0
  - 2.8 GHz Intel Core i7, 16GB RAM

#### 可变参数

- mesh.hpp代码中有`DFS_T_SEARCH`、`FLIP_COST`和`SHARP_COST`三个宏，注释可关闭其作用
  - 其中若`DFS_T_SEARCH`开启，则在选择Pair过程中会执行DFS，会根据当前距离进行近似剪枝，如果注释会严格计算（也做了一些优化，但还是很慢），DFS的速度非常快，而且二者效果相当
  - 如果`FLIP_COST`开启，则在选择的Pair造成其他面片翻转（法向量变化夹角余弦的绝对值小于0.2时认为翻转）时为代价函数乘一个惩罚系数
  - 如果`SHARP_COST`开启，则在尖锐边缘（夹角小于45度）的Pair会乘一个惩罚系数，这样可以尽量保持细节

#### 备注

- 具体算法说明和细节请参见报告