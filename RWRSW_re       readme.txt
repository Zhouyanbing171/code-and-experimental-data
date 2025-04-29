数据文件解释：
    原始数据文件：
        ec.tab              ec分组数据
        gene.gaf            本体-基因的注释文件
        net.txt             原始基因网络数据（用于做netHashMap）
        onto.obo            基因本体术语总库
    中间数据文件：
        geneList.txt        需要用的基因名称列表  来自于net.txt  此文件经过了进一步的筛选处理
        netMatrix.txt       原始基因网络数据矩阵  来自于net.txt
        pMatrix.txt         归一化的概率转移矩阵  来自于RandWalk类的计算
        child.txt           后代术语集合
        son.txt             直接后代术语集合
    结果数据文件：
        result0.9.mat       启动参数为0.9的随机游走矩阵计算结果
        lfc.txt             算法在各个ec分组下的得分计算结果

代码文件解释：
    Annotation  注释信息数据模型
    FunctionNet 本体术语网络数据模型
    Term        本体术语数据模型
    Reader      特殊类型文件读取 辅助工具类
    Matrix      矩阵 辅助工具类
    RandWalk    随机游走算法 类
    Calculator  核心计算类，根据算法流程图划分为4个部分，并向外提供必要的三个接口
    Main        算法启动类，组合三个接口，串连整个算法
