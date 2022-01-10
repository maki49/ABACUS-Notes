<font face="Comic Sans MS">

## Background
在局域轨道基组下并行计算时，H，S，density matrix等(nbasis*nbasis)矩阵有两种分法：
- 按轨道（basis）分成2D块：用于对角化求解KS本征值方程
- 按格点分：用于格点积分，通过density matrix求解charge density

## Method
#### BLAS: Basic Linear Algebra Subprograms
- Level 1: vector-vector
- Level 2: matrix-vector
- Level 3: matrix-matrix

#### 矩阵并行分配的考量标准
- 负载均衡
- 内存连续，可做Level 2,3 BLAS

#### 两个标准之间的 Trade-off

1. 1D块
    - 依次连续分配$\lceil N/P\rceil$个块给每个核。
（内存连续）
    - 缺点：不利于负载平衡，第一个连续的block算完后，process 0 空闲了

2. 1D 列循环
    - 第k列 被分到进程 $((k-1) \mod p)$
    - 解决负载平衡，但因为矩阵被切成列向量分在了不同进程，无法做Level 2,3 BLAS

![pic1](#./pic/2D1.png "1D-block")

![pic2](#./pic/2D2.png "1D-col-cyclic")

3. 1D 块循环
    - block size: NB
    - 第k列被分到进程 $((k-1)/NB \mod p)$
    - 1, 2 中的分法分别对应NB=4和NB=1的特例
    - Serial bottleneck: 列块的BLAS算法在同一个处理器上进行，而这部分其实也可以并行处理——切割行。

4. 2D 块循环
    - P个进程视为 $P_r\times P_c$的逻辑结构，它们的index为 $(p_r, p_c)$ (对应ABACUS中的coord[0], coord[1])
    - 矩阵切成MB*NB的2D块，图中MB=NB=2
    - 行对应行，列对应列

![pic3](#./pic/2D3.png "1D-block-cyclic")
![pic4](#./pic/2D4.png "2D-block-cyclic")

4对应于ABACUS中`DSIZE=1`的情况。

问题：如果`DSIZE`大于1呢?
- "进程2D块"，ABACUS中dim0, dim1是“进程2D块”的大小，根据DSIZE算出
- comm_2D $\subsetneqq$ DIAG_WORLD
- 尝试一下
## Code
### `class Parallel_Orbitals : Pdiag_Double : Pdiag_Basic`
- divide_HS_2d
    1. 进程（核）2D块相关变量
        - `DIAG_WORLD`: 全局通信子
        - `comm_2D`: 2D块通信子
        - **`DSIZE = NPOOL` //=dim0*dim1，一个进程2D块的大小, default: 1**
        - `DRANK`  default: -1, comm_2D在DIAG_WORD中的index？
        - `dim0, dim1`: the number of **processes** in each dimension   （所有CPU核被切分成dim0*dim1大小的二维矩形）
            - dim0 //小于sqrt(DSIZE)的最大约数，DSIZE是质数时，dim0=1
        - `coord[0], coord[1]`，分别是当前进程（核）在进程块comm_2D中所属的行和列坐标
    2. basis 2D块相关变量
        - `NLOCAL`: 基组大小（原子轨道总数）。H, S是NLOCAL*NLOCAL的矩阵。
        - `nb`: 1/32/64/GlobalV::NB2D   //每块包含的轨道数为nb，矩阵元数nb*nb
        - `block=NLOCAL/nb`   行/列 所分的块数（若有剩余，block+=1）
    3. basis 2D块向进程分发 相关变量
        - <u>问题：发放的范围是一个[进程2D块]？也就是说每个2D块包含了所有矩阵元的信息？（计算时就不用在不同2D块之间通信）</u>
        - `LM.row_b=block/dim0`		当前处理器被分到了多少行的basis 2D块
            - if (coord[0]<block%dim[0])  LM.row_b++
        - `nrow = LM.row_num = LM.row_b*nb`，当前进程（核）有多少行的basis（注意最后一个basis 2D块包含的轨道行数可能小于nb）
            - if(coord[0]==end_id) LM.row_num=(LM.row_b-1)*nb+(M_A-(block-1)*nb);
        - `LM.row_set[LM.row_num]` 映射：当前进程的行轨道index -> NLOCAL中的index
        - `LM.col_b ~ dim1`
        - `ncol = LM.col_num`
        - `LM.col_set[LM.col_num]`
        - `nloc=nrow*ncol`	当前进程的轨道组（矩阵元）个数
        - trace_loc_row[NLOCAL] 映射：global index -> local index (row)
        - trace_loc_col[NLOCAL] 映射：global index -> local index (col)
    - set_parameters
        - loc_size=NBANDS/DSIZE;	每个核被分到的总轨道数(余数依次加到每个核中)
        - **这个函数是否gamma_only没有区别，却写了两份完全一样的代码**
    - mpi_creat_cart
        - 【MPI_Cart_create】DIAG_WORLD分割成 dim0*dim1个小块：comm_2D
    - mat_2d
        - 【MPI_Cart_get】得到coord[0], coord[1]
        - 开始分割NLOCAL*NLOCAL的轨道矩阵（对行和列）：
            - if (dim0>block)	有的核分不到矩阵块，要设置更小的nb，切成更多blocks
            - 计算分发相关变量（LM.xxx, see 3.)
    - set_trace
     <!-- cart2blacs(comm_2D, dim0, dim1, NLOCAL, nb, nrow, desc, mpi_comm_rows, mpi_comm_cols)
        - descinit_(desc, &N, &N, &nblk, &nblk, &ISRC, &ISRC, &my_blacs_ctxt, &lld, &info);
            - N = NLOCAL，元素总数
            - nblk = nb，块总数
            - ISRC = 0，矩阵第一行/列所在进程号？
            - lld: local leading dimension -->

- diago_double_begin(ik, wfc_2d, h_mat, s_mat, ekb)
- diago_complex_begin
- readin    //not used, for hpseps?



### class `LCAO_Matrix`
- divide_HS_in_frag
    - po.divide_HS_2d
    - po.set_trace
    - allocate_HS_gamma (gamma_only)
    - allocate_HS_k (multi-k)

【待补充】
格点分法；basis 2D块index**和格点index的对应关系**
