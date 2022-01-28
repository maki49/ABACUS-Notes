## ABACUS-LCAO并行对角化前后的波函数处理

### Pdiag_Double::diago_complex_begin
(代码包括了elpa和hpseps两种KS_SOLVER，这里只整理elpa)
1. 调用`elpa_generalized_eigenvectors_dc`函数对角化，波函数输出到`wfc_2d`中，本征值输出到`wfc.ekb[ik]`中
2. 波函数2d转格点
    - MPI_Cart_rank通过进程坐标获取进程号
    - 将`wfc_2d` copy到缓冲区`work`中，然后在通信子`comm_2D`内Bcast
    - `q2WFC_WFCAUG_complex`
        - (TDDFT用的work2WFC函数和这个函数功能相同，可以合并)
        - 获取当前核上的2D矩阵元在所有行、列中的globalIndex （示意图）![glindex](./pic/glindex.jpg#pic_center "glindex")
            - naroc=[ParaO.nrow, ParaO.ncol]: 2D块分法下，当前进程有多少行/列的basis
            - iprow, ipcol: 即coord[0], coord[1], 当前进程在comm_2D中的坐标
            - （更多2D块相关变量参考[2D_Block_Cyclic.md](./2D_Block_Cyclic.md)）
        - 用`GridT.trace_lo`和`LOWF.trace_aug`获取格点分法的index
        - copy到格点分法的波函数：`WFC`和`WFC_AUG`
            - WFC[NBANDS][gt.lgd]
            - WFC_AUG[NBANDS][NLOCAL-gt.lgd]

