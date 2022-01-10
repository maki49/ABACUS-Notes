<font face="Comic Sans MS">

# Notes: LCAO-Line Wave Function ##
## Introduction
### 波函数相关变量
- LOWF.WFC_GAMMA
- LOWF.WFC_K    （和WFC_K_POOL公用内存）
- SGO.totwfc    (已弃用)
### 波函数的生成
1. 初始化
    - LCAO-line
        - in LOOP_cell, call   `wf.allocate_ekb_wg(kv.nks)`
            - 根据Ecut计算平面波的个数（用于计算v_loc和v_of_rho）
            - 根据NBANDS分配内存
        - in LOOP_elec, call  [`LOC.allocate_dm_wfc(GridT)`](#class-Local_Orbital_Charge)，分配局域轨道波函数和密度矩阵的内存
            - 不需要初猜
            - if wf.start_wfc="atomic", 什么都不做
            - if wf.start_wfc="file"    （[WF_Local](#namespace-WF_Local)）
                - read_lowf (gamma_only)
                - read_lowf_complex (multi-k)
                - distri_lowf(_new) (涉及并行分发)
                - <u>问题：会被对角化结果覆盖掉；是否只用于nscf？</u>
    - PW-line 
        - 需要初猜，因为求解KS方程用的是迭代算法。根据[wf.start_wfc](#wavefunc)选取初始化方式：
            1. atomic(default) //基于非局域赝势投影子$\left|\chi\right>$生成试探波函数
                - $\left<(log)r|\chi\right>$，径向插值表，类似双中心积分？
            2. random
            3. file (未实现，或只用于lcao-line)
2. 求解$H(k)\psi(k)=ES\psi(k)$
- Diago_LCAO_Matrix::solve_complex_matrix
    - using_LAPACK_complex
        - hm.diagoH_LAPACK(GlobalV::NLOCAL, GlobalV::NBANDS, Htmp, Stmp, GlobalV::NLOCAL, en, hvec);
        - Htmp = LM.Hloc2
        - Stmp = LM.Sloc2   //这两个是k空间的吗?确认一下
        - en: 能量本征值
        - hvec(iw, ib): 本征矢量，第ib条能带的本征波函数在第iw个轨道基上的投影系数
        - LOWF.WFC_K[ik][ib][iw]=hvec(iw, ib)
        - if NSPIN==4
            - LOWF.WFC_K[ik][ib][iw] = hvec(iw * NPOL, ib);
            - LOWF.WFC_K[ik][ib][iw + NLOCAL / NPOL] = hvec(iw * NPOL + 1, ib);
            - //NSPIN 和 NPOL?
    - using_HPSEPS_complex
3. Evolve 
- time dependent Sch.eq
- ELEC_evolve::using_LAPACK_complex    //used in TDDFT
    - Evolve_LCAO_Matrix::using_LAPACK_complex
    - ccc[j] += U_operator(j, k)*c_init[i][k]
    - LOWF.WFC_K[ik][i][j] = ccc[j]
### 参与的计算
- [测试中直接获取](#class-unkOverlap_lcao)
- 密度矩阵
    - 计算DM_gamma/[DM_k](#dmatom)
    - Exx_Abfs::DM::cal_DMk_raw
- 计算Pulay Force
    - Force_LCAO_k::set_EDM_k
- TDDFT的波函数初始化，在ELEC_scf.cpp中，需要重构
    - WFC_init=LOWF.WFC_K
- 其他
    - Chi0_hilbert::Cal_b_lcao
### 并行分配
- src_parallel: class [SubGrid_oper](#class-SubGrid_oper) (弃用)
- src_io: namespace [WF_Local](#namespace-WF_Local)  (仅在从文件读入波函数后用)
### 销毁
- Local_Orbital_wfc::~Local_Orbital_wfc()
#
[Note Tips]
- 调用其他类的函数： call xxx，链接
- 调用this->子函数：列表缩进一级    
## src_parallel
### Summarize
- common
    - parallel_common
    - parallel_global
    - parallel_reduce
    - parallel_grid     // Pgrid.init
    - parallel_kpoints   //k点并行, called in Driver::reading
- pw-line
    - parallel_pw
- lcao-line
    - parallel_orbitals
    - subgrid_oper
### `class SubGrid_oper`
#### Variables
- trace_lo_tot[NLOCAL]  <!--到轨道i为止的总占据数-->
    - used in totwfc, 对应于WFC_GAMMA中的trace_lo ?
- lgd=trace_lo_tot[NLOCAL-1]    <!--总占据数?-->
- totwfc[NSPIN][NBANDS][lgd]
#### Functions
- cal_totwfc    
    - occupy=(GridT.trace_lo>=0)
    - reduce_int_grid(occupy, NLOCAL);  //get full occupation
        - MPI_Allreduce
    - update totwfc
- dis_subwfc    分发波函数至grid group中的每个进程  **(已经弃用，被替代为`WF_local::distri_lowf和WF_local::distri_lowf_new`)**
    - bcast wf.ekb in GRID_WORLD from proc 0
    - for i in GSIZE:
        - if (GRANK=0)
            - 【 **`LOWF.WFC_GAMMA = totwfc ` was commented out**】
            - MPI_Recv lgd2, trace_lo2 from proc i
            - MPI_Send tot_wfc to proc i    (size: NBANDS*lgd2)
        - else if (GRANK=i)
            - MPI_Send lgd2, trace_lo to proc 0
            - MPI_Recv totwfc from proc 0  (size: NBANDS*GridT.lgd)
            - 【 **`LOWF.WFC_GAMMA = totwfc ` was commented out**】

## src_io
### Summarize
1. basis division: Grid & 2D
    - SGO.lgd, trace_lo2 是2D块分法的当前核轨道数?
    - GridT.lgd, trace_lo 是格点分法的当前核轨道数?
    - lgd2 被同时用作两个lgd的代称，比较迷惑
    - [SGO](#class-SubGrid_oper): totwfc弃用了，但lgd却还在用，导致这个类无法删除?
2. newdm or not
    - newdm值区分了波函数分发的两种实现方式：**是否将所有进程视为2D块结构**（cartsian location，即MPI_Cart_xx） （basis 2D块的切法没有区别）
    - newdm的优势：遍历进程2D块，然后在2D块内部Bcast（而不是一个个遍历进程，节省通信时间？）
### `namespace WF_Local`
-  `distri_lowf` (ctot, SGO.totwfc[0])
    - mu_local = SGO.trace_lo_tot
    - for i in DSIZE    //注意是DSIZE，在DIAG_WORD内通信
        - if DRANK==0   //当前进程
            - i=0: c[ib][mu_local]=ctot[ib][iw]
            - else
                - MPI_Recv trace_lo2, lgd2 from proc i
                - csend[mu_local*NBANDS+ib] = ctot[ib][iw];
                - MPI_Send ctot to proc i;    size: (NBANDS*lgd2)
        - else
            - if i==DRANK
                - MPI_Send SGO.trace_lo_tot, SGO.lgd to proc 0
                - MPI_Recv ctot from proc 0; size:(NBANDS*SGO.lgd)
    - MPI_Barrier(DIAG_WORLD)
- distri_lowf_complex
- `distri_lowf_new` (ctot, is)  (when newdm = 1)  //copy from ctot to [wfc_gamma]
    - for iprow in ParaO.dim0, ipcol in ParaO.dim1  //遍历所有[进程2D块](./2D_Block_Cyclic.md)
        - MPI_Cart_Rank(ParaO.comm_2D, [dim0, dim1], src_rank)  获取当前进程在comm_2D中的rank，然后 MPI_Bcast
        - info=CTOT2q_c  copy from ctot(global index) to work_buffer(local index)
        - MPI_Bcast work_buffer in comm_2D (在2D块内部bcast)
        - copy work_buffer to LOC.[wfc_dm_2d](#class-WFC_DM_2d).wfc_k
- distri_lowf_complex_new
- read_lowf //从文件读入平面波
    - LOWF_K.dat
    - LOWF_GAMMA_S[is].dat
- read_lowf_complex
- write_lowf

### `class unkOverlap_lcao`
- 这个类似乎是为了测试 center_2_orb 而写的
- get_lcao_wfc_global_ik(lcao_wfc_global[ik],GlobalC::LOWF.WFC_K[ik]) //获取每个cpu核的原子轨道系数
    - 和WF_Local中的功能类似



## `src_lcao`
### `class Local_Orbital_wfc`
#### Varaibles
- WFC_K[NK, NBANDS, NLOCAL]
    - how to construct dm?
- WFC_K_POOL [NK* NBANDS* NLOCAL]
- WFC_GAMMA_aug [NSPIN, NBANDS, daug]
- ORB_control orb_con

#### Functions
- allocate_k    //allocate WFC_K and WFC_K_POOL （二者共用内存）
    - if(wf.start_wfc=="file") call [WF_Local](#namespace-WF_Local)::read_lowf_complex
- set_trace_aug //not used

### `class WFC_DM_2d`
[2D Block-Cyclic](./2D_Block_Cyclic.md) Density Matrix 
#### Variables
- wfc_gamma[is](ib, iw)
- wfc_k[ik](ib, iw)     ParaO.ncol * ParaO.nrow     //<u>为什么行数ncol, 列数是nrow?</u>
- dm_gamma[is](iw1, iw2)    
- dm_k[ik](iw1, iw2)
#### Functions
- init
- cal_dm([wf.wg](#wg-ref)   //called in [Local_Orbital_Charge](#class-Local_Orbital_Charge)::sum_bands
    - for ib_global < wg.nc
        - ib_local=ParaO.trace_loc_col[ib_global]
        - if(ib_global>=0)  wg_local=wg(is,ib_global)
    - for (ir<wg_wfc.nr) LapackConnector::scal

### `class Local_Orbital_Charge`
#### Variables
- DM_pool[NSPIN][lgd*lgd]   为了让DM内存连续?
- DM[NSPIN][lgd][lgd]
- DM_R[NSPIN, [LNNR](#lnnr).nnrg]
- sender_2D_index  //DM矩阵元的buffer-index(i)到2D-block-index(idx)的映射
- sender_size
- sender_size_process[nprocs]
- receiver_2D_index
- receiver_size
- receiver_size_process[nprocs]
- receiver_local_index[receiver_size]

- lgd_now   //sub-FFT-mesh orbitals number in this step
- lgd_last   //sub-FFT-mesh orbitals number in previous step
- nnrg_now   //sub-FFT-mesh orbitals number in this step, with k
- nnrg_last   //sub-FFT-mesh orbitals number in previous step, with k
#### Functions
- allocate_dm_wfc(GridT)
    - allocate_gamma    //(actually allocate_DM_gamma?)
        - DM[is][i] = &DM_pool[is][i*lgd];
        - **`setAlltoallvParameter`**(ParaO.comm_2D, ParaO.blacs_ctxt, ParaO.nb);
            - receiver_size=lgd_now*lgd_now
            - nblk=ParaO.nb
            - trace_2D_row[iLocalGrid]=`localIndex`(iGlobal, nblk, nprows, p)
                - // 局域格点到局域2D块的映射?
                - p=int((iGlobal%(nblk*nprows))/nblk);
                - return int(globalIndex/(nblk*nprocs))*nblk+globalIndex%nblk;
    - call [LOWF](#class-Local_Orbital_wfc).allocate_k
    - allocate_DM_k
        - LNNR.nnrg //?  每个核上的实空间格点矩阵元个数?
        - wfc_dm_2d.init()
        - kpt_file(GridT)
- kpt_file  读入波函数文件
    - [WF_Local](#namespace-WF_Local)::read_lowf_complex  读入“LOW_K_(ik+1).dat”文件
- sum_bands //计算电子密度（同时也计算能量?）
    - en.eband += wf.ekb[ik][ib] * wf.wg(ik, ib)
    
    //if gamma_only: 
    - call [wfc_dm_2d](#class-WFC_DM_2d).cal_dm
    - cal_dk_gamma_from_2D()    //**2D block-cyclic DM -> grid distributed DM**
        1. **send**
        - idx = sender_2D_index[i]
            - $idx=icol\times NLOCAL + irow$
        - icol = idx%NLOCAL
        - irow = (idx-icol)/NLOCAL
        - sender_buffer[i]=wfc_dm_2d.dm_gamma[is](icol, irow); // sender_buffer is clomun major
        ![2D1](./2D-1.jpg#pic_center "2D-1" )
        <!-- <img src=./2D-1.jpg#pic_center width=60% /> -->
        - `MPI_Alltoallv`(sender_buffer, sender_size_process, sender_displacement_process, MPI_DOUBLE, receiver_buffer, receiver_size_process, receiver_displacement_process, MPI_DOUBLE, GlobalC::ParaO.comm_2D);
            - send_buffer: 缓冲区（起始地址）
            - sender_size_process[i]: 往第i个进程发多少个sendtype的数据
            - sender_displacement_process[i]: 被发往第i个进程的数据相对缓冲区起始地址的位移量
            - MPI_DOUBLE: 第四个参数, send_type
            - ParaO.comm_2D: 通信子，MPI_Comm类的实例
        2. **receive**
        - idx = receiver_local_index[i]
        - icol = idx%lgd_now
        - irow = (idx-icol)/lgd_now
        - DM[is][irow][icol] = receiver_buffer[i]
            - //GridT.trace_lo 给出了全局 grid index (max: NLOCAL) 到局域 grid index (max: lgd_now) 的映射
    - cal_dk_gamma()    //used in hpseps
    

    //else (multi-k)
    - cal_dk_k(GridT)
        - //**只用了GridT.in_this_processor 一个变量**
        - Record adj RA;    RA.for_grid(GridT)     //记录轨道沿z方向的格点切分之后，每个核所分到格点的近邻原子的信息?
        - for each atoms:
            - gstart = LNNR.nlocstartg[iat]
            - ng = [LNNR](#lnnr).nlocdimg[iat]
            - <span id="dmatom">DM_ATOM[NSPIN][ng]</span>
            - call cal_DM_ATOM  //inline in DM_k.cpp
                - nw1: 当前原子(ia1)的轨道数
                - WFC_PHASE[NBANDS*nw1] = 
                - wfc = [LOWF](#class-Local_Orbital_wfc).WFC_K[ik]
                - for ik in kv.nks  //当前pool（处理器?）的k点个数??
                    - for each adjacent atom ia2:
                        <!-- WFC_PHASE[iline+iw1] = $e^{2\pi i \mathbf{k}_{direct}\cdot \mathbf{r_{12}}}$ * wf.wg(ik, ib) * conj(wfc[ib][iw1_lo+iw1]) -->
                        - zgemm_: 计算 $\rho^{ia1, ia2}_{\mu\nu} = \sum_{ia2}^{grid\_adj(ia1)}\sum_{ib}c^{ia1*}_{\mu, ib}c^{ia2}_{ib, \nu }\cdot f_{ik, ib}e^{2\pi i \mathbf{k}_{direct}\cdot \mathbf{r_{12}}}$ (DM_ATOM中的[nw1*nw2]块)

    - call [wfc_dm_2d](#class-WFC_DM_2d).cal_dm



## src_pw
### `class wavefunc: WF_atomic : WF_igk`
#### Variables
- (WF_igk)
    - npwx: 最大平面波数
    - npw
    - igk[nks, npw_max]
    - g2kin[npwx], 当前k点的动能
    - ggpsi: planewave cut off for wavefunctions, unit(2pi/a0)
    - ggwfc=4*ggpsi
    - ggwfc2=4*ggpsi??
- (WF_atomic)
    - tabel_local   //one-dimension tables for NAOs
    - **`evc`** //wavefunctions in PW basis
    - wanf2
- (wavefunc)
    - ekb[nks][NBANDS]?
    - start_wfc //string, [波函数生成方式](#波函数的生成)
    - <span id="wg-ref">[wg(ik, ib)]</span>每个k点和能带的weight
#### Functions

- (WF_igk)
    - get_sk
    - setupIndGk
- (WF_atomic)
    - atomic_wfc
    - random
    - atomicrandom
    - check_psi
    - init_at_1 //径向赝波函数的插值表?
- <span id="wavefunc">(wavefunc)</span>
    - (PW-line) allocate
    - **(LCAO-line) allocate_ekb_wg(nks)**      //called in LOOP_cell::opt_cell **（为什么LCAO-line要用这个函数？因为要计算v_loc和v_of_rho？）
        - npwx=setupIndGk(pw, nks)**
            - 计算kv.ngk[ik]    //每个k点的平面波个数
                - 算了两次！"not a smart job "
                    - 第一次：get npw_max, allocate `igk`
                    - 第二次：for each igk, set ig
            - for ig < pwb.ngmw 
                - ngmw  // num of G vectors within ggwfc2 in each proc
                - if (sqrt(pwb.gg[ig]) > sqrt(pwb.ggpsi) + sqrt(k2) &&  (gk2 <= pwb.ggpsi)), ng++
                    - pwb.gg    // $G^2$
                    - pwb.ggpsi     //平面波 $G^2/\frac{2\pi}{a_0}$ cutoff
                    - gk2   //$|G+k|^2$
        - ekb[nks][NBANDS]  //band energies
        - wg(nks, NBANDS)      //weight of each k-point and band
    - (PW-line / lcao_in_pw) wfcinit
        - wfcinit_k
            - (PW-line) diago_PAO_in_pw_k(ik, wf.evc[ik])
                - hm.hpw.init_k(ik)
                    - get_starting_nw
                    - file: not implemented?
                    - atomic: ucell.natomwfc*(atomic) + (NBANDS-ucell.natomwfc)*(random)
                    - random: all NBANDS random
                1.  atomic_wfc
                    - gk[ig] = $\frac{2\pi}{a_0}$ |gcar[ig]+kvec_c[ik]|
                    - get_sk: 相因子$e^{-i(\mathbf{G}+\mathbf{k})\cdot\tau}$
                    - flq = radial fourier transform of atomic orbitals chi
                    - wfcatom(index, ig) = $i^l$ * sk [ig] * ylm(lm, ig) * flq[ig];
                2. atomicrandom

                3. random
    - (LCAO_in_PW) LCAO_in_pw_k; LCAO_in_pw_k_q

    - Chi0_hilbert ??? (about spectrum?)
    - prepare_k
    
## Others
-  `lcao_in_pw`: 局域轨道用平面波展开（ Q: 仅在lcao_in_pw时用，还是lcao-line都要用?）
    - LCAO_in_pw_k
- GSIZE=DSIZE, 格点2D块?
- ELEC_scf.cpp line #180    TDDFT的波函数初始化应该移走?
    - WFC_init = LOWF.WFC_K[ik][ib][i]
- GridT.trace_lo    轨道index，全局(NLOCAL)到局域(lgd_now)的映射
- <span id="lnnr">LNNR</span>
    - nnrg  //当前处理器上，实空间格点的矩阵元个数 (g means "grid")
    - nlocdimg
    - nlocstartg
    - nnr   //当前处理器上，2D块overlap矩阵元的个数
    - nlocdim[iat]  //某个特定原子的nnr
    - nlocstart[iat] //某个特定原子在当前处理器上的2D块的第一个矩阵元在nnr
